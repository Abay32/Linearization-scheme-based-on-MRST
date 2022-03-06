function [state] = simpleIncompTPFAModefiedold(state, G, rock, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   state = simpleIncompTPFA(state, G, hT, fluid)
%   state = simpleIncompTPFA(state, G, hT, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (sources,
%   and boundary conditions).
%
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir solution structure either properly
%            initialized from function 'initResSol'
%
%   G, hT  - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
%   LinSolve - Handle to linear system solver software to which the
%            fully assembled system of linear equations will be passed.
%            Assumed to support the syntax
%
%              x = LinSolve(A, b)
%
%            in order to solve a system Ax=b of linear equations.
%            Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   state - Update reservoir solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'bc' and 'src' are empty and there are no effects of gravity, then the
%   input value 'xr' is returned unchanged and a warning is printed in the
%   command window. This warning is printed with message ID
%
%           'incompTPFA:DrivingForce:Missing'
%
%
% SEE ALSO:
%   computeTrans, addBC, addSource, addWell, initSingleFluid, initResSol,
%   initWellSol.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   opt = struct('bc', [], 'src', [], ...
                'LinSolve',     @mldivide, ...
                'pc_form'     , 'wetting');
            
   opt = merge_options(opt, varargin{:});
   % Sanity checks on grid and gravity term
   assert (1 <= G.griddim && G.griddim < 4);
   gvec = gravity();
   gvec = gvec(1:G.griddim);
   assert (all(size(gvec) == [1,G.griddim]));
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ~(norm(gvec) > 0)]),
       if ~isfield(fluid,'pc')
           warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
       end
   end

   % Preliminaries: set various constants and maps
   nf     = G.faces.num;%number of faces 
   nc     = G.cells.num;%number of cells
   cf     = G.cells.faces(:,1);% faceses of a given cells
   nhf    = numel(cf);%number of repeated faces of cells.(faces willbe repeated excepts the boundary faces)
   hf2cn  = gridCellNo(G);
   
   iface  = all(G.faces.neighbors ~= 0, 2);
   eface  = ~iface;
   ni     = G.faces.neighbors(iface,:)   ;
 
   % Define effective face transmissibility as harmonic average of
   % viscosity-weighted one-sided transmissibilities.
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s,state);
   % 
  
   if isfield(fluid, 'pc'),
       pc = fluid.pc(state);
       po = state.pressure + pc;
   else 
       po = state.pressure;
   end
   %upsteaming the internal faces mobilities,
   %krUp = kr(iface,:);

   mob    = bsxfun(@rdivide, kr, mu);
   avmob  = (mob(ni(:,1),1) + mob(ni(:,2),1))/2;%arthimetic averaged mobilities at internal faces
   totmob = sum(mob,2);   
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;
   halfTran = simpleComputeTransModefiedold(G, rock,totmob) ;
   
   assert(numel(halfTran.hTM) == numel(hf2cn), ...
      ['Expected one one-sided transmissibility for each ', ...
      'half face (=%d), but got %d.'], numel(hf2cn), numel(halfTran.hTM));
   TM = halfTran.hT.*totmob(hf2cn);
   T  = 1 ./ accumarray(cf, 1 ./TM, [G.faces.num, 1]); % full transmisibility with mobility 
   TP = 1 ./ accumarray(cf, 1 ./halfTran.hT, [G.faces.num, 1]); % full transmisibility without mobility 
   
   % Compute gravity contribution to right-hand side   
   cvec = G.faces.centroids(cf, :) - G.cells.centroids(hf2cn, :);   
   gp   = rho(1) .* (cvec * gvec.')  ;
   clear cvec;
   
   % Initiatlize right-hand side
   rhs1 = zeros(nhf, 1);
   rhs2 = zeros(nc,  1);
   rhs3 = zeros(nf,  1);
   
   % Source terms
   src = opt.src;
   if ~isempty(src),
      % Compatibility check on cell numbers for source terms
      assert (max(src.cell) <= nc && min(src.cell>0), ...
         'Source terms refer to cell not existant in grid.');
      % Sum source terms inside each cell and add to rhs
      s  = accumarray(src.cell, src.rate);
      ii = accumarray(src.cell, 1)> 0;
      
      rhs2(ii) = rhs2(ii) + s(ii);
   end
   % Boundary conditions
   Dface    = false([nf, 1]);
   DfaceVal = [];
   bc       = opt.bc;
   
   %
   if ~isempty(bc),
       % Compatibility checks
      assert (max(bc.face) <= nf && min(bc.face) > 0, ...
         'Boundary condition refer to face not existant in grid.');
      assert (all(accumarray(bc.face, 1, [nf, 1]) <= 1), ...
         'There are repeated faces in boundary condition.');
      
      % Pressure (Dirichlet) boundary conditions.
      % Extract the faces marked as defining pressure conditions. Define a
      % local numbering (map) of the face indices to the pressure condition
      % values.
      is_press = strcmpi('pressure', bc.type);
      pface    = bc.face(is_press);
      DfaceVal = bc.value(is_press);
      map      = sparse(double(pface), 1, 1:numel(pface));
     
      % Mark the faces as having pressure boundary conditions.  This
      % information will be used to eliminate known pressures from the
      % resulting system of linear equations.
      Dface(pface) = true;

      % Enter Dirichlet conditions into system right hand side. Relies
      % implictly on boundary faces being mentioned exactly once in
      % G.cells.faces(:,1).
      ind       = Dface(cf);      
      
      rhs1(ind) = -DfaceVal(map(G.cells.faces(ind,1)));
      %
      clear ind
      % Reorder Dirichlet conditions according to SORT(pface) so that we
      % later can set 'X(Dface) = DfaceVal' even when ISLOGICAL(Dface).
      DfaceVal = DfaceVal(map(Dface));

      % Flux (Neumann) boundary conditions.
      % Note negative sign due to bc.value representing INJECTION flux.
      is_flux = strcmpi('flux', bc.type);  
      rhs3(bc.face(is_flux)) = -bc.value(is_flux);
   end
   assert(~any(DfaceVal < 0), 'Pressure conditions should always be non-negative');
   
   % Add gravity contribution to all internal faces and faces with
   % Dirichlet boundary conditions  
 
   sgn = 2*(G.faces.neighbors(cf, 1) == hf2cn) - 1;
   j   = iface(cf) | Dface(cf);  
   
  
   faceGrav  = accumarray(cf(j), gp(j).*sgn(j), [nf, 1]); 
   %% Upwinding   
   dpo  = po(ni(:,2)) - po(ni(:,1)) - faceGrav(iface);
   upco = (double(dpo)<=0);   
   mup  = Simpleupstream(G,upco, mob(:,2));%#OK., upwided oil mobility   
   pw   = state.pressure;
   dpw  = pw(ni(:,2)) - pw(ni(:,1)) - faceGrav(iface);
   upcw = (double(dpw)<=0); 
   
   wmup = Simpleupstream(G,upcw, mob(:,1));
   Uptotmob = mup + wmup;   
   %% Right hand side
   % no problem on T as long as we have homogeneous <DBC>
   rhs = accumarray(hf2cn, -T(cf).*(sgn.*faceGrav(cf) + rhs1), [nc, 1]) + ...
      rhs2 + accumarray(hf2cn, -rhs3(cf), [nc, 1]);
   
   clear rhs1 rhs2 sgn;   
   %% Capillar pressure
   if isfield(fluid, 'pc'),
%        ijk      = gridLogicalIndices(G);
%        upInd    = (ijk{1} > nc/2);
%        pc = zeros(nc,1);
%        if any(~upInd)
%            pc(~upInd) = fluid.pcl(state.s(~upInd,1))
%        end
           
         
       

       po  = state.pressure + pc;
       dpo = po(ni(:,2)) - po(ni(:,1)) - faceGrav(iface);  
       
       upco = (double(dpo)<=0);   
       mup = Simpleupstream(G,upco, mob(:,2));       
       if any(abs(pc) > 0),
           PcFace =  accumarray(find(iface),...
               mup.*TP(iface).*(pc(ni(:,2)) - pc(ni(:,1))), [nf, 1]);
           %PcFace(eface) = mob(cNo,2).*pc(cNo);
           neighbors = getNeighbourship(G, 'Topological', true);
           internal = all(neighbors~=0, 2);
           ic1  = neighbors(internal,1);
           ic2  = neighbors(internal,2);           
           CapPressure  = accumarray([ic1; ic2], [PcFace(internal); -PcFace(internal)], size(po))     ;     
           %CapPressure(cNo) = CapPressure(cNo) + sgn.*mob(cNo,2).*TP(bc.face).*( - pc(cNo));
           clear mup ic1 ic2 cNo sgn;
           rhs = rhs + CapPressure;
       end
   end   
   %% Add up internal face transmissibilities plus Dirichlet pressure faces for each cell.
   reshape(G.faces.neighbors(iface,:), [], 1)
   d  = accumarray(reshape(G.faces.neighbors(iface,:), [], 1), repmat(Uptotmob.*TP(iface), [2,1]),  [nc, 1]) 
   % Assemble coefficient matrix for internal faces.  Boundary conditions
   % may introduce additional diagonal entries. 
   
   I  = [G.faces.neighbors(iface,1); G.faces.neighbors(iface,2); (1:nc)'];
   J  = [G.faces.neighbors(iface,2); G.faces.neighbors(iface,1); (1:nc)'];  
   V  = [-Uptotmob.*TP(iface); -Uptotmob.*TP(iface); d]; 
   
   if any(Dface)
       cNo = sum(G.faces.neighbors(bc.face(is_press),:), 2);
       %   .*halfTran.hT(Dface(cf)    
       B = sparse(cNo,cNo,totmob(cNo).*TP(pface),nc,nc);      
       A  = sparse(double(I), double(J), V, nc, nc); %accumulate V in I,J indecies
       A = A+B;
   elseif ~any(Dface)        
        
       A  = sparse(double(I), double(J), V, nc, nc); %accumulate V in I,J indecies
             
       % If there are no Dirichlet boundary conditions, do a fix to ensure that
       % we have a solvable system
       if A(1) > 0,
           A(1) = 2*A(1);
       else
         [j, j] = max(diag(A));  %#ok
         A(j,j) = 2*A(j,j);
       end
   end
   
   clear I J V d;
   % Solve the flow problem
   n     = size(ni,1) ;% 
   C     = sparse([(1:n)'; (1:n)'], ni, ones(n,1)*[1 -1], n, G.cells.num);
   grad  = -C;
   div   =  C'; 
   pflux = -bsxfun(@times, Uptotmob.*TP(iface), grad);
   Ap    = div*pflux;
   %
   if any(Dface)        
       cNo  = sum(G.faces.neighbors(bc.face(is_press),:), 2);
       pB   = sparse(cNo, cNo,   0.5*totmob(cNo).*TP(pface), nc, nc);
        
       Ap   = Ap   + pB;%pressure coefficient matrix at the pressure euation 
       %Ap_o = Ap_o + paxB; % a materix from the phase oil transport equation
   end
   if Ap(1) > 0,
       Ap(1) = 2*Ap(1);
   else
       [j, j] = max(diag(Ap));  %#ok
       Ap(j,j) = 2*Ap(j,j);
   end
    
   p = opt.LinSolve(A, rhs);
   clear A rhs;

   % Reconstruct face pressures and fluxes.
   fpress =  ...
      accumarray(cf, (p(hf2cn) + gp).*TM, [G.faces.num,1])./ ...
      accumarray(cf, TM, [G.faces.num,1]);

   % Recompute face pressure at Neumann faces
   
   b         = any(G.faces.neighbors==0, 2);
   fpress(b) = fpress(b) - rhs3(b)./T(b);
   % Reset correct values at Dirichlet faces
   fpress(Dface) = DfaceVal;

   % Compute face fluxes for internal faces
   %ni   = G.faces.neighbors(iface,:);
   
   pw   = p;      
   dpw  = pw(ni(:,2)) - pw(ni(:,1)) - faceGrav(iface);  
   
   upcw = (double(dpw)<=0);   
   mup  = Simpleupstream(G,upcw, mob(:,1));
   
   
   flux = -accumarray(find(iface),...
      mup.*TP(iface).*(pw(ni(:,2)) - pw(ni(:,1)) - faceGrav(iface)), [nf, 1]);

   %Compute fluxes for external faces using Darcy's law   
   sgn         = 2*(G.faces.neighbors(eface,2)==0)-1;
   cNo         = sum(G.faces.neighbors(eface,:), 2); % cell number corrosponding to outer faces
   
   faceGrav    = accumarray(cf, gp, [nf, 1]);
   flux(eface) = -sgn.*mob(cNo,1).*TP(eface).*(fpress(eface) - pw(cNo) - faceGrav(eface));
   
    
   clear sgn;
   %% ---------------------------------------------------------------------
   % adding up the boundary fluxes to the system 
   if ~isempty(bc)
       sgn = 2*(G.faces.neighbors(bc.face,2)==0)-1   ;
       HT  = TP(bc.face);
       
       bcFace = G.faces.neighbors(bc.face,:); %boundary faces neighbors   
       assert(~any(all(bcFace > 0, 2)),'bc on internal boundary');
       ss = sum(bc.sat, 2);
       % Values should either sum to zero or one
       assert(all(ss - 1 < sqrt(eps) | ss < sqrt(eps)));
       
       BCcells = sum(bcFace, 2);%Cells related to the boundary faces
       nbc = numel(bc.face);
       
       cellToBCMap = sparse((1:nbc)', BCcells, 1, nbc, G.cells.num);
       BCTocellMap = cellToBCMap';
       qBcw = zeros(nbc,1);
       %extract boundary cell pressures and other preperties
       pBc   = cellToBCMap*pw;
       sBc   = cellToBCMap*state.s;
       mobBc = cellToBCMap*mob ;
       %rhoBc = rho(BCcells,:);
       
       sat   = bc.sat;
       NoSat = all(sat == 0, 2);
       hasNoSat = any(NoSat);
              
       if any(strcmpi(G.type, 'topSurfaceGrid'))
           dzbc = model.gravity(3) * (G.cells.z(BCcells) - G.faces.z(bc.face));
       else
           g = gravity();
           g = g(1:G.griddim);
           dz = G.cells.centroids(BCcells, :) - G.faces.centroids(bc.face,:);
           dzbc = dz*g';
       end
       % If no saturations are defined, we explicitly set it to mirror the
       % cell values on the other side of the interface
       if hasNoSat
           sBC = double(sBc);
           sat(NoSat, :) = sBC(NoSat,:);
       end
       isP = reshape(is_press,[],1)   ;
       dP = bc.value(isP) - pBc(isP) + rho(1).*dzbc(isP);   
       
       % Determine if pressure bc are injecting or producing
       injDir = dP > 0; %if satisfied, the flow is from the boundary to the cell else from the cell to the boundary 
       
       injP = isP;
       injP(isP) = injDir;
       if any(~injDir)
           % Write out the flux equation over the interface
           subs = isP & ~injP;      
           qBcw(subs) = -sgn(subs).*mobBc(subs,1).*HT(subs).*dP(~injDir);
           clear subs
       end
       %
       mobBcF  = bsxfun(@rdivide, fluid.relperm(bc.sat ,state), mu);
       totMob = sum(double(mobBcF),2);
       if any(injDir)
           % In this case, pressure drives flow inwards, we get the injection rate
           % determined by the sat field            
           subs = isP & injP;      
           
           qBcw(subs) = -sgn(subs).*totMob(subs).*HT(subs).*dP(injDir);   
           
           clear subs
       end
       %-----------------------------------------------------------------------
       % Treat flux / Neumann BC
       injNeu = bc.value > 0;     
       subs = ~isP &  injNeu;
       if any(subs)
           % Injection
           %q_s(subs) = bc.value(subs).*sat(subs, i);
           qBcw(subs) = bc.value(subs);%.*sat(subs, 1);%./bBC(subs);
       end
       subs = ~isP & ~injNeu;
       if any(subs)
           % Production fluxes, use fractional flow of total mobility to
           % estimate how much mass will be removed.
           f   = mobBcF(subs)./totMob(subs);
           tmp = f.*bc.value(subs);
           
           %q_s(subs) = tmp;
           qBcw(subs) = tmp;%./bBC(subs);
       end
       flux(bc.face) = qBcw;
       % Extract the corresponding state variables   
       state.pressure  = p ;
       state.flux(:)   = flux;
   end
end
 