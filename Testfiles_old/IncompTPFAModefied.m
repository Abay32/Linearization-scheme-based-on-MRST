function [state] = IncompTPFAModefied(state, G, rock, fluid,  varargin)
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
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   % Preliminaries: set various constants and maps
   nf     = G.faces.num;%number of faces 
   nc     = G.cells.num;%number of cells
   cf     = G.cells.faces(:,1);% faceses of a given cells
   nhf    = numel(cf);%number of repeated faces of cells.(faces willbe repeated exepts the boundary faces)
   hf2cn  = gridCellNo(G);
   
   iface  = all(G.faces.neighbors ~= 0, 2);
   eface  = ~iface;
   ni     = G.faces.neighbors(iface,:);
   
 
   % Define effective face transmissibility as harmonic average of
   % viscosity-weighted one-sided transmissibilities.
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s,state);
   
   
   
   if isfield(fluid, 'pc'),
       pc = fluid.pc(state);
       po = state.pressure + pc;
   else 
       po = state.pressure;
   end
   %upsteaming the internal faces mobilities,
   %krUp = kr(iface,:);

   mob    = bsxfun(@rdivide, kr, mu);
   %avmob  = (mob(ni(:,1),1) + mob(ni(:,2),1))/2;%arthimetic averaged mobilities at internal faces
   totmob = sum(mob,2);
   
   
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;  
   halfTran = simpleComputeTransModefiedold(G, rock,totmob);
 
   assert(numel(halfTran.hTM) == numel(hf2cn), ...
      ['Expected one one-sided transmissibility for each ', ...
      'half face (=%d), but got %d.'], numel(hf2cn), numel(halfTran.hTM));
   TM = halfTran.hT.*totmob(hf2cn);  
   T  = 1 ./ accumarray(cf, 1 ./ TM, [G.faces.num, 1]); % full transmisibility with mobility 
   TP = 1 ./ accumarray(cf, 1 ./ halfTran.hT, [G.faces.num, 1]); % full transmisibility without mobility 
  % TP = TP./2;
   % Compute gravity contribution to right-hand side
   
   cvec = G.faces.centroids(cf, :) - G.cells.centroids(hf2cn, :);   
   gp   = rho(1) .* (cvec * gvec.');
   
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
      ind  = Dface(cf);      
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
   
   dpo = po(ni(:,2)) - po(ni(:,1)) - faceGrav(iface);
   upco = (double(dpo)<=0);   
   mup = Simpleupstreamold(G,upco, mob(:,2));%upwided oil mobility
   
   pw  = state.pressure;
   dpw = pw(ni(:,2)) - pw(ni(:,1)) - faceGrav(iface);
   upcw = (double(dpw)<=0); 
   wmup = Simpleupstream(G,upcw, mob(:,1));
   Uptotmob = mup + wmup;   
   %% Right hand side
   % no problem on T as long as we have homogeneous <DBC>
   rhs = accumarray(hf2cn, - T(cf).*(sgn.*faceGrav(cf) + rhs1), [nc, 1]) + ...
      rhs2 + accumarray(hf2cn, -rhs3(cf), [nc, 1]);

   clear rhs1 rhs2 sgn;   
   %% Capillar pressure
   if isfield(fluid, 'pc'),
       sgn         = 2*(G.faces.neighbors(bc.face,2)==0)-1;
       cNo         = sum(G.faces.neighbors(bc.face,:), 2); % cell number corrosponding to outer faces
       pc  = fluid.pc(state);
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
           CapPressure  = accumarray([ic1; ic2], [PcFace(internal); -PcFace(internal)], size(po));          
           %CapPressure(cNo) = CapPressure(cNo) + sgn.*mob(cNo,2).*TP(bc.face).*(pc(bc.sat(:,1)) - pc(cNo));
           clear ic1 ic2 cNo sgn;
           rhs = rhs + CapPressure;
       end
   end   
   %% Add up internal face transmissibilities plus Dirichlet pressure faces for each cell.
   n   = size(ni,1);
   C   = sparse( [(1:n)'; (1:n)'], ni, ones(n,1)*[1 -1], n, G.cells.num);
   grad = -C;
   div  = C';
   Flux = -bsxfun(@times, Uptotmob.*TP(iface), grad);
   A = div*Flux;
   
   if any(Dface)       
        cNo = sum(G.faces.neighbors(bc.face(is_press),:), 2);
        for i = 1:length(cNo)
            A(cNo(i),cNo(i)) =  bsxfun(@plus,A(cNo(i),cNo(i)), totmob(cNo).*TP(pface));
        end
   end
  
   clear I J V;

   % If there are no Dirichlet boundary conditions, do a fix to ensure that
   % we have a solvable system
%    if ~any(Dface)
%       if A(1) > 0,
%          A(1) = 2*A(1);
%       else
%          [j, j] = max(diag(A));  %#ok
%          A(j,j) = 2*A(j,j);
%       end
%    end   
   % Solve the flow problem
   p = opt.LinSolve(A, rhs);
   clear A rhs;

%    % Reconstruct face pressures and fluxes.
%    fpress =  ...
%       accumarray(cf, (p(hf2cn) + gp).*halfTran.hT, [G.faces.num,1])./ ...
%       accumarray(cf, halfTran.hT, [G.faces.num,1]);
% 
%    % Recompute face pressure at Neumann faces
%    
%    b         = any(G.faces.neighbors==0, 2);
%    fpress(b) = fpress(b) - rhs3(b)./T(b);
%    % Reset correct values at Dirichlet faces
%    fpress(Dface) = DfaceVal;

   % Compute face fluxes for internal faces
   %ni   = G.faces.neighbors(iface,:);
   
   pw  = p;      
   dpw = pw(ni(:,2)) - pw(ni(:,1)) - faceGrav(iface);  
   upcw = (double(dpw)<=0);   
   mup  = Simpleupstream(G,upcw, mob(:,1));   
   
   flux = -accumarray(find(iface),...
      mup.*TP(iface).*(pw(ni(:,2)) - pw(ni(:,1)) - faceGrav(iface)), [nf, 1]);

%    %Compute fluxes for external faces using Darcy's law   
%    sgn         = 2*(G.faces.neighbors(eface,2)==0)-1;
%    cNo         = sum(G.faces.neighbors(eface,:), 2); % cell number corrosponding to outer faces
%    
%    faceGrav    = accumarray(cf, gp, [nf, 1]);
%    flux(eface) = -sgn.*mob(cNo,1).*TP(eface).*(fpress(eface) - pw(cNo) - faceGrav(eface));
   
    
   clear sgn;
   %% ---------------------------------------------------------------------
   %adding up the boundary fluxes to the system 
   sgn         = 2*(G.faces.neighbors(bc.face,2)==0)-1;
   
   HT = TP(bc.face);
   
   bcFace = G.faces.neighbors(bc.face,:) ;%boundary faces neighbors   
   assert(~any(all(bcFace > 0, 2)),'bc on internal boundary');
   ss = sum(bc.sat, 2);
   % Values should either sum to zero or one
   assert(all(ss - 1 < sqrt(eps) | ss < sqrt(eps)));
   
   BCcells = sum(bcFace, 2);
   nbc = numel(bc.face);
   cellToBCMap = sparse((1:nbc)', BCcells, 1, nbc, G.cells.num);
   %BCTocellMap = cellToBCMap';
   qBcw = zeros(nbc,1);
   %extract boundary cell pressures and other preperties
   pBc   = cellToBCMap*pw;
   sBc   = cellToBCMap*state.s;
   mobBc = cellToBCMap*mob ;
   %rhoBc = rho(BCcells,:);
   
   sat   = bc.sat;
   NoSat = all(sat == 0, 2);
   hasNoSat = any(NoSat);
   %   
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
   
   isP = reshape(is_press,[],1);   
   dP = bc.value(isP) - pBc(isP) + rho(1).*dzbc(isP);   
   
   % Determine if pressure bc are injecting or producing
   injDir = dP > 0;
   
   injP = isP;
   injP(isP) = injDir;
   if any(~injDir)
      % Write out the flux equation over the interface
      subs = isP & ~injP;      
      qBcw(subs) = -sgn(subs).*mobBc(subs,1).*HT(subs).*dP(~injDir);
      clear subs
   end
   %
   mobBcF  = fluid.relperm(bc.sat,state);
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
   
   state.pressure(1:nc) = p(1:nc);
   state.flux(:)        = flux;
   
end