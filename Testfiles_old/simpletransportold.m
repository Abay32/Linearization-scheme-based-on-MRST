function state = simpletransportold(G,state,rock,dT,varargin)
% opt = struct('verbose', mrstVerbose, 'gravity', gravity(), ...
%                 'wells', [], 'src', [], 'bc', [], 'Trans', [],'dhfz',[]);
% opt = merge_options(opt, varargin{:});
% 
% assert ((size(state.s,2) < 3) || all(state.s(:,3) == 0), ...
%            'Function ''%s'' is for two-phase flow only', mfilename);
% 
% assert(all(isfinite(state.flux)), 'Passed state contained non-finite fluxes');
% % All source terms, i.e., boundary conditions, source terms and wells.
% 
% compi = {'use_compi', true};
% q = computeTransportSourceTerm(state, G, opt.wells, ...
%                                   opt.src, opt.bc, compi{:});
% q = assembleTransportSource(q, G.cells.num, compi{:});
% assert(all(isfinite(q)))
opt = struct(...
      'verbose'  , false      , ...
      'onlygrav' , false      , ...
      'computedt', true       , ...
      'max_dt'   , inf        , ...
      'dt_factor', 0.5        , ...
      'wells'    , []         , ...
      'src'      , []         , ...
      'bc'       , []         , ...
      'dt'       , tf         , ...
      'Trans'    ,[]          , ...
      'gravity'  , gravity()  , ...
      'satwarn'  , sqrt(eps));

opt = merge_options(opt, varargin{:});

assert ((size(state.s,2) < 3) || all(state.s(:,3) == 0), ...
    'Function ''%s'' is for two-phase flow only', mfilename);
       
assert(all(isfinite(state.flux)), 'Passed state contained non-finite fluxes');
% All source terms, i.e., boundary conditions, source terms and wells.
compi = { 'use_compi', true };

q = addfluxfromSourceandBC(state, G, opt.src, opt.bc, compi{:});
q = assemblefluxfromSourceandBC(q, G.cells.num, compi{:});
assert(all(isfinite(q)))

source = max(q,0) + min(q,0);

iface  = all(G.faces.neighbors ~= 0, 2);
eface  = ~iface;
ni     = G.faces.neighbors(iface,:)   ;

n   = size(ni,1);
C   = sparse( [(1:n)'; (1:n)'], ni, ones(n,1)*[1 -1], n, G.cells.num);
%grad = -C;
div  = C';

Flux = div*state.flux(iface) - source;
pv   = poreVolume(G, rock);

state.s(:,1) = state.s(:,1) - (dT./pv).*Flux;
state.s(:,1)       = correct_saturations(state.s(:,1), opt.satwarn);
state.s(:,2) = 1-state.s(:,1);

end
%-------------------------------------------
function s = correct_saturations(s, satwarn)
   % Display error if s > 1+satwarn
   % or s < 0 - satwarn
   i = find(s(:,1) > 1 + satwarn);
   if ~isempty(i),
      disp('Saturation exceeds 1 in cells:')
      fprintf('%5d %g\n', [i, s(i,1)] .');
      error('simpleTransport: Error larger than satwarn')
   end

   i = find(s(:,1) < -satwarn);
   if ~isempty(i),
      disp('Saturation less than 0 in cells:')
      fprintf('%5d %g\n', [i, s(i,1)] .');
      error('simpleTransport: Error larger than satwarn')
   end
   % Correct numerical errors
   s(s(:,1) > 1, 1) = 1 - eps;
   s(s(:,1) < 0, 1) = 0;
end



