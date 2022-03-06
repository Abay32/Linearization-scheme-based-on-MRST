close all
clear all
clc
nx = 50;
[ny,nz] = deal(1,1);

G = cartGrid([nx,ny], [4,1]);
G = computeGeometry(G);

%% Set rock and fluid data
% rock.perm = repmat(3*milli*darcy, [G.cells.num, 1]);
% rock.poro = repmat(0.1, [G.cells.num, 1]);
per = [4.2025*ones(nx/2,1);0.5625*ones(nx/2,1)];
rock = struct('perm',    per, ...
              'poro',   1*ones(G.cells.num, 1))%100*milli*

x = linspace(0, 1, 500).';
y = linspace(1, 0, 500).';

pc_form = 'nonwetting';
cap_scale = 1;
sor = 0.0; swi = 0.0; krwout = 0.5; kroout = 0.8;
son = (1-x-sor)/(1-swi-sor);
swn = (x-swi)/(1-swi-sor);
pc  =  x.^-2;
n = 2;
m = 3\2;
pcl =  (x.^-(1/m) - 1).^(1-m);
pcr =  (x.^-(1/m) - 1).^(1-m);
% n = 2;
% m = 1-1/n;
relperm = [x.^2,  y.^2];
relperm = [sqrt(x).*(1-(1-x.^(1/m)).^m).^2,  sqrt(1-x).*(1-(x).^(1/m)).^(2*m)];
[kr, pcl]= tabulatedSatFunc([x,relperm, pcl.*cap_scale]);
[kr, pcr]= tabulatedSatFunc([x,relperm, pcr.*cap_scale]);
%, 
props   = constantProperties([1 ,  1], ...
                           [1080, 600]);
fluid   = struct('properties', props                  , ...
                  'saturation'  , @(x, varargin)    x.s  , ...
                  'relperm'   , kr  );

fluid.pcl = @(x) sqrt(x.^-1 - 1).*cap_scale*barsa;
fluid.pcr = @(x) sqrt(x.^-1 - 1).*cap_scale*barsa;

fluid_c = struct('properties', props                     , ...
                  'saturation', @(x, varargin) x.s     , ...
                  'relperm'   , kr                     , ...                  
                  'pcl'        , @(x, varargin) pcl(x.s), ...
                  'pcr'        , @(x, varargin) pcr(x.s));           
%% Compute transmissibility and initialize reservoir state
% hT = simpleComputeTrans(G, rock);
%[mu,rho] = fluid.properties();
ijk      = gridLogicalIndices(G);
upInd    = (ijk{1} > nx/2);
state    = initState(G, [], 0, [1,0]);
state.s(upInd,1) = 0;state.s(upInd,2)=1;
state.s(~upInd,:)

%
rSol_pc  = initState(G, [], 50*barsa, [0.2,0.8]);
%
%% Impose boundary conditions
% Our flow solvers automatically assume no-flow conditions on all outer
% (and inner) boundaries; other type of boundary conditions need to be
% specified explicitly.
[bc,src] = deal([]);
% c   = (47:1:53).';
% src = addSource([], c, 0.01*ones(size(c))*kilogram ./ day(), 'sat', [0.0,1]);
% display(src);


pv = poreVolume(G, rock);

injRate = -sum(pv)/(500*day);
% bc = pside(bc, G, 'xmin', 40*barsa , 'sat', [1, 0]);
% bc = pside(bc, G, 'xmax', 5*barsa,   'sat', [0.7, 0.3]);
bc = fluxside(bc, G, 'LEFT', -injRate, 'sat',[1,0]);
bc = pside(bc, G, 'RIGHT', 0*barsa, 'sat', [0 1]);


gravity off
verbose = false;
tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, ...
                                                fluid, 'bc', bc,'src',src, ...
                                                'verbose', verbose);
T      = 2;
dTplot = 80*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(G,rock);
dT     = T/4000;
t      = 0;
plotNo = 1;
h1 = 'No pc - '; H2 = 'Linear pc -';
e = []; 
p_org = []; 
p_pc  = [];
figure;
%
while t < T    
    state  = simpleIncompTPFAModefied(state, G, rock, fluid,  dT, 'bc', bc, 'src', src);
    %rSol_pc = simpleIncompTPFAModefied(rSol_pc, G, rock, fluid_pc,  'bc', bc, 'src', src);
    %
    %state  = simpletransport(G, state, rock, dT,'src',src,'bc',bc);    
    state   = tsolve(state,dT,fluid);    
    t      = t + dT;
end
%
xvals = linspace(0,100, G.cells.num);
figure(1)
 
plot(state.s(:,1),'--','linewidth',2)
legend('CO2 saturation after 250 days')
xlabel('Distance [m]')
ylabel('Saturtaion Distribution')


figure(2)
if isfield(fluid,'pc')
    po = state.pressure + fluid.pc(state);
else
    po = state.pressure;
end
%
plot([state.pressure(:,1),po]/1e+6,'--','linewidth',2)
legend('CO2 Pressure after 250 days' )
xlabel('Distance [m]')
ylabel('Pressure   [MPa]')