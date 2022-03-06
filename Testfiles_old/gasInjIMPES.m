close all
clear all
clc
nx = 300;
[ny,nz] = deal(1,1);

G = cartGrid([nx,ny], [100,1]*meter);
G = computeGeometry(G);

%% Set rock and fluid data
% rock.perm = repmat(3*milli*darcy, [G.cells.num, 1]);
% rock.poro = repmat(0.1, [G.cells.num, 1]);
rock = struct('perm',   milli*darcy*ones(G.cells.num, 1), ...
              'poro',   0.1*ones(G.cells.num, 1));%100*milli*

x = linspace(0, 1, 500).';
y = linspace(1, 0, 500).';

pc_form = 'nonwetting';
cap_scale = 20;
sor = 0.0; swi = 0.0; krwout = 0.5; kroout = 0.8;
son = (1-x-sor)/(1-swi-sor);
swn = (x-swi)/(1-swi-sor);
pc  =  x.^-2;
n = 2;
m = 1-1/n;
Pc =  (x.^-(1/m) - 1).^(1/n);
  
relperm = [sqrt(x).*(1-(1-x.^(1/m)).^m).^2,  sqrt(1-x).*(1-(x).^(1/m)).^(2*m)];
[kr, pc]= tabulatedSatFunc([x, relperm, Pc.*cap_scale*barsa]);
%, 
props   = constantProperties([.86 ,  9e-3]*centi*poise, ...
                             [1080, 600].*kilogram/meter^3);
% fluid   = struct('properties', props                  , ...
%                   'saturation'  , @(x, varargin)    x.s  , ...
%                   'relperm'   , kr  );

fluid.pc = @(x) sqrt(x.^-1 - 1).*cap_scale*barsa;

fluid = struct('properties', props                  , ...
                'saturation', @(x, varargin)    x.s  , ...
                'relperm'   , kr                     , ...
                'pc'        , @(x, varargin) pc(x.s));           
%% Compute transmissibility and initialize reservoir state
% hT = simpleComputeTrans(G, rock);
%[mu,rho] = fluid.properties();
state    = initState(G, [], 5*barsa, [0.7, 0.3]);

%% Impose boundary conditions
% Our flow solvers automatically assume no-flow conditions on all outer
% (and inner) boundaries; other type of boundary conditions need to be
% specified explicitly.
[bc,src] = deal([]);
% c   = (47:1:53).';
% src = addSource([], c, 0.01*ones(size(c))*kilogram ./ day(), 'sat', [0.0,1]);
% display(src);
if isfield(fluid,'pc')
    Pc =  @(x)(x.^-(1/m) - 1).^(1/n);
    PR = 5 + Pc(0.7);
else
    PR = 0;
end
 
bc = pside(bc, G, 'xmin', 40*barsa , 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', PR*barsa,   'sat', [0.7, .3]);

gravity off
verbose = false;
tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, ...
                                                fluid, 'bc', bc,'src',src, ...
                                                'verbose', verbose);
T      = 45*day();
dTplot = 80*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(G,rock);
dT     = T/8000;
t      = 0;
plotNo = 1;
h1 = 'No pc - '; H2 = 'Linear pc -';
e = []; 
p_org = []; 
p_pc  = [];
figure;
%
Th = computeTrans(G,rock);
theta = 0;
while t < T    
 %   state = incompTPFAphaseFormulation(state, G, rock, Th, fluid, dT, theta, 'bc', bc, 'src', src)   ;
    state  = simpleIncompTPFAModefiedold(state, G, rock, fluid, 'bc', bc, 'src', src);
    state  = simpletransportold(G, state, rock, dT,'src',src,'bc',bc);   
    t      = t + dT;
end
%
xvals = linspace(0,100, G.cells.num);
% plot(xvals,pw_pc,'--','linewidth',2)
% legend('CO2 pressure[mPa]')
% 
% 
figure(1) 
plot(xvals,state.s(:,1),'--','linewidth',2)
legend('H_2O saturation after 45 days')
xlabel('Distance [m]')
ylabel('Saturtaion Distribution')
axis([0 100 0.65 1.05])

figure(2)
if isfield(fluid,'pc')
    po = state.pressure + fluid.pc(state);
else
    po = state.pressure;
end
%
plot(xvals, [state.pressure(:,1),po]/1e+6,'--','linewidth',2)
legend('H_2O Pressure after 45 days' ,'CO_2 Pressure after 45 days' )
xlabel('Distance [m]')
ylabel('Pressure   [MPa]')
axis([0 100 0 4])

figure(3)
plot(state.flux,'--','linewidth',2);