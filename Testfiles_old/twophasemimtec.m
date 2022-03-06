%nx = 10; ny = 10; nz  = 1;

G = cartGrid([300, 1], [1000, 1]*meter);
G = computeGeometry(G);

rock.perm = repmat(milli.*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3        , [G.cells.num, 1]);
hT = computeTrans(G, rock);

x = linspace(0, 1, 500).';
y = linspace(1, 0, 500).';
pc_form = 'wetting';
cap_scale = 20;
pc = sqrt(x.^-2 - 1);
%pc  =  x.^-2;
n = 2;
m = 1-1/n;
relperm = [x.^2,  y.^2];
relperm = [sqrt(x).*(1-(1-x.^(1/m)).^m).^2,  sqrt(y).*(1-(1-y).^(1/m)).^(2*m)];
% krw = @(x) sqrt(x).*(1-(1-x.^(1/m)).^m).^2
% krn = @(x) sqrt(1-x).*(1-(x).^(1/m)).^(2*m)
%  
% plot(x,relperm,'linewidth',2)
[kr, pc]  = tabulatedSatFunc([x,relperm, pc.*cap_scale*barsa]);
props     = constantProperties([1,3].*centi*poise, ...
                           [1080, 600].*kilogram/meter^3);
fluid = struct('properties', props                  , ...
               'saturation', @(x, varargin)    x.s  , ...
               'relperm'   , kr);
           
fluid = struct('properties', props                  , ...
                  'saturation', @(x, varargin)    x.s  , ...
                  'relperm'   , kr                     , ...
                  'pc'        , @(x, varargin) pc(x.s));

% rate = 0.5*meter^3/day;
% bhp  = 1*barsa;

% W = verticalWell([], G, rock, 1, 1, 1:nz,          ...
%                  'Type', 'rate', 'Val', rate, ...
%                  'Radius', .1, 'Name', 'I', 'Comp_i', [1 0]);
% W = verticalWell(W, G, rock, nx, ny, 1:nz,     ...
%                  'Type','bhp', 'Val', bhp, ...
%                  'Radius', .1, 'Dir', 'x', 'Name', 'P', 'Comp_i', [0 1]);
% rSol    = initState(G, W, 0, [0.2, 0.8]);
% rSol_pc = initState(G, W, 0, [0.2, 0.8]);
rSol         = initState(G,[], 5*barsa, [1,0]);
rSol_pc      = initState(G,[], 50*barsa,[0,1]);


% bc = pside([], G, 'xmin', 400*barsa , 'sat', [1, 0]);
% bc = pside(bc, G, 'xmax', 50*barsa,   'sat', [0, 0]);

pv      = poreVolume(G, rock);
injRate = -sum(pv)/(500*day);
bc      = fluxside([], G, 'xmin', -injRate, 'sat', [1,0]);
bc      = pside(bc, G, 'xmax', 0*barsa, 'sat', [0,1]);

% bc = pside([], G, 'xmin', 40*barsa ,   'sat', [1, 0]);
% bc = pside(bc, G, 'xmax', 5*barsa, 'sat', [0.7, .3]);


gravity off
verbose = false;

S  = computeMimeticIP(G, rock, 'Verbose', verbose,'InnerProduct','ip_tpf');
psolve  = @(state, fluid) incompMimetic(state, G, S, fluid, 'bc',bc);
tsolve  = @(state, dT, fluid) explicitTransport(state, G, dT, rock, ...
                                                fluid, 'bc', bc, ...
                                                'verbose', verbose);
rSol    = psolve(rSol, fluid);
%rSol_pc = psolve(rSol_pc, fluid_pc);
T       = 45*day();
dT      = T/4000;
dTplot  = 100*day();  % plot only every 100th day
N       = fix(T/dTplot);
pv      = poreVolume(G,rock);
t       = 0; 
plotNo  = 1;
h1      = 'No pc - '; 
H2      = 'Linear pc - ';
e       = []; 
p_org   = []; 
p_pc    = [];
figure;


while t < T,
   % TRANSPORT SOLVE
   rSol    = psolve(rSol,    fluid);
   rSol    = tsolve(rSol, dT, fluid);
   %rSol_pc = tsolve(rSol_pc, dT, fluid_pc);

   % Check for inconsistent saturations
   s = [rSol.s(:,1); rSol_pc.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   
   %rSol_pc = psolve(rSol_pc, fluid_pc);
   t = t + dT;
end

xvals = linspace(0,100, G.cells.num);
figure(1)
 
plot(xvals, rSol.s,'--','linewidth',2)
legend('Water saturation after 45 days','CO2 saturation after 45 days with Cap. Pressure')
xlabel('Distance [m]')
ylabel('Saturtaion Distributions')
grid on

figure(2)
po = rSol.pressure + fluid.pc(rSol);
plot(xvals,  [rSol.pressure(:,1),po]/1e+6, '--','linewidth',2)
legend('Water pressure after 45 days','CO2 pressure after 45 days with Cap. Pressure')
xlabel('Distance [m]')
ylabel('Pressure Distributions [MPa]')
grid on