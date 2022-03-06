mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up model
% We start by generating a model object that describes the reservoir. To
% construct this model, we need three more fundamental structures: 'G'
% represents the grid with reservoir geometry, 'rock' holds the
% petrophysical properties, and 'fluid' gives the fluid properties. In
% addition to the model object, we need a 'state' struct that defines the
% initial state in the reservoir (pressure and fluid saturations and
% compositions).

% Construct 3D grid with 50 cells in the x-direction
G = cartGrid([200, 1], [100, 1]*meter);
G = computeGeometry(G);
gravity off
% Homogenous rock properties
rock = struct('perm', 1*milli*darcy*ones(G.cells.num, 1), ...
              'poro', 0.1*ones(G.cells.num, 1));%100*milli*

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'WO', ... 
                           'mu',  [7.86,9e-3]*centi*poise,...
                           'rho', [1080,400]*kilogram/meter^3,...
                           'n',   [2.0 2.0]);%, ...
                           %'c',   [0,0]/barsa

n = 2.0;
m = 1-1/n;
PcO = @(s) 3e6*((s).^-(1/m) -1).^(1/n); 
fluid.pcOW = PcO;

% Set up model and initial state.
model          = TwoPhaseOilWaterModel(G, rock, fluid);
%

state0         = initResSol(G, 43*barsa, [0.7, 0.3]);%2*barsa
state0.wellSol = initWellSolAD([], model, state0);

% x = 0:0.01:1;
% y = 1:-0.01:0;
% [krW, krO] = model.evaluateRelPerm({x, y})
% plot(x,krO)
% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
pv      = poreVolume(G, rock);
injRate = sum(pv)/(500*day);

%bc = fluxside([], G, 'xmin', injRate, 'sat', [1, 0]);
a = 75*barsa + fluid.pcOW(0.7)/4 ;
bc = pside([], G, 'xmin', 82.652*barsa ,   'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 43*barsa, 'sat', [0.0, 1]);

%% Simulate 1 PVI using a manual loop
% There are several ways to run a simulation. A simple approach is to use a
% manual loop, where you explicitly call a nonlinear solver to solve each
% implicit time step
solver = NonLinearSolver();

dT        = 1*day;
n         =46;
states    = cell(n+1, 1);
states{1} = state0;
for i = 1:n
    state       = solver.solveTimestep(states{i}, dT, model, 'bc', bc);
    states{i+1} = state;
end

%% Plot the result
% We set up a plot using plotToolbar from the mrst-gui module. Since the
% problem is essentially one dimensional, we plot the water saturation
% (first column of the "s" field in state) as a 1D plot.

% % close all
% % plotToolbar(G, states, 'field', 's:2', 'plot1d', true, ...
% %                        'lockCaxis', true, 'startplayback', true);

%% Repeat simulation with general solver
% To use the general solver, we first need to set up a schedule that
% describes the time steps and the drive mechanisms (wells, boundary
% conditions, and source terms) that are active in each time step. In
% addition, one can specify various forms of time-step control. Here,
% however, we simply rely on the default setup
schedule    = simpleSchedule(repmat(dT,1,n), 'bc', bc);
[~,sstates] = simulateScheduleAD(state0, model, schedule);

pressureG   = sstates{n}.pressure - fluid.pcOW(sstates{n}.s(:,1));

close all
% figure(1)
% plotToolbar(G, sstates,'field', 's:1','lockCaxis',true), 
% caxis([0 1])%, view(10,10)
% colorbar

%plot(G,sstates{25}.pressure)
%% Repeat simulation with visualization
% The general solver has a hook, that enables you to visualize the progress
% of the simulation (and stop it and continue running it in 'debug' mode).
%close all    
% % fn = getPlotAfterStep(state0, model, schedule, ...
% %     'plotWell', false, 'plotReservoir', true, 'field', 's:1', ...
% %     'lockCaxis',true, 'plot1d', true);

%[~,sstates,report] = ...
%   simulateScheduleAD(state0, model, schedule);%,'afterStepFn', fn

pressures  = [sstates{n}.pressure,pressureG]/2e6;

figure(1)
xvals = linspace(0,100,model.G.cells.num);
plot(xvals,pressures,'--','linewidth',2)
legend('CO2 pressure[mPa]','Water pressure[mPa]')
xlabel('Distance [m]')
ylabel('Pressure   [MPa]')
axis([0 100 0 4])
% % 
% % n = [10,15,20,25];figure(1)
% % for i = 1:length(n)    
figure(2)
plot(xvals,sstates{n}.s(:,1),'--','linewidth',2)
% %     hold on
legend('CO2 saturation at 45 days')
xlabel('Distance [m]')
ylabel('saturation')
axis([0 100 0.65 1.05])
% % end
%    figure(3)
%pcOW = fluid.pcOW(sstates{n}.s(:,1));
%plot(xvals,pcOW,'linewidth',2)