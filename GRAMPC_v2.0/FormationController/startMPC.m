function [vec,grampc,figNr] = startMPC(figNr,compile,varargin)
% This function runs the BallOnPlate example. It compiles the c-Code,
% initializes the grampc struct, runs the simulation and plots the results.
%
% Input arguments are:
% 1) figNr - number of the first plot
% 2) compile - flag whether to compile the whole toolbox or/and the problem function: 
%              1: only the problem function is compiled
%              2: the whole toolbox and the problem function is compiled
%              else or empty input: nothing is compiled
% 3 - end) - flags for the compilation (e.g. 'debug' or 'verbose') see 
%            make.m in the matlab folder for more details
%
%
% This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
%
% GRAMPC -- A software framework for embedded nonlinear model predictive
% control using a gradient-based augmented Lagrangian approach
%
% Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
% Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0). 
% Developed at the Institute of Measurement, Control, and Microtechnology,
% Ulm University. All rights reserved.
%
% GRAMPC is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as 
% published by the Free Software Foundation, either version 3 of 
% the License, or (at your option) any later version.
%
% GRAMPC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public 
% License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
%

%% Check input arguments
if nargin < 1 || isempty(figNr)
    figNr = 1;
end
if nargin < 2 || isempty(compile)
    compile = 0;
end

%% Parameters
% path to grampc root
grampc_root_path = '../';
addpath([grampc_root_path 'matlab/mfiles']);
% name of problem function
probfct = 'probfct_formation-control.c';

% plot predicted trajectories
PLOT_PRED = 1;
% plot solution trajectories
PLOT_TRAJ = 1;
% plot optimization statistics
PLOT_STAT = 1;
% update plots after N steps
PLOT_STEPS = 100;
% pause after each plot
PLOT_PAUSE = 0;

% Options for the reference simulation
odeopt =  [];


%% Compilation
% compile toolbox
if compile > 1 || ~exist([grampc_root_path 'matlab/bin'], 'dir')
    grampc_make_toolbox(grampc_root_path, varargin{:});
end
% compile problem
if compile > 0 || ~exist('+CmexFiles', 'dir')
    grampc_make_probfct(grampc_root_path, probfct, varargin{:});
end


%% Initialization
% init GRAMPC and print options and parameters
[grampc,Tsim] = initData;
CmexFiles.grampc_printopt_Cmex(grampc);
CmexFiles.grampc_printparam_Cmex(grampc);

% init solution structure
vec = grampc_init_struct_sol(grampc, Tsim);

% init plots and store figure handles
if PLOT_PRED
    phpP = grampc_init_plot_pred(grampc,figNr);
    figNr = figNr+1;
end
if PLOT_TRAJ
    phpT = grampc_init_plot_sim(vec,figNr);
    figNr = figNr+1;
end
if PLOT_STAT
    phpS = grampc_init_plot_stat(vec,grampc,figNr);
    figNr = figNr+1;
end

% initalization values for the Lyapunov controller (only necessary for the plots)
xdes_lyap = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
lyapCosts(1:length(vec.t)) = 0;
mpcCosts(1:length(vec.t)) = 0;
Vdot_lyap(1:length(vec.t)) = 0;
Vdot_mpc(1:length(vec.t)) = 0;
P = 0.001*[1077.86996565576		 -2.18720481141238e-14	464.598059590277		262442.240474333		-2.64696126889320e-12	387485.312310805		;
          -2.18720481141238e-14	 795.061209829976		-7.18623514600261e-14	-1.31435150012107e-11	292134.658374445		-2.54898786454275e-11	;
          464.598059590277		 -7.18623514600261e-14	1569.39820749626		-40035.2497803359		-1.28583463802585e-10	689473.441363119		;
          262442.240474333		 -1.31435150012107e-11	-40035.2497803359		116295882.947240		1.56598682560876e-08	37586413.9960404		;
          -2.64696126889320e-12	 292134.658374445		-1.28583463802585e-10	1.56598682560876e-08	296169690.670132		-5.03244672061450e-08	;
          387485.312310805		 -2.54898786454275e-11	689473.441363119		37586413.9960404		-5.03244672061450e-08	397874483.646451	   ];
Pdot = -0.001*eye(6);

%% MPC loop
i = 1;
while 1
    % set current time and current state
    grampc = CmexFiles.grampc_setparam_Cmex(grampc,'t0',vec.t(i));
    grampc = CmexFiles.grampc_setparam_Cmex(grampc,'x0',vec.x(:,i));
    
    % run MPC and save results
    [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex(grampc);
    vec = grampc_update_struct_sol(grampc, vec, i);
    
    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status,'Error');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
    
    % check for end of simulation
    if i+1 > length(vec.t)  || vec.t(i+1) >Tsim
        break;
    end
    
    % reset Lyapunov controller states, if prediction horizon is reached
    if (mod(vec.t(i),vec.T(i)) == 0)
        vec.x(7:12,i) = vec.x(1:6,i);
    end    
    
    % simulate system
    [~,xtemp] = ode45(@CmexFiles.grampc_ffct_Cmex,vec.t(i)+[0 double(grampc.param.dt)],vec.x(:,i),odeopt,vec.u(:,i),vec.p(:,i),grampc.userparam);
    vec.x(:,i+1) = xtemp(end,:);
        
    % for u_dot constraints, p is set to the last value of u:
    if i > 1
        vec.p(:,i) = vec.u(:,i-1);
    end
    
    % evaluate time-dependent constraints to obtain h(x,u,p) instead of max(0,h(x,u,p))
    vec.constr(:,i) = CmexFiles.grampc_ghfct_Cmex(vec.t(i), vec.x(:,i), vec.u(:,i), vec.p(:,i), grampc.userparam);
    
    % calculate Lyapunov costs V(x) and V_dot(x)
    deltaX = vec.x(7:12, i) - xdes_lyap';
    lyapCosts(i) = (deltaX' * P * deltaX);
    Vdot_lyap(i) = (deltaX' * Pdot * deltaX);
    deltaX2 = vec.x(1:6, i) - xdes_lyap';
    mpcCosts(i) = (deltaX2' * P * deltaX2);
    Vdot_mpc(i) = (deltaX2' * Pdot * deltaX2);
    
    % update iteration counter
    i = i + 1;
    
    % plot data
    if mod(i,PLOT_STEPS) == 0 || i == length(vec.t)
        if PLOT_PRED
            grampc_update_plot_pred(grampc,phpP);
        end
        if PLOT_TRAJ
            grampc_update_plot_sim(vec,phpT);
        end
        if PLOT_STAT
            grampc_update_plot_stat(vec,grampc,phpS);
        end
        drawnow
        if PLOT_PAUSE
            pause;
        end
    end
end


%% Plots for the Lyapunov controller

figure(4)
%subplot(2,2,[1,3]);
subplot(2,2,1);
plot(mpcCosts)
hold on
plot(lyapCosts)
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
title('Comparison of Lyapunov function V(x)')

subplot(2,2,3);
plot(Vdot_mpc)
hold on
plot(Vdot_lyap)
hold off
legend({'$\bf \dot{V}(x_{MPC})$', '$\bf \dot{V}(x_{Lyap})$'}, 'Interpreter','latex')
xlabel('Time in [dt]')
title('Comparison of derivative of V(x)')


subplot(2,2,2);
plot(vec.x(7,:))
hold on
plot(vec.x(8,:))
hold on
plot(vec.x(9,:))
hold off
legend({'x', 'y', 'z'})
xlabel('Time in [dt]')
ylabel('relative position in [m]')
title('Lyapunov controller: Positions')

subplot(2,2,4);
plot(vec.x(10,:))
hold on
plot(vec.x(11,:))
hold on
plot(vec.x(12,:))
hold off
%legend({'xdot (m/s)', 'ydot (m/s)', 'zdot (m/s)'})
xlabel('Time in [dt]')
ylabel('relative velocity in [m/s]')
legend({'$\bf \dot{x}$', '$\bf \dot{y}$', '$\bf \dot{z}$'}, 'Interpreter','latex')
title('Lyapunov controller: Velocities')

%% Plot comparison of the states (x_MPC vs. x_Lyap)
figure(5)
subplot(3,2,1)
plot(vec.x(1,:))
hold on
plot(vec.x(7,:))
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
ylabel('x-position in [m]')
title('Comparison of relative position in x-direction')

subplot(3,2,2)
plot(vec.x(4,:))
hold on
plot(vec.x(10,:))
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
ylabel('x-velocity in [m/s]')
title('Comparison of relative velocity in x-direction')

subplot(3,2,3)
plot(vec.x(2,:))
hold on
plot(vec.x(8,:))
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
ylabel('y-position in [m]')
title('Comparison of relative position in y-direction')

subplot(3,2,4)
plot(vec.x(5,:))
hold on
plot(vec.x(11,:))
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
ylabel('y-velocity in [m/s]')
title('Comparison of relative velocity in y-direction')

subplot(3,2,5)
plot(vec.x(3,:))
hold on
plot(vec.x(9,:))
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
ylabel('z-position in [m]')
title('Comparison of relative position in z-direction')

subplot(3,2,6)
plot(vec.x(6,:))
hold on
plot(vec.x(12,:))
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
ylabel('z-velocity in [m/s]')
title('Comparison of relative velocity in z-direction')

figure(6)
subplot(1,2,1);
plot(mpcCosts)
hold on
plot(lyapCosts)
hold off
legend({'MPC', 'Lyapunov'})
xlabel('Time in [dt]')
title('Comparison of Lyapunov function V(x)')

subplot(1,2,2);
plot(Vdot_mpc)
hold on
plot(Vdot_lyap)
xlabel('Time in [dt]')
hold off
legend({'$\bf \dot{V}(x_{MPC})$', '$\bf \dot{V}(x_{Lyap})$'}, 'Interpreter','latex')
title('Comparison of derivative of V(x)')
