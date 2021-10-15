function [grampc,Tsim,grampc_sdata] = initData()
% This function initializes a grampc struct in MATLAB and sets parameters 
% and options. In case of three output arguments the struct grampc_sdata for
% the use in Simulink is created as well. Define all options and parameters
% for the use of GRAMPC in MATLAB here.
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

%% Parameter definition

hor = 20;	       % MPC prediction and control horizon %%11 is almost half an orbit
Nsim = 4.0;        % number of simulated orbits
deltaT = 240.0;    % MPC sample time

% Initial values and setpoints of the states
user.param.x0    = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, ...  % MPC
                    1.0, 1.0, 1.0, 0.0, 0.0, 0.0];     % Lyapunov

user.param.xdes  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  ...  % MPC
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0];      % Lyapunov

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0, 0.0, 0.0];
user.param.udes  = [0.0, 0.0, 0.0];
user.param.umax  = [6.15384615384615e-06, 6.15384615384615e-06, 6.15384615384615e-06];
user.param.umin  = [-6.15384615384615e-06, -6.15384615384615e-06, -6.15384615384615e-06];

user.param.p0    = [0.0, 0.0, 0.0];

% Time variables
user.param.dt    = deltaT;          % Sampling time dt --> 240 sec due to slow actuators 
user.param.t0    = 0.0;             % time at the current sampling step
user.param.Thor  = deltaT * hor;    % Prediction horizon

%% Option definition

% Scaling values for the states and control inputs
user.opt.ScaleProblem = 'on';    
  
user.opt.xScale  = [1.00, 1.00, 1.00, 0.001, 0.001, 0.001, ...    % MPC
                    1.00, 1.00, 1.00, 0.001, 0.001, 0.001];       % Lyapunov  
user.opt.uScale  = [0.0001, 0.0001, 0.0001];   

% Solver options
user.opt.Nhor        = 100;              % Number of steps for the system integration (necessary for 'heun' and 'euler' integrator) (default = 30)
user.opt.MaxGradIter = 30;               % Maximum number of gradient iterations (default = 2)
user.opt.Integrator = 'ruku45';          % alternatives: 'heun', 'euler' for fixed Nhor or 'ruku45', 'rodas' for variable Nhor (default = 'heun')
user.opt.MaxMultIter = 5;                % default: 1

% Constraint options
user.opt.InequalityConstraints = 'on';
user.opt.ConstraintsHandling = 'auglag';
user.opt.ConstraintsAbsTol = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00];     % all zeros = hard constraints
user.opt.MultiplierMax =10;

%% User parameter definition 

Q_x     = 10*[1.0, 1.0, 1.0];    % weights on x, y, z
Q_xdot	= 0*[1.0, 1.0, 1.0] ;    % weights on xdot, ydot, zdot
         
R = [1.0, 1.0, 1.0];     % weights on u_x, u_y, u_z

stateConstraints = [-2000, 2000, -2000, 2000, -2000, 2000, ...      % min-value and max-value of x, y, z
                    -2, 2, -2, 2, -2, 2];                           % min-value and max-value of xdot, ydot, zdot

userparam = [Q_x, Q_xdot, R, stateConstraints];
         
%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin (optional)
%grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
a = 6378 + 500;                                      % semi-major axis a = earth's radius + orbit hight
TOrbit = sqrt((4*pi^2*a^3)/(5.9722e24*6.6743e-20));  % orbital period in seconds
if nargout>1
    Tsim = Nsim*TOrbit;  
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
