%% Glucose-Insulin Dynamics Simulation (Nonlinear and Linearized)
clear; clc; close all;

%% Parameters (Nominal)
p1 = 0.03;   % 1/min
p2 = 0.02;   % 1/min
p3 = 0.01;   % 1/min
n  = 0.1;    % 1/min
Gb = 100;    % mg/dL
Ib = 10;     % mU/L

% Simulation time (24 hours = 1440 minutes)
tfinal = 1440;   % minutes
tspan = [0 tfinal];

% Basal insulin delivery u(t). For simplicity, let's assume constant basal insulin input:
u_basal = n*Ib;  % This keeps insulin at Ib in steady state if no other disturbances.

% Define the nonlinear model as a function handle for ODE45
% State vector: x = [G; X; I]
% D(t): meal disturbance function defined below
% u(t) = u_basal (constant)
ode_nonlinear = @(t,x) nonlinear_ode(t,x,p1,p2,p3,n,Gb,Ib,u_basal);

% Solve the nonlinear system
x0 = [Gb; 0; Ib];  % initial conditions at equilibrium
[tnl, xnl] = ode45(ode_nonlinear, tspan, x0);

% Extract nonlinear results
G_nl = xnl(:,1);
X_nl = xnl(:,2);
I_nl = xnl(:,3);

%% Linearized System
% Linearization around (G=Gb, X=0, I=Ib)
A_lin = [ -p1,     -Gb,     0;
           0,       -p2,    p3;
           0,        0,     -n ];

E_d = [1; 0; 0]; 

x0_lin = [0;0;0];  % deviation from equilibrium
ode_linear = @(t,x) A_lin*x + E_d*D_meal(t);

[tlin, xlin] = ode45(ode_linear, tspan, x0_lin);

% Reconstruct actual states for linear model:
G_lin = xlin(:,1) + Gb;
X_lin = xlin(:,2) + 0;
I_lin = xlin(:,3) + Ib;

%% Plot Results
figure;

% Left column: Nonlinear
subplot(3,2,1)
plot(tnl, G_nl, 'b','LineWidth',2);
xlabel('Time (min)'); ylabel('Glucose (mg/dL)');
title('Nonlinear Glucose Response');

subplot(3,2,3)
plot(tnl, X_nl, 'b','LineWidth',2);
xlabel('Time (min)'); ylabel('X (1/min)');
title('Nonlinear Insulin Action');

subplot(3,2,5)
plot(tnl, I_nl, 'b','LineWidth',2);
xlabel('Time (min)'); ylabel('Insulin (mU/L)');
title('Nonlinear Insulin Concentration');

% Right column: Linearized
subplot(3,2,2)
plot(tlin, G_lin, 'r--','LineWidth',2);
xlabel('Time (min)'); ylabel('Glucose (mg/dL)');
title('Linearized Glucose Response');

subplot(3,2,4)
plot(tlin, X_lin, 'r--','LineWidth',2);
xlabel('Time (min)'); ylabel('X (1/min)');
title('Linearized Insulin Action');

subplot(3,2,6)
plot(tlin, I_lin, 'r--','LineWidth',2);
xlabel('Time (min)'); ylabel('Insulin (mU/L)');
title('Linearized Insulin Concentration');

%% Nested Functions

function dx = nonlinear_ode(t,x,p1,p2,p3,n,Gb,Ib,u_basal)
    G = x(1); X = x(2); I = x(3);
    
    D_val = D_meal(t);
    u = u_basal; 
    
    dGdt = -p1*(G - Gb) - X*G + D_val;
    dXdt = -p2*X + p3*(I - Ib);
    dIdt = -n*I + u;
    
    dx = [dGdt; dXdt; dIdt];
end

function D_val = D_meal(t)
    % Adjust the number and spacing of meals:
    % We want 4 meals over 24 hours (1440 min), so a meal every 360 minutes.
    period = 360;       % interval between meals
    meal_duration = 30; % duration of each meal spike
    D_amplitude = 1;    % increase amplitude to emphasize spikes
    
    time_in_period = mod(t, period);
    if time_in_period < meal_duration
        D_val = D_amplitude;
    else
        D_val = 0;
    end
end
