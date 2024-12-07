%% Nonlinear Glucose-Insulin Simulation with State Feedback and Luenberger Observer
clear; clc; close all;

%% Parameters
p1 = 0.03;   % 1/min
p2 = 0.02;   % 1/min
p3 = 0.01;   % 1/min
n  = 0.1;    % 1/min
Gb = 100;    % mg/dL
Ib = 10;     % mU/L

% Basal insulin infusion
u_basal = n * Ib;

%% Linearized Matrices (from previous derivation)
A = [ -p1    -Gb     0;
       0     -p2     p3;
       0      0      -n ];

B = [0;0;1];
C = [1 0 0];
E = [1;0;0];

%% Design Gains (Example)
% Choose controller poles (just an example)
controller_poles = [-0.05; -0.06; -0.08];
K = place(A,B,controller_poles);

% Choose observer poles (faster than controller)
observer_poles = [-0.1; -0.12; -0.15];
L = place(A',C',observer_poles)';

%% Simulation Setup
tfinal = 1440;  % 24 hours in minutes
tspan = [0 tfinal];

% Initial conditions
% Actual states: start at equilibrium (G=Gb, X=0, I=Ib)
G0 = Gb; X0 = 0; I0 = Ib;
x0 = [G0; X0; I0];

% Observer initial conditions: start with no knowledge of states
% The observer states represent deviations: [G-Gb; X; I-Ib]
xhat0 = [0;0;0];

% Augmented state vector: z = [G; X; I; xhat1; xhat2; xhat3]
z0 = [x0; xhat0];

options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[time, Z] = ode45(@(t,z) nonlinear_closed_loop_ode(t,z,p1,p2,p3,n,Gb,Ib,u_basal,A,B,C,K,L,@D_meal), tspan, z0, options);

% Extract results
G_sim = Z(:,1);
X_sim = Z(:,2);
I_sim = Z(:,3);

xhat_sim = Z(:,4:6);  % [xhat1; xhat2; xhat3]
xhat1 = xhat_sim(:,1); % estimated G-Gb
xhat2 = xhat_sim(:,2); % estimated X
xhat3 = xhat_sim(:,3); % estimated I-Ib

% Compute control input
u_pid = - (K * xhat_sim')'; % Δu = -K xhat
u = u_basal + u_pid;

%% Plot Results
figure;
subplot(3,1,1)
plot(time,G_sim,'LineWidth',2);
xlabel('Time (min)'); ylabel('Glucose (mg/dL)');
title('Nonlinear System Response: Glucose');

subplot(3,1,2)
plot(time,X_sim,'LineWidth',2);
xlabel('Time (min)'); ylabel('X (1/min)');
title('Insulin Action');

subplot(3,1,3)
plot(time,I_sim,'LineWidth',2);
xlabel('Time (min)'); ylabel('Insulin (mU/L)');
title('Insulin Concentration');

figure;
subplot(3,1,1)
plot(time,G_sim - Gb,'b','LineWidth',2); hold on;
plot(time,xhat1,'r--','LineWidth',2);
legend('Actual (G-Gb)','Estimated'); ylabel('G-Gb');
title('State and Observer Estimates');

subplot(3,1,2)
plot(time,X_sim,'b','LineWidth',2); hold on;
plot(time,xhat2,'r--','LineWidth',2);
legend('Actual X','Estimated X'); ylabel('X');

subplot(3,1,3)
plot(time,I_sim - Ib,'b','LineWidth',2); hold on;
plot(time,xhat3,'r--','LineWidth',2);
legend('Actual I-Ib','Estimated I-Ib'); ylabel('I - Ib');
xlabel('Time (min)');

figure;
plot(time,u,'k','LineWidth',2);
xlabel('Time (min)'); ylabel('u(t) (mU/min)');
title('Control Input (Insulin Infusion Rate)');


%% Nested Functions
function dz = nonlinear_closed_loop_ode(t,z,p1,p2,p3,n,Gb,Ib,u_basal,A,B,C,K,L,D_func)
    % z = [G; X; I; xhat1; xhat2; xhat3]

    G = z(1);
    X = z(2);
    I = z(3);
    xhat = z(4:6);

    % Disturbance
    D_val = D_func(t);

    % Output
    y = G; % measured output is glucose directly

    % Compute control input from observer
    u_pid = -K*xhat;
    u = u_basal + u_pid;

    % Nonlinear plant dynamics
    dGdt = -p1*(G - Gb) - X*G + D_val;
    dXdt = -p2*X + p3*(I - Ib);
    dIdt = -n*I + u;

    dx = [dGdt; dXdt; dIdt];

    % Observer dynamics
    % Observer is linear and uses deviations x = [G-Gb; X; I-Ib]
    % Let's define:
    x_dev = [G - Gb; X; I - Ib];
    y_dev = y - Gb; % y = G, so y_dev = G-Gb

    dxhat = A*xhat + B*(u - u_basal) + L*(y_dev - (C*xhat));

    dz = [dx; dxhat];
end

function D_val = D_meal(t)
    % Meal disturbance every 4 hours (240 min), lasting 30 min
    period = 240;   
    meal_duration = 30;
    D_amplitude = 1;
    time_in_period = mod(t, period);
    if time_in_period < meal_duration
       D_val = D_amplitude * sin(pi * time_in_period / meal_duration)^2;
    else
        D_val = 0;
    end
end
