%% Glucose-Insulin Dynamics: State Feedback with Luenberger Observer
clear; clc; close all;

%% Parameters (Nominal)
p1 = 0.03;   % 1/min
p2 = 0.02;   % 1/min
p3 = 0.01;   % 1/min
n  = 0.1;    % 1/min
Gb = 100;    % mg/dL (baseline glucose)
Ib = 10;     % mU/L (baseline insulin)

% Basal insulin infusion to keep I=Ib at steady state
u_basal = n*Ib;  

%% Equilibrium Analysis:
% At equilibrium: G=Gb, X=0, I=Ib, D=0, and u = u_basal.
% Define state vector deviations: x = [G - Gb; X; I - Ib]

%% Linearized State-Space Model:
% From derivations:
% x1 = G - Gb
% x2 = X
% x3 = I - Ib
% Inputs: Δu = u - u_basal
% Disturbance: D(t)

A = [ -p1,        -Gb,      0;
       0,          -p2,     p3;
       0,           0,      -n ];

B = [0;0;1];      % Input affects I (and thus x3)
E = [1;0;0];      % Disturbance input affects G-equation
C = [1 0 0];      % Output y = G = Gb + x1

% The output equation: y = Gb + x1, but for state-feedback and observer
% design, we use y - Gb = x1. Since adding Gb just shifts the measured output,
% it doesn't affect observer design.

%% Design the State-Feedback Controller K
% We want to place the closed-loop poles for A - B*K somewhere in the left half-plane.
% Choose desired closed-loop poles (faster than natural dynamics):
desired_poles_controller = [-0.05; -0.06; -0.08]; % adjust as needed
K = place(A,B,desired_poles_controller);

%% Design the Luenberger Observer L
% We know y = Cx + Gb. For the observer, we only measure y (or equivalently x1).
% Choose observer poles:
desired_poles_observer = [-0.1; -0.12; -0.15]; % faster than controller poles
L = place(A',C',desired_poles_observer)';  % L is transpose of solution from place(A',C',...)

% Now we have K for feedback and L for observation.

%% Closed-Loop with Observer Setup
% The combined observer-controller system:
% Controller: Δu = -K * x_hat
% Observer: \dot{x_hat} = A x_hat + B Δu + L(y - C x_hat)
%
% Disturbance D(t): Suppose periodic meals every 4 hours (240 min),
% producing a disturbance in G. Let's implement a simple square pulse.

% Simulation parameters
tfinal = 1440;  % 24 hours in minutes
tspan = [0 tfinal];

% Initial states:
% Actual system initial condition near equilibrium:
x0 = [0; 0; 0];   % start at equilibrium (x1=0, x2=0, x3=0)
% Observer initial condition:
xhat0 = [0;0;0];  % observer starts with no initial state knowledge error
% Combine actual and observer states into one augmented system:
% Augmented state: z = [x; x_hat], dimension 6

z0 = [x0; xhat0];

%% Simulate the Closed-Loop System with Observer
% We'll define a function that implements the closed-loop dynamics:
% Actual system: dx/dt = A x + B(-K x_hat) + E D(t)
% Observer: dx_hat/dt = A x_hat + B(-K x_hat) + L(y - C x_hat) with y = C x + Gb
%
% Remember u = u_basal + Δu, but Δu = -K x_hat. The baseline u_basal doesn't affect states directly in the linearized form except shifting equilibrium.

options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[time, Z] = ode45(@(t,z) closed_loop_ode(t,z,A,B,C,E,K,L,Gb,u_basal,@D_meal), tspan, z0, options);

% Extract results
X = Z(:,1:3);       % Actual states [x1; x2; x3]
X_hat = Z(:,4:6);   % Estimated states
G_state = X(:,1) + Gb;   % Actual glucose
I_state = X(:,3) + Ib;   % Actual insulin
U = u_basal - (K*X_hat')'; % Control input

%% Plot Results
figure;
subplot(3,1,1);
plot(time, G_state,'LineWidth',2); hold on;
xlabel('Time (min)'); ylabel('Glucose (mg/dL)');
title('Glucose Response with State Feedback + Observer');
ylim([95 125]); 
subplot(3,1,2);
plot(time, X(:,2),'LineWidth',2); 
xlabel('Time (min)'); ylabel('X (1/min)');
title('Insulin Action');
subplot(3,1,3);
plot(time, I_state,'LineWidth',2);
xlabel('Time (min)'); ylabel('Insulin (mU/L)');
title('Insulin Concentration');


figure;
subplot(2,1,1)
plot(time,X(:,1),'b','LineWidth',2); hold on;
plot(time,X_hat(:,1),'r--','LineWidth',2);
legend('x_1 (actual)','x_1 hat (estimated)');
ylabel('x_1');
title('State and Observer Estimates');


subplot(2,1,2)
plot(time,U,'k','LineWidth',2);
xlabel('Time (min)'); ylabel('u(t) (mU/min)');
title('Control Input (Insulin Infusion Rate)');



%% Nested Functions
function dz = closed_loop_ode(t,z,A,B,C,E,K,L,Gb,u_basal,D_func)
    % z = [x; x_hat]
    x = z(1:3);
    x_hat = z(4:6);

    % Disturbance
    D_val = D_func(t);

    % Output
    y = C*x + Gb;

    % Control input
    u_pid = -K*x_hat;
    u = u_basal + u_pid;  % total insulin

    % Actual system dynamics
    dx = A*x + B*(u - u_basal) + E*D_val;

    % Observer dynamics
    % y - C x_hat = measurement residual
    dx_hat = A*x_hat + B*(u - u_basal) + L*(y - (C*x_hat + Gb));

    dz = [dx; dx_hat];
end

function D_val = D_meal(t)
    % A simple meal disturbance every 240 min lasting 30 min
    period = 360;
    meal_duration = 30;
    D_amplitude = 2;  % Smaller amplitude for demonstration
    time_in_period = mod(t, period);
    if time_in_period < meal_duration
       D_val = D_amplitude * sin(pi * time_in_period / meal_duration)^2;
    else
        D_val = 0;
    end
end
