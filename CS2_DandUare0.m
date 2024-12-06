%% ESE 441 Case Study 2
% Keeler Tardiff and Tyler White
%% The autonomous (u(t) and D(t) = 0) system.

% Define symbolic variables
syms G(t) X(t) I(t) p1 p2 p3 n Gb Ib 
syms G_e X_e I_e % Equilibrium variables

% Define the nonlinear equations
eq1 = diff(G(t), t) == -p1 * (G(t) - Gb) - X(t) * G(t);
eq2 = diff(X(t), t) == -p2 * X(t) + p3 * (I(t) - Ib);
eq3 = diff(I(t), t) == -n * I(t);

% Substitute equilibrium conditions into the equations
eq1_eq = subs(eq1, [diff(G(t), t), G(t), X(t), I(t)], [0, G_e, X_e, I_e]);
eq2_eq = subs(eq2, [diff(X(t), t), G(t), X(t), I(t)], [0, G_e, X_e, I_e]);
eq3_eq = subs(eq3, [diff(I(t), t), G(t), X(t), I(t)], [0, G_e, X_e, I_e]);

% Solve for equilibrium points
[sol_G_e, sol_X_e, sol_I_e] = solve([eq1_eq, eq2_eq, eq3_eq], [G_e, X_e, I_e]);

% Display equilibrium points in readable format
disp('Equilibrium Points:');
disp('G_e =');
pretty(sol_G_e);
disp('X_e =');
pretty(sol_X_e);
disp('I_e =');
pretty(sol_I_e);

% Define state variables and inputs
syms G X I u

% Rewrite the system in state-space form
f1 = -p1 * (G - Gb) - X * G; % Glucose dynamics
f2 = -p2 * X + p3 * (I - Ib); % Insulin action dynamics
f3 = -n * I; % Insulin dynamics

% Define the state vector
states = [G; X; I];

% Compute the Jacobian matrix (A matrix)
A = jacobian([f1; f2; f3], states);

% Substitute equilibrium points into A
A_linearized = subs(A, [G, X, I], [sol_G_e, sol_X_e, sol_I_e]);

% Display the linearized state matrix in readable format
disp('Linearized A Matrix:');
pretty(A_linearized);

% Numerical parameter values (given in the problem)
p1_val = 0.03; % min^-1
p2_val = 0.02; % min^-1
p3_val = 0.01; % min^-1
n_val = 0.1; % min^-1
Gb_val = 100; % mg/dL
Ib_val = 10; % mU/L

% Substitute the numerical values into the linearized A matrix
A_numeric = double(subs(A_linearized, ...
    [p1, p2, p3, n, Gb, Ib], ...
    [p1_val, p2_val, p3_val, n_val, Gb_val, Ib_val]));

% Display the numeric A matrix
disp('Numeric Linearized A Matrix:');
disp(A_numeric);

% Compute the eigenvalues of the numeric A matrix
eigenvalues = eig(A_numeric);

% Display the eigenvalues
disp('Eigenvalues of the Numeric Linearized A Matrix:');
disp(eigenvalues);

%% Simulate the linearized autonomous system
% Define the time span and initial conditions
tspan = [0 100]; % Simulation time from 0 to 100 minutes
x0_linear = [1; 0; 0]; % Small perturbation in glucose, no initial perturbation in X or I

% Linearized system dynamics
linearized_system = @(t, x) A_numeric * x;

% Solve the linearized system using ode15s
[t_linear, x_linear] = ode15s(linearized_system, tspan, x0_linear);

% Add baseline to linearized system results
G_linear = x_linear(:, 1) + Gb_val; % Add baseline glucose
X_linear = x_linear(:, 2); % Insulin action deviation
I_linear = x_linear(:, 3) + Ib_val; % Add baseline insulin

%% Simulate the nonlinear autonomous system
% Initial conditions for the nonlinear system
x0_nonlinear = [Gb_val + 1; 0; Ib_val]; % Slight perturbation in glucose, baseline for others

% Nonlinear system dynamics
nonlinear_system = @(t, x) [
    -p1_val * (x(1) - Gb_val) - x(2) * x(1); % Glucose dynamics
    -p2_val * x(2) + p3_val * (x(3) - Ib_val); % Insulin action dynamics
    -n_val * x(3)                            % Plasma insulin dynamics
];

% Solve the nonlinear system using ode15s
[t_nonlinear, x_nonlinear] = ode15s(nonlinear_system, tspan, x0_nonlinear);

% Extract nonlinear system results
G_nonlinear = x_nonlinear(:, 1); % Glucose
X_nonlinear = x_nonlinear(:, 2); % Insulin action
I_nonlinear = x_nonlinear(:, 3); % Plasma insulin

%% Combined Plots
figure;

% Glucose
subplot(3, 1, 1);
plot(t_linear, G_linear, 'b--', 'LineWidth', 1.5); hold on;
plot(t_nonlinear, G_nonlinear, 'r', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Glucose (G)');
legend('Linearized', 'Nonlinear');
title('Glucose Dynamics');
grid on;

% Insulin Action
subplot(3, 1, 2);
plot(t_linear, X_linear, 'b--', 'LineWidth', 1.5); hold on;
plot(t_nonlinear, X_nonlinear, 'r', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Insulin Action (X)');
legend('Linearized', 'Nonlinear');
title('Insulin Action Dynamics');
grid on;

% Plasma Insulin
subplot(3, 1, 3);
plot(t_linear, I_linear, 'b--', 'LineWidth', 1.5); hold on;
plot(t_nonlinear, I_nonlinear, 'r', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Plasma Insulin (I)');
legend('Linearized', 'Nonlinear');
title('Plasma Insulin Dynamics');
grid on;

%% Separate Plots for Individual Systems (Optional)
figure;

% Nonlinear system
plot(t_nonlinear, G_nonlinear, 'r', 'LineWidth', 1.5); hold on;
plot(t_nonlinear, X_nonlinear, 'g', 'LineWidth', 1.5);
plot(t_nonlinear, I_nonlinear, 'b', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Nonlinear States');
legend('G (Glucose)', 'X (Insulin Action)', 'I (Plasma Insulin)');
title('Nonlinear Autonomous System');
grid on;

figure;

% Linearized system
plot(t_linear, G_linear, 'r--', 'LineWidth', 1.5); hold on;
plot(t_linear, X_linear, 'g--', 'LineWidth', 1.5);
plot(t_linear, I_linear, 'b--', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Linearized States');
legend('G (Glucose)', 'X (Insulin Action)', 'I (Plasma Insulin)');
title('Linearized Autonomous System');
grid on;
