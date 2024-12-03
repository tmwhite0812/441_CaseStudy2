% Define system parameters
p1 = 0.03; p2 = 0.02; p3 = 0.01; n = 0.1; % Example values
Gb = 100; Ib = 10; % Baseline glucose and insulin
D = 10; % Disturbance (e.g., meal glucose influx)

% State-space matrices (linearized system)
A = [-p1, -Gb, 0;
      0, -p2, p3;
      0,  0, -n];
B = [0; 0; 1];
E = [1; 0; 0];
C = [1, 0, 0];

% Design state feedback gain (K) using pole placement
desired_poles_controller = [-0.5, -0.6, -0.7]; % Desired closed-loop poles
K = place(A, B, desired_poles_controller);

% Design observer gain (L) using pole placement
desired_poles_observer = [-1.5, -1.6, -1.7]; % Faster poles for observer
L = place(A', C', desired_poles_observer)';

% Combined system matrices
A_cl = [A - B*K, -B*K;
        zeros(size(A)), A - L*C];
B_cl = [E; zeros(size(E))];
C_cl = [C, zeros(size(C))];
D_cl = 0;

% Simulate the system
sys_cl = ss(A_cl, B_cl, C_cl, D_cl);

% Time vector
t = 0:0.1:100;

% Disturbance (step input)
D_input = D * ones(size(t));

% Initial conditions for [x; e]
x0 = [120; 0; 0]; % Initial glucose, insulin action, plasma insulin
e0 = [0; 0; 0];   % Initial observer error
initial_state = [x0; e0];

% Simulate response
[y, t, x] = lsim(sys_cl, D_input, t, initial_state);

% Plot results
figure;
subplot(2, 1, 1);
plot(t, y, 'LineWidth', 2);
title('Glucose Concentration (G)');
xlabel('Time (min)');
ylabel('G (mg/dL)');
grid on;

subplot(2, 1, 2);
plot(t, x(:, 1), 'LineWidth', 2);
title('States (x)');
xlabel('Time (min)');
ylabel('State Variables');
legend('G', 'X', 'I');
grid on;

%%
clear; clc;

% Define system parameters
p1 = 0.03; p2 = 0.02; p3 = 0.01; n = 0.1; % Example values
Gb = 100; Ib = 10; % Baseline glucose and insulin
D = 10; % Disturbance (e.g., meal glucose influx)

% Equilibrium points (symbolically derived earlier)
u_eq = 0; % Example equilibrium input
I_eq = u_eq / n;
X_eq = (p3 / p2) * (I_eq - Ib);
G_eq = (D + p1 * Gb) / (p1 + X_eq);

% Linearized matrices (evaluate Jacobian at equilibrium)
A = [-p1 - X_eq, -G_eq, 0;
     0, -p2, p3;
     0, 0, -n];
B = [0; 0; 1];
E = [1; 0; 0];
C = [1, 0, 0];

% Design state feedback gain (K) using pole placement
desired_poles_controller = [-0.5, -0.6, -0.7]; % Desired closed-loop poles
K = place(A, B, desired_poles_controller);

% Design observer gain (L) using pole placement
desired_poles_observer = [-1.5, -1.6, -1.7]; % Faster poles for observer
L = place(A', C', desired_poles_observer)';

% Combined system matrices
A_cl = [A - B*K, -B*K;
        zeros(size(A)), A - L*C];
B_cl = [E; zeros(size(E))];
C_cl = [C, zeros(size(C))];
D_cl = 0;

% Simulate the system
sys_cl = ss(A_cl, B_cl, C_cl, D_cl);

% Time vector
t = 0:0.1:100;

% Disturbance (step input)
D_input = D * ones(size(t));

% Initial conditions for [x; e]
x0 = [120; 0; 0]; % Initial glucose, insulin action, plasma insulin
e0 = [0; 0; 0];   % Initial observer error
initial_state = [x0; e0];

% Simulate response
[y, t, x] = lsim(sys_cl, D_input, t, initial_state);

% Plot results
figure;
subplot(2, 1, 1);
plot(t, y, 'LineWidth', 2);
title('Glucose Concentration (G)');
xlabel('Time (min)');
ylabel('G (mg/dL)');
grid on;

subplot(2, 1, 2);
plot(t, x(:, 1), 'LineWidth', 2);
title('States (x)');
xlabel('Time (min)');
ylabel('State Variables');
legend('G', 'X', 'I');
grid on;

