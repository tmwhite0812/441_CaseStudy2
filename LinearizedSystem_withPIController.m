%% Glucose-Insulin Dynamics Simulation (Nonlinear and Linearized with PI Control)
clear; clc; close all;

%% Parameters (Nominal)
p1 = 0.03;   % 1/min
p2 = 0.02;   % 1/min
p3 = 0.01;   % 1/min
n  = 0.1;    % 1/min
Gb = 100;    % mg/dL
Ib = 10;     % mU/L

% Basal insulin
u_basal = n*Ib;

% Simulation time
tfinal = 1440;
tspan = (0:.1:tfinal);

% Gain ranges
Kp_vals = linspace(0.0005, 0.05, 5); % Test 5 values of Kp
Ki_vals = linspace(0.00001, 0.001, 5); % Test 5 values of Ki
Kd = 0; % Keep Kd = 0 for PI control

% Initialize storage for results
responses = struct();

%% Iterate through Kp and Ki
for i = 1:length(Kp_vals)
    for j = 1:length(Ki_vals)
        % Current gains
        Kp = Kp_vals(i);
        Ki = Ki_vals(j);

        % Linearized Model
        A_lin = [ -p1 - (p3 * (Ib * n)) / (n * p2), -n * p2 * (Gb * p1) / (p3 * Ib * n + n * p1 * p2),  0;
                   0,  -p2, p3;
                   0,   0,  -n];


        B_u = [0;0;1];  % Input affects I directly
        E_d = [1;0;0];

        % Augmented state: x_aug = [x_1; x_2; x_3; x_i]
        x0_lin = [0;0;0;0]; % Start with no deviation and no integral error

        % Solve the linear system with PI control
        [tlin, xlin_aug] = ode45(@(t,x_aug) linear_pi_ode(t,x_aug,A_lin,B_u,E_d,Kp,Ki,u_basal,Gb), tspan, x0_lin);

        % Extract glucose and insulin states
        x1_lin = xlin_aug(:,1); % G-Gb
        x2_lin = xlin_aug(:,2); % Insulin action
        x3_lin = xlin_aug(:,3); % I-Ib

        G_lin = x1_lin + Gb;    % Glucose
        X_lin = x2_lin;         % Insulin action
        I_lin = x3_lin + Ib;    % Insulin concentration
        
        % Store results
        responses(i, j).Kp = Kp;
        responses(i, j).Ki = Ki;
        responses(i, j).time = tlin;
        responses(i, j).glucose = G_lin;
        responses(i, j).insulin_action = X_lin;
        responses(i, j).insulin_concentration = I_lin;
    end
end


%% Plot Results
% Plot Glucose Dynamics
figure;
for i = 1:length(Kp_vals)
    for j = 1:length(Ki_vals)
        subplot(length(Kp_vals), length(Ki_vals), (i-1)*length(Ki_vals) + j);
        plot(responses(i, j).time, responses(i, j).glucose, 'r', 'LineWidth', 1.5);
        title(['Kp = ', num2str(responses(i, j).Kp), ', Ki = ', num2str(responses(i, j).Ki)]);
        xlabel('Time (min)');
        ylabel('Glucose (mg/dL)');
        ylim([80 140]);
    end
end
sgtitle('Glucose Dynamics');

% Plot Insulin Action Dynamics
figure;
for i = 1:length(Kp_vals)
    for j = 1:length(Ki_vals)
        subplot(length(Kp_vals), length(Ki_vals), (i-1)*length(Ki_vals) + j);
        plot(responses(i, j).time, responses(i, j).insulin_action, 'b', 'LineWidth', 1.5);
        title(['Kp = ', num2str(responses(i, j).Kp), ', Ki = ', num2str(responses(i, j).Ki)]);
        xlabel('Time (min)');
        ylabel('Insulin Action (X, 1/min)');
    end
end
sgtitle('Insulin Action Dynamics');

% Plot Insulin Plasma Concentration
figure;
for i = 1:length(Kp_vals)
    for j = 1:length(Ki_vals)
        subplot(length(Kp_vals), length(Ki_vals), (i-1)*length(Ki_vals) + j);
        plot(responses(i, j).time, responses(i, j).insulin_concentration, 'g', 'LineWidth', 1.5);
        title(['Kp = ', num2str(responses(i, j).Kp), ', Ki = ', num2str(responses(i, j).Ki)]);
        xlabel('Time (min)');
        ylabel('Insulin Conc. (I, mU/L)');
    end
end
sgtitle('Insulin Plasma Concentration Dynamics');


%% Function Definitions
function dx_aug = linear_pi_ode(t, x_aug, A_lin, B_u, E_d, Kp, Ki, u_basal, Gb)
    % x_aug = [x_1; x_2; x_3; x_i]
    x_1 = x_aug(1); % Glucose deviation
    x_2 = x_aug(2); % Insulin action
    x_3 = x_aug(3); % Insulin concentration deviation
    x_i = x_aug(4); % Integral of error

    % Error: e = Gb - G = -x_1 since G = Gb + x_1
    e = -x_1;

    % Anti-windup for integral term
    x_i_clamped = min(max(x_i, -20), 20); % Clamp integral term

    % PI control with clamped integral term
    u_pid = min(max(Kp * e + Ki * x_i_clamped, 0), 15); % Clamp between 0 and 15 mU/L% Clamp u_pid
    u = u_basal + u_pid; % Total insulin input

    % Disturbance
    D_val = D_meal(t);

    % Linearized dynamics
    x_lin = [x_1; x_2; x_3];
    dx_lin = A_lin * x_lin + E_d * D_val + B_u * (u - u_basal);

    % Integral state dynamics
    dx_i = e;

    dx_aug = [dx_lin; dx_i];
end

function D_val = D_meal(t)
    period = 360;             % Period of meals (minutes)
    meal_duration = 30;       % Approximate meal duration (minutes)
    D_amplitude = 100;        % Maximum amplitude of the disturbance
    time_in_period = mod(t, period); % Time within the current period

    % Use a smoothed function (e.g., sine wave) for gradual rise and fall
    if time_in_period < meal_duration
        % Gradual rise and fall during the meal
        D_val = D_amplitude * sin(pi * time_in_period / meal_duration)^2;
    else
        % Return to baseline outside of meal periods
        D_val = 0;
    end
end

