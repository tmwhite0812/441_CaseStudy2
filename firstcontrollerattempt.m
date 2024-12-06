%% PI Control with State Observer Implementation
function main()
    % Parameters
    A = [-0.03, -100, 0;
          0, -0.02, 0.01;
          0, 0, -0.1]; % Example plant dynamics
    B = [0; 0; 1];     % Input matrix
    C = [1, 0, 0];     % Output matrix

    % Check controllability and observability
    if rank(ctrb(A, B)) ~= size(A, 1)
        error('System is not controllable');
    end
    if rank(obsv(A, C)) ~= size(A, 1)
        error('System is not observable');
    end

    % Desired closed-loop poles (for PI and observer)
    desired_PI_poles = [-5, -6, -7];       % Controller poles
    desired_observer_poles = [-8, -9, -10]; % Observer poles

    % Design state feedback gains (Kp)
    Kp = place(A, B, desired_PI_poles);

    % Add integral action (augment the system)
    Ai = [A, zeros(size(A, 1), 1);
         -C, 0];
    Bi = [B; 0];
    desired_augmented_poles = [-2, -3, -4, -5]; % PI control poles
    Ki = place(Ai, Bi, desired_augmented_poles);
    Kp_aug = Ki(1:end-1); % Extract state feedback gains
    Ki = Ki(end);         % Integral gain

    % Design observer gains (L)
    L = place(A', C', desired_observer_poles)';
    
    % Simulation parameters
    tspan = [0, 10];       % Simulation time
    x0 = [0.1; 0.1; 0.1];  % Initial states
    m0 = 0;                % Initial integral state
    xhat0 = [0; 0; 0];     % Initial state estimates
    r = 1;                 % Reference signal (step input)

    % Simulate the system
    [t, states] = ode45(@(t, x) system_dynamics(t, x, A, B, C, Kp_aug, Ki, L, r), tspan, [x0; m0; xhat0]);

    % Extract results
    x = states(:, 1:3);         % True states
    m = states(:, 4);           % Integral state
    xhat = states(:, 5:7);      % Estimated states
    y = C * x';                 % Output
    u = -Kp_aug * xhat' - Ki * m'; % Control input

    % Plot results
    figure;
    subplot(3, 1, 1);
    plot(t, y, 'b', 'LineWidth', 1.5);
    hold on;
    plot(t, r * ones(size(t)), 'r--', 'LineWidth', 1.5);
    title('Output Tracking');
    xlabel('Time (s)');
    ylabel('y(t)');
    legend('Output', 'Reference');

    subplot(3, 1, 2);
    plot(t, u, 'g', 'LineWidth', 1.5);
    title('Control Input');
    xlabel('Time (s)');
    ylabel('u(t)');

    subplot(3, 1, 3);
    plot(t, xhat(:, 1), 'r--', 'LineWidth', 1.5);
    hold on;
    plot(t, x(:, 1), 'b', 'LineWidth', 1.5);
    title('State Estimation');
    xlabel('Time (s)');
    ylabel('State x_1');
    legend('Estimated', 'True');
end

function dx = system_dynamics(t, x, A, B, C, Kp, Ki, L, r)
    % Extract states
    x_true = x(1:3);       % True states
    m = x(4);              % Integral state
    xhat = x(5:7);         % Estimated states

    % Compute control input
    u = -Kp * xhat - Ki * m;

    % True system dynamics
    dx_true = A * x_true + B * u;

    % Integral action dynamics
    dm = r - C * x_true;

    % Observer dynamics
    dxhat = A * xhat + B * u + L * (C * x_true - C * xhat);

    % Combine dynamics
    dx = [dx_true; dm; dxhat];
end
