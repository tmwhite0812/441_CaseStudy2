function main()
    % Parameters
    p1 = 0.03; % min^-1
    p2 = 0.02; % min^-1
    p3 = 0.01; % min^-1
    n = 0.1;   % min^-1
    Gb = 100;  % mg/dL
    Ib = 10;   % mU/L

    % Simplified Linearized System Matrix A
    A = [
        -p1, -Gb, 0;
        0, -p2, p3;
        0, 0, -n
    ];

    % Input and Output Matrices
    B = [0; 0; 1]; % Example input matrix (adjust as needed)
    C = [1, 0, 0]; % Output matrix

    % Display the Simplified Linearized System Matrix
    disp("Simplified Linearized System Matrix A:");
    disp(A);

    % Check Controllability
    disp("The Rank of Controllability Matrix:");
    controllability(A, B);

    % Check Observability
    disp("The Rank of Observability Matrix:");
    observability(A, C);
end

function controllability(A, B)
    % Compute controllability matrix
    ctrl_mtrx = [B, A * B, A^2 * B];
    rank_controllability = rank(ctrl_mtrx);
    disp(rank_controllability);
    if rank_controllability == size(A, 1)
        disp("The controllability matrix is full rank.");
        disp(ctrl_mtrx);
    else
        disp("The matrix is not full rank.");
    end
end

function observability(A, C)
    % Compute observability matrix
    obs_mtrx = [C; C * A; C * A^2];
    rank_observability = rank(obs_mtrx);
    disp(rank_observability);
    if rank_observability == size(A, 1)
        disp("The observability matrix is full rank.");
        disp(obs_mtrx);
    else
        disp("The matrix is not full rank.");
    end
end
