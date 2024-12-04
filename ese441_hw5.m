function main()
    gamma = 2;
    A1 = [0 1 0; 
        10*(1 + gamma^2) -1 200; 
        0 -1 -4];

    A2 = [0 1 0; 
        10*(gamma^2-1) -1 200; 
        0 -1 -4];

    B = [0; 0; 20];
    C = [0 1 0];
    
    disp("The Rank of Controllability matrix for x_e1 is");
    controllability(A1, B);
    
    disp("The Rank of Observability matrix for x_e1 is");
    observability(A1, C);
    
    disp("Solutions for K gains for x_e1:");
    solveK(A1, gamma);
    
    disp("Solutions for L gains for x_e1:");
    solveL(A1, C, gamma);
    
    disp("The Rank of Controllability matrix for x_e2 is");
    controllability(A2, B);
    
    disp("The Rank of Observability matrix for x_e2 is");
    observability(A2, C);
    
    disp("Solutions for K gains for x_e2:");
    solveK(A2, gamma);
    
    disp("Solutions for L gains for x_e2:");
    solveL(A2, C, gamma);
end

function controllability(A, B)
    cnt_mtrx = [B, A*B, A^2*B];
    rank_controllability = rank(cnt_mtrx);
    disp(rank_controllability);
    if rank_controllability == size(A, 1)
        disp("The controllability matrix is full rank.");
        disp(cnt_mtrx);
    else
        disp("The matrix is not full rank.");
    end
end

function observability(A, C)
    obs_mtrx = [C; C*A; C*A^2];
    rank_observability = rank(obs_mtrx);
    disp(rank_observability);
    if rank_observability == size(A, 1)
        disp("The observability matrix is full rank.");
        disp(obs_mtrx);
    else
        disp("The matrix is not full rank.");
    end
end

function solveK(A, gamma)
    syms k1 k2 k3 s;
    
    A_sym = sym(A);

    G = A_sym;
    G(3, 1) = -20*k1;
    G(3, 2) = -1 - 20*k2;
    G(3, 3) = -4 - 20*k3;
    
    % characteristic equation is s^3 + 9s^2 + 27s + 27
    poly_fb = charpoly(G, s);
    coeff_fb = coeffs(poly_fb, s, 'All');
    a2_fb = coeff_fb(end-2);
    a1_fb = coeff_fb(end-1);
    a0_fb = coeff_fb(end);

    eqn1 = a2_fb == 9;
    eqn2 = a1_fb == 27;
    eqn3 = a0_fb == 27;
    
    solution_fb = solve([eqn1, eqn2, eqn3], [k1, k2, k3]);
    disp(solution_fb);
    
    G_verified = subs(G, [k1, k2, k3], [solution_fb.k1, solution_fb.k2, solution_fb.k3]);
    disp(G_verified);
    eigenvalues = eig(G_verified);
    disp('Eigenvalues of the verified closed-loop matrix:');
    disp(eigenvalues);
    
    if all(abs(double(eigenvalues + 3)) < 1e-6)
        disp('The eigenvalues match the desired value of -3.');
    else
        disp('The eigenvalues do NOT match the desired value of -3.');
    end
end


function solveL(A, C, gamma)
    syms L1 L2 L3 s;
    
    A_sym = sym(A);
    
    % 
    M = A_sym;
    M(1, 2) = 1 - L1;
    M(2, 2) = -1 - L2;
    M(3, 2) = -1 - L3;
    
    % characteristic equation is s^3 + 30s + 300s +1000
    poly_obs = charpoly(M, s);
    coeff_obs = coeffs(poly_obs, s, 'All');
    a2_obs = coeff_obs(end-2);
    a1_obs = coeff_obs(end-1);
    a0_obs = coeff_obs(end);
    
    eqn1_obs = a2_obs == 30;
    eqn2_obs = a1_obs == 300;
    eqn3_obs = a0_obs == 1000;
    
    solution_obs = solve([eqn1_obs, eqn2_obs, eqn3_obs], [L1, L2, L3]);
    disp(solution_obs);
    
    % verification of polynomial 
    M_obs_verified = subs(M, [L1, L2, L3], [solution_obs.L1, solution_obs.L2, solution_obs.L3]);
    disp(M_obs_verified);
    poly_obs_verified = charpoly(M_obs_verified, s);
    disp('Characteristic polynomial of the verified observer matrix:');
    disp(poly_obs_verified);
    
    desired_poly_obs = s^3 + 30*s^2 + 300*s + 1000;
    if simplify(poly_obs_verified - desired_poly_obs) == 0
        disp('The polynomial matches the desired polynomial s^3 + 30s^2 + 300s + 1000.');
    else
        disp('The polynomial does NOT match the desired polynomial.');
    end
end

