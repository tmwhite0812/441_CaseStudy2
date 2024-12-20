%% Tyler and Keeler Case Study 2
%% Script to Check Minimality of Transfer Function
clear; clc; close all;

%all base values
p1 = 0.03; 
p2 = 0.02; 
p3 = 0.01; 
n = 0.1;  
Gb = 100;  
Ib = 10;   

% Linearized Matrices
A = [-p1, -Gb, 0;
      0, -p2, p3;
      0,   0, -n];
B = [0; 0; 1];
E = [1; 0; 0];
C = [1, 0, 0];

% Controller and Observer Poles
desired_poles_controller = [-0.05; -0.06; -0.08];
desired_poles_observer = [-0.1; -0.12; -0.15];

% pole placement
K = place(A, B, desired_poles_controller);
L = place(A', C', desired_poles_observer)';

% Augmented Matrices
A_aug = [A - B*K, zeros(3, 3);
         L*C, A - B*K - L*C];
E_aug = [E; zeros(3, 1)];
C_aug = [C, zeros(1, 3)];

% Convert to State-Space and Transfer Function
sys_aug = ss(A_aug, E_aug, C_aug, 0);
G_tf = tf(sys_aug);

disp('Transfer Function (Numeric):');
disp(G_tf);

% Check Minimality
[z, p, k] = zpkdata(G_tf, 'v'); % Extract zeros, poles, and gain
disp('Poles of the Transfer Function:');
disp(p);
disp('Zeros of the Transfer Function:');
disp(z);

% Check for Pole-Zero Cancellations
tolerance = 1e-6; % tolerance for decimal points on p/z
cancellations = [];
for i = 1:length(p)
    for j = 1:length(z)
        if abs(p(i) - z(j)) < tolerance
            cancellations = [cancellations; p(i)]; 
        end
    end
end

if isempty(cancellations)
    disp('The transfer function is minimal (no pole-zero cancellations).');
else
    disp('The transfer function is NOT minimal (pole-zero cancellations found):');
    disp(cancellations);
end

