clc; clear; close all;

% Symbol declarations
syms ia ib il ua ub uc p1 p2 p3 us n L R C real
syms duc dil s  % duc = d(uc)/dt, dil = d(il)/dt

% Equation system
eq1 = ia == -n*ib;            
eq2 = ub == n*ua;             
eq3 = ua == p3 - p2;         
eq4 = ub == uc - p2;          
eq5 = (p1-us)/(2*R) + (p1-(uc+L*dil))/R + (p1-p3)/R == 0;  
eq6 = ((uc+L*dil)-p1)/R + il == 0;                        
eq7 = (p3-p1)/R + ia == 0;                                
eq8 = p2/R - ia - ib == 0;                                 
eq9 = uc/(2*R) - il - C*duc == 0;

% Solve for all variables including derivatives
sol = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9], ...
            [ia, ib, p1, p2, p3, ua, ub, duc, dil]);

% Extract state derivatives
duc_dt = sol.duc;
dil_dt = sol.dil;

% State-space representation:
% States: x1 = uc, x2 = il
% Input: u = us
% Output: y = il

% The state derivatives are already in terms of uc, il, and us
% We need to express them in standard state-space form:
% dx/dt = A*x + B*u

% Collect terms for state matrix A
A11 = diff(duc_dt, uc);
A12 = diff(duc_dt, il);
A21 = diff(dil_dt, uc);
A22 = diff(dil_dt, il);

A = [A11, A12;
     A21, A22];

% Collect terms for input matrix B
B1 = diff(duc_dt, us);
B2 = diff(dil_dt, us);
B = [B1; B2];

% Output matrix (we're measuring il)
C = [0 1];
D = 0;

% Simplify matrices
A = simplify(A);
B = simplify(B);

% Display results
disp('State derivatives:');
disp('duc/dt = ');
pretty(duc_dt);
disp('dil/dt = ');
pretty(dil_dt);

disp('State matrix A:');
pretty(A);

disp('Input matrix B:');
pretty(B);

disp('Output matrix C:');
disp(C);

disp('Feedthrough matrix D:');
disp(D);

% Calculate transfer function il/us
[numerator, denominator] = numden(C*((s*eye(2) - A)\B));
transfer_function = numerator/denominator;
disp('Transfer function H(s) = il/us:');
pretty(simplify(transfer_function));

% Optional: Substitute numerical values for a specific case
% n_val = 1; R_val = 1000; L_val = 1e-3; C_val = 1e-6;
% A_num = double(subs(A, [n, R, L, C], [n_val, R_val, L_val, C_val]));
% B_num = double(subs(B, [n, R, L, C], [n_val, R_val, L_val, C_val]));
% sys = ss(A_num, B_num, C, D);
% tf(sys)
