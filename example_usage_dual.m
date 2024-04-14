% Example use of solve_dual_lp function
A = [1 -13.5 0; 3 0 8; 0 -1 5 ; -1.5 36 -7.5];  % A is rearrenged by [gt;lt]
b = [2 ; 5 ; 1.5 ; 4.5];
f = [1 -3 0];
n_gt = 3;                                 %number of >= constraints
n_lt = 1; 

[y_opt, z_dual_opt] = simplex_dual_ineq(f, A, b, n_gt, n_lt);

disp('Optimal Dual Solution (y):');
disp(y_opt);
disp('Optimal Dual Objective Value (z_dual):');
disp(z_dual_opt);
