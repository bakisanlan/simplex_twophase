% -----------------------------------------------------------------------
% This scripts is created as example of using two phase simplex method.
% -----------------------------------------------------------------------

% Constructing of problem constraints and objective function.
A = [1 -13 0; 3 0 8; 0 -1 5 ];  % A is rearrenged by [gt;lt]
b = [2 ; 5 ; 1.5];
f = [-1 1 0];


[x_opt, z_opt] = simplex_ineq(f,A,b,3,0);

disp(['Optimal x values are: ', num2str(x_opt)]);
disp(['Optimal z value is: ', num2str(z_opt)]);




