% -----------------------------------------------------------------------
% This scripts is created as example of using two phase simplex method.
% -----------------------------------------------------------------------

%%
% A = [1 -13 0; 3 0 8; 0 -1 5 ; -6 9 -7];  % A is rearrenged by [gt;lt]
% b = [2 ; 5 ; 1.5 ; 4.5];
% f = [-1 1 0];
% n_gt = 3;                                 %number of >= constraints
% n_lt = 1;  
%%
% [x_opt, z_opt] = simplex_ineq(f,A,b,n_gt,n_lt);
% A = [1 -13.5 0; 3 0 8; 0 -1 5 ; -1.5 36 -7.5];  % A is rearrenged by [gt;lt]
% b = [2 ; 5 ; 1.5 ; 4.5];
% f = -[1 -3 0];
% n_gt = 3;                                 %number of >= constraints
% n_lt = 1; 

%%
A = [1 0 ;1 1 ; 4 10 ];
b = [3 ; 7 ; 40];
f = [30 100];
n_gt = 1;
n_lt = 2;
[x_opt, z_opt] = simplex_ineq(f,A,b,n_gt,n_lt);


disp(['Optimal x values are: ', num2str(x_opt)]);
disp(['Optimal z value is: ', num2str(z_opt)]);




