function [x_opt,z_opt] = simplex_ineq(f,A,b,n_gt,n_lt)

% This function takes maximizing objective function as f, both greater 
% than zero and less then zero inequality constraints matrix as A, 
% right hand side of constraints as b. And also number of >= inequalities
% as n_gt and number of <= inequalities should be given as input.
% Then gives output as [x_opt,z_opt]

% ------------------------------------------------------------------
% Example use of simplex_eq function  
% A = [1 -13 0; 3 0 8; 0 -1 5 ; -6 9 -7];  % A is rearrenged by [gt;lt]
% b = [2 ; 5 ; 1.5 ; 4.5];
% f = [-1 1 0];
%n_gt = 3;                                 %number of >= constraints
%n_lt = 1;                                 %number of <= constraints
% [x_opt,z_opt] = simplex_eq(f,A,b)
% ------------------------------------------------------------------

%% 
% A = [1 -13 0; 3 0 8; 0 -1 5 ; -6 9 -7];  % A is rearrenged by [gt;lt]
% b = [2 ; 5 ; 1.5 ; 4.5];
% f = [-1 1 0];

% A = [1 -13 0; 3 0 8; 0 -1 5 ; -6 9 -7; 1 0 0 ; 0 1 0 ; 0 0 1];  
% b = [2 ; 5 ; 1.5 ; 4.5 ; 3 ; 3 ; 3];
% f = -[1 -1 0];

% n_gt = 3; % number of greater than/equal zero cons.
% n_lt = 1; % number of leAss than/equal zero cons.
% n_gt = 3; % number of greater than/equal zero cons.
% n_lt = 4; % number of leAss than/equal zero cons.

[A_row, A_col] = size(A);

S_var_matrix = diag([-ones(1,n_gt), ones(1,n_lt)]); % S variable matrix

A_var_matrix = [eye(n_gt) ; zeros(n_lt,n_gt)];      % A variable matrix

w_matrix = [-sum(A(1:n_gt,:),1) , ones(1,n_gt) , zeros(1,n_lt) , zeros(1,n_gt)];

RHS = [b ; -sum(b(1:n_gt))];

tableau = [[A , S_var_matrix, A_var_matrix ; w_matrix] , RHS];  % tableau without basis and original cost obj. function

basis = [A_col+A_row+1: A_col+A_row+n_gt, (A_col+n_gt+1:A_col+n_gt+n_lt)]; % basis index of each const Equations

% initial x_opt z_opt
x_opt = zeros(1,A_col);
z_opt = 0;

%% Phase 1
while any(tableau(end,1:end-1) < 0)
    [tableau,basis] = pivotter(tableau,basis);
end

%checking are there any negative RHS column element
tableau = RHS_posneg(tableau);

%% % phase 2
% checking are there any basis var as artificial var
% Case 1 --> w==0 and there are no basis artificial vars in objective row
if ~any(ismember(basis,A_col+A_row+1:A_col+A_row+n_gt)) && tableau(end,end) < 10e-12 
    fprintf("Optimal Solution Found.\nw=0 and there are no basis artificial vars in objective row.\n")

    % replacing w with original f objective function
    % delete artificial var column
    tableau(:,A_col+A_row+1:A_col+A_row+n_gt) = [];
    tableau(end,:) = [-f, zeros(1,A_row+1)];

    % forcing pivot column rule to added f row
    for j=basis
        if tableau(end,j) ~= 0
            row_factor = tableau(end,j);
            tableau(end,:) = tableau(end,:) - row_factor * tableau(find(basis==j),:);
        end
    end

    % simplex method algorithm, finding most negative elements in f
    % objective row
    while any(tableau(end,1:end-1) < 0)
        [tableau,basis] = pivotter(tableau,basis);
    end

    % finding x_opt values from basis and tableau
    non_zero_opt_var = [basis', tableau(1:end-1,end)];
    non_zero_opt_x_var = non_zero_opt_var(ismember(non_zero_opt_var(:,1),1:A_col),:);
    x_opt(non_zero_opt_x_var(:,1)') =  non_zero_opt_x_var(:,2)';

    % finding z_opt
    z_opt = tableau(end,end);

% Case 2 --> w==0 and there is at least one basis artificial var in
% objective row
elseif tableau(end,end) < 10e-12
    fprintf('Optimal Solution Found.\nw=0 and there is basis artificial vars in objective row.\n')

    % storing basis artificial var.
    art_var = A_col+A_row+1:A_col+A_row+n_gt ;
    basis_art_var_idx = ismember(basis,art_var);
    basis_art_var = basis(basis_art_var_idx);

    % storing non basis artifical var.
    non_basis_art_var_idx = ~ismember(art_var,basis_art_var);
    non_basis_art_var = art_var(non_basis_art_var_idx);

    % dropping(replace with zeros) non basis artifical var columns
    if ~isempty(non_basis_art_var)
        tableau(:,non_basis_art_var) = zeros(A_row,1);
    end

    % dropping(replace with zeros) non basis original var with positive coef. in objective function 
    coef_of_org_var_in_obj = tableau(end,1:A_col);
    org_var_pos_coef_in_obj_logical= (coef_of_org_var_in_obj > 0);
    org_var_pos_coef_in_obj = find(org_var_pos_coef_in_obj_logical == 1);
    tableau(:,org_var_pos_coef_in_obj) = zeros(size(tableau(:,org_var_pos_coef_in_obj)));

    % replacing w with original f objective function
    tableau(end,:) = [-f, zeros(1,A_row+n_gt+1)];

    % forcing pivot column rule to added f row
    for j=basis
        if tableau(end,j) ~= 0
            row_factor = tableau(end,j);
            tableau(end,:) = tableau(end,:) - row_factor * tableau(find(basis==j),:);
        end
    end

    % simplex method algorithm, finding most negative elements in f
    % objective row
    while any(tableau(end,1:end-1) < 0)
       [tableau,basis] = pivotter(tableau,basis);
    end 

    % finding x_opt values from basis and tableau
    non_zero_opt_var = [basis', tableau(1:end-1,end)];
    non_zero_opt_x_var = non_zero_opt_var(ismember(non_zero_opt_var(:,1),1:A_col),:);
    x_opt(non_zero_opt_x_var(:,1)') =  non_zero_opt_x_var(:,2)';

    % finding z_opt
    z_opt = tableau(end,end);

else
    fprintf('Optimal Solution is not found.\nw ~= 0.\n')
end

end

