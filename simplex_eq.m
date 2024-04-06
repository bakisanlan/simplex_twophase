function [x_opt,z_opt] = simplex_eq(f,Aeq,beq)
% This function takes maximizing objective function as f, equality
% constraints matrix as Aeq, right hand side of constraints as beq.
% Then gives output as [x_opt,z_opt]

% ------------------------------------------------------------------
% Example use of simplex_eq function  
% Aeq= [1 -1 0 0 2 0 ; -2 1 0 0 -2 0 ; 1 0 1 0 1 -1 ; 0 1 1 1 2 1];
% beq = [0 ; 0 ; 3 ; 4];
% f = [40 10 0 0 7 14];
% [x_opt,z_opt] = simplex_eq(f,Aeq,beq)
% ------------------------------------------------------------------

[A_row, A_col] = size(Aeq);

n_eq = A_row;
%S_var_matrix = diag([-ones(1,n_gt), ones(1,n_lt)]); % S variable matrix

% A variable matrix
A_var_matrix = eye(n_eq);        

% fake objective function for maximazing of -(sum(artificial vars)) 
w_matrix = [-sum(Aeq,1) , zeros(1,n_eq)];

% RHS column
RHS = [beq ; -sum(beq)];

% tableau without basis and original cost obj. function
tableau = [[Aeq , A_var_matrix ; w_matrix] , RHS];  
basis = [A_col+1: A_col+A_row]; % basis index of each const Equations

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
if ~any(ismember(basis,A_col+1:A_col+A_row)) && tableau(end,end) < 10e-12 
    fprintf("Optimal Solution Found.\nw=0 and there are no basis artificial vars in objective row.\n")

    % replacing w with original f objective function
    % delete artificial var column
    tableau(:,A_col+1:A_col+A_row) = [];
    tableau(end,:) = [-f, 0];

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
    fprintf('Optimal Solution Found.\nw=0 and there is at least one basis artificial var in objective row.\n')

    % storing basis artificial var.
    art_var = A_col+1:A_col+A_row ;
    basis_art_var_idx = ismember(basis,art_var);
    basis_art_var = basis(basis_art_var_idx);

    % storing non basis artifical var.
    non_basis_art_var_idx = ~ismember(art_var,basis_art_var);
    non_basis_art_var = art_var(non_basis_art_var_idx);

    % dropping(replace with zeros) non basis artifical var columns
    tableau(:,non_basis_art_var) = zeros(A_row+1,length(non_basis_art_var));

    % dropping(replace with zeros) non basis original var with positive coef. in objective function 
    coef_of_org_var_in_obj = tableau(end,1:A_col);
    org_var_pos_coef_in_obj_logical= (coef_of_org_var_in_obj > 0);
    org_var_pos_coef_in_obj = find(org_var_pos_coef_in_obj_logical == 1);

    % replacing w with original f objective function
    tableau(end,:) = [-f, zeros(1,A_row+1)];

    % dropping(replace with zeros) non basis original var with positive coef. in objective function 

    tableau(:,org_var_pos_coef_in_obj) = zeros(A_row+1,length(org_var_pos_coef_in_obj));

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
    x_opt = [];
    z_opt = [];
end


end

