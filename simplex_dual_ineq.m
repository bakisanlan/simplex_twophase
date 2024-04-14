function [y_opt, z_dual_opt] = simplex_dual_ineq(f, A, b, n_gt, n_lt)
    % Generate Dual LP Problem
    
    % finding new n_gt and n_lt
    n_dual_const = length(f);
    if n_dual_const < n_gt + n_lt
        n_lt = min(n_dual_const,n_gt);  % n_gt turns to n_lt
        n_gt = n_dual_const - n_lt;
    end

    A_dual = A';  % Transpose of Aeq for dual problem [lt;gt]
    f_dual = -b';  % Negative of the primal objective function for dual problem
    b_dual = f'; % converting dual problem to maximize

    % forcing [lt;gt] --> [gt;lt]
    gt = A_dual(end-n_gt+1:end,:);
    lt = A_dual(1:n_lt,:);
    A_dual = [gt ; lt];

    gt_b = b_dual(end-n_gt+1:end);
    lt_b = b_dual(1:n_lt,:);
    b_dual = [gt_b ; lt_b];
    
    % forcing RHS to be positive
    for i=1:length(b_dual)
        if b_dual(i) < 0
            A_dual(i,:) = -A_dual(i,:); 
            b_dual(i) = -b_dual(i);


            % rearrenging rows as [gt ; lt]
            temp_A_row = A_dual(i,:);
            A_dual(i,:) = [];

            temp_b_row = b_dual(i);
            b_dual(i) = [];
            % lt to gt
            if i > n_gt
                A_dual = [temp_A_row ; A_dual];
                b_dual = [temp_b_row ; b_dual];
                n_gt = n_gt + 1;
                n_lt = n_lt - 1;            
            % gt to lt
            else
                A_dual = [A_dual ; temp_A_row];
                b_dual = [b_dual ; temp_b_row];
                n_gt = n_gt - 1;
                n_lt = n_lt + 1;
            end
        end
    end



    % Transform Dual LP Problem to Standard LP Form
    % m = size(A_dual, 1);
    % slack_vars = eye(m);  % Slack variables for converting inequalities to equalities
    % A_dual = [A_dual, slack_vars];  % Append slack variables to A_dual
    % %tableau_dual = [A_dual, b'; f_dual, 0];  % Create tableau for dual LP

    % Solve Dual LP Problem using Simplex Method
    [y_opt, z_dual_opt] = simplex_ineq(f_dual, A_dual, b_dual, n_gt,n_lt);

    % Output the optimal solution and objective value of the dual LP
end