function [tableau, basis] = pivotter(tableau,basis)

% This function takes tableu and basis as inputs.
% Function rearrange tableu elements and basis based on pivotting rules of
% simplex method. 

[row,~] = size(tableau);

[~,pivot_col_idx] = min(tableau(end,1:end-1));

% finding RHS/coef rate for pivot column
rhs_coef_rate = tableau(:,end) ./ tableau(:,pivot_col_idx);
% force non positive ratio to be inf 
rhs_coef_rate(rhs_coef_rate<=0) = inf;  
% search through const Rows, then take min RHS/coef rate row
[~,pivot_row_idx] = min(rhs_coef_rate(1:end-1)); 

 % assign old basis var with new basis var
basis(pivot_row_idx) = pivot_col_idx;

% pivot row update through forcing row basis element coef as 1
pivotting_factor = tableau(pivot_row_idx,pivot_col_idx);
tableau(pivot_row_idx,:) = tableau(pivot_row_idx,:) / pivotting_factor;  % entering row update

% gauss elimination method for making non pivot elements in column as 0
for i=1:row
    if i ~= pivot_row_idx
        row_factor = tableau(i,pivot_col_idx) / tableau(pivot_row_idx,pivot_col_idx);
        tableau(i,:) = tableau(i,:) - row_factor * tableau(pivot_row_idx,:);
    end
end

end

