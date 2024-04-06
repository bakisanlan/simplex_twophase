function tableau = RHS_posneg(tableau)

%checking are there any negative RHS column element
if any(tableau(:,end) < 0)
    neg_RHS_row = find(tableau(:,end-1) < 0);

    % select arbitrary positive RHS row for adding to negative RHS row
    pos_RHS_row = find(tableau(:,end-1) > 0);
    pos_RHS_row = pos_RHS_row(1);

    for k=neg_RHS_row
        adding_factor = tableau(k,end) / tableau(pos_RHS_row,end);
        tableau(k,:) = tableau(k,:) - adding_factor * tableau(pos_RHS_row,:);
    end

else 
    tableau = tableau;
end
end

