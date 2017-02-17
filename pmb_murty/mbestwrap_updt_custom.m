function [assignments,costs]= mbestwrap_updt_custom(P0,m,w)

if m==0
    assignments= [];
    costs= [];
    return;
end

l = length(w);
temp = ones(l)*Inf;
temp(logical(eye(l))) = 0;
blk1 = -log(diag(w)) + temp;

P0 = [P0' blk1];

% Make costs non-negative (required by 'assignmentoptimal')
x = min(min(P0));
P0 = P0 - x;

% Murty
[assignments, costs] = murty_custom(P0,m);

% Restore correct costs to assignments
costs = costs + (x.*sum(assignments>0,2))';


end
