function [assignments,costs]= mbestwrap_updt_custom(P0,m,w)

if m==0
    assignments= [];
    costs= [];
    return;
end

n = size(P0,2);

% Padding blocks for dummy variables
% blk1 = -log(repmat(w(:,1),1,n1));
l = length(w);
temp = ones(l)*Inf;
temp(logical(eye(l))) = 0;
blk1 = -log(diag(w)) + temp;
% blk1 = -log(diag(w(:,1))+eps);

P0 = [P0 blk1];

% Make costs non-negative (required by 'assignmentoptimal')
x = min(min(P0));
P0 = P0 - x;

% Murty
[assignments, costs] = murty_custom(P0,m);

% Restore correct costs to assignments
costs = costs + (x.*sum(assignments>0,2))';

assignments(assignments>n) = 0;

end

