function [ x_est ] = state_extract(r,x)

r(r==1) = 1-eps;
ss = false(size(r));
pcard = prod(1-r)*poly(-r./(1-r));
[~,n] = max(pcard);
[~,o] = sort(-r);
n = n - 1;
ss(o(1:n)) = true;
x_est = x(:,ss);

end

