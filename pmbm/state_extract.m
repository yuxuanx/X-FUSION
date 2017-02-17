function [ x_est ] = state_extract( w_update, r_update, x_update)

% extract number of MB components
num = length(w_update);
pcard = cell(num,1);
len = zeros(num,1);

% for each MB component, obtain cardinality pmf
for i = 1:num
    r = r_update{i};
    r(r==1) = 1-eps;
    pcard{i} = prod(1-r)*poly(-r./(1-r));
    len(i) = length(pcard{i});
end

% add 0 to make each MB component has equal length of cardinality pmf
len_max = max(len);
pcard_mb = zeros(num,len_max);
for i = 1:num
    if len(i) < len_max
        pcard_mb(i,:) = w_update(i)*[pcard{i} zeros(1,len_max-length(pcard{i}))];
    else
        pcard_mb(i,:) = w_update(i)*pcard{i};
    end
end

% the MBM cardinality pmf is equal to the sum of the MB weights multiplied 
% with the MB cardinality PMFs
pcard_mbm = sum(pcard_mb);

% MAP cardinality estimate
[~,n] = max(pcard_mbm);
C_max = n-1;

% for each MB hypothesis in the MBM, find the MAP cardinality estimate
C_mb = zeros(num,1);
for i = 1:num
    C_mb(i) = sum(r_update{i}>=0.5);
end

% take the highest weight MB component satisfy C_max = C_mb
w_update(C_mb~=C_max) = 0;
[~,idx] = max(w_update);

% state extraction
idices = r_update{idx} >= 0.5;
x_est = x_update{idx}(:,idices);

end