function [ lambdau,xu,Pu,r_update,x_update,p_update ] = ...
    recycling_v2( w_update,r_update,x_update,p_update,lambdau,xu,Pu,threshold )
% Recycling using variational approximation

len = length(w_update);
r_recycle = cell(len,1);
x_recycle = cell(len,1);
p_recycle = cell(len,1);
M = zeros(len,1);
for i = 1:len
    idx = r_update{i} < threshold;
    M(i) = length(idx);
    r_recycle{i} = r_update{i}(idx);
    x_recycle{i} = x_update{i}(:,idx);
    p_recycle{i} = p_update{i}(:,:,idx);
    idx = r_update{i} > threshold;
    r_update{i} = r_update{i}(idx);
    x_update{i} = x_update{i}(:,idx);
    p_update{i} = p_update{i}(:,:,idx);
end


h_r = zeros(0,1);
h_x = zeros(4,0);
h_p = zeros(4,4,0);
for i = 1:len
    h_r = [h_r;r_recycle{i}];
    h_x = [h_x x_recycle{i}];
    h_p = cat(3,h_p,p_recycle{i});
end
[~,idx_unique,~] = unique(h_x(1,:));
h_r = h_r(idx_unique)';
h_x = h_x(:,idx_unique);
h_p = h_p(:,:,idx_unique);

if isempty(h_r)
    return;
end

maxLen = max(M);
H = length(h_r);
phi = zeros(H,maxLen);

for i = 1:H
    for j = 1:len
        for k = 1:length(r_recycle{j})
            if isequal(h_x(1,i),x_recycle{j}(1,k)) && isequal(h_r(i),r_recycle{j}(k))
                phi(i,k) = phi(i,k) + w_update(j);
            end
        end
    end
end

ph = sum(phi,2);
pn = sum(phi,1)';

[C,r_hat,x_hat,P_hat] = cost(phi,h_r,h_x,h_p);

[Cmin,phi] = LP_transport(C,pn,ph);
temp = Cmin;
maxIterations = 1e2;
iteration = 0;
while(1)
    [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p);
    [Cmin,phi] = LP_transport(C,pn,ph);
    iteration = iteration + 1;
    if (temp - Cmin < threshold && temp >= Cmin) || (iteration > maxIterations)
        r_hat = r_temp;
        x_hat = x_temp;
        P_hat = P_temp;
        break;
    else
        temp = Cmin;
    end
end

lambdau = [lambdau;r_hat];
xu = [xu x_hat];
Pu = cat(3,Pu,P_hat);


end

