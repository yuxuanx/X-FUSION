function [ lambdau,xu,Pu,r_update,x_update,p_update ] = ...
    recycling( w_update,r_update,x_update,p_update,lambdau,xu,Pu,threshold )

len = length(w_update);
for i = 1:len
    idx = r_update{i} < threshold;
    lambdau = [lambdau;w_update(i)*r_update{i}(idx)];
    xu = [xu x_update{i}(:,idx)];
    Pu = cat(3,Pu,p_update{i}(:,:,idx));
    idx = r_update{i} > threshold;
    r_update{i} = r_update{i}(idx);
    x_update{i} = x_update{i}(:,idx);
    p_update{i} = p_update{i}(:,:,idx);
end

% Merge components within close Mahalanobis distance
[lambdau,xu,Pu]= gaus_merge(lambdau,xu,Pu,1);


end

