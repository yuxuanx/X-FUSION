function [ Error ] = gospa_dist( X,Y,c,p,alpha )

%Calculate sizes of the input point patterns
n = size(X,2);
m = size(Y,2);

if isempty(Y)
    dist_error = 0;
    miss_error = n*c^p/2;
    false_error= 0;
    error = n*c^p/2;
    Error = [error dist_error miss_error false_error];
    return;
end

%Calculate cost/weight matrix for pairings - fast method with vectorization
XX= repmat(X,[1 m]);
YY= reshape(repmat(Y,[n 1]),[size(Y,1) n*m]);
D = reshape(sqrt(sum((XX-YY).^2)),[n m]);

block = c^p/2*ones(n);
D_hat = [D block];

[assignment,~]= Hungarian(D_hat);

idx1 = assignment(:,1:m)==1;
false_est = m - sum(idx1(:));
idx2 = assignment(:,m+1:end)==1;
miss_est = sum(idx2(:));

false_error = false_est*c^p/alpha;
miss_error = miss_est*c^p/alpha;

temp = D.*assignment(:,1:m);
dist_error = sum(temp(:));

error = dist_error + miss_error + false_error;

Error = [error dist_error miss_error false_error];

end

