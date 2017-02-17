function [ rr, xx, PP ] = make_assign(bestAssign, rupd, xupd, Pupd, ...
    rnew, xnew, Pnew, n, m ,model)
%Make measurement to track assignment
dim = model.x_dim;
M = size(bestAssign,1);
rr = zeros(n+m,M);
xx = zeros(dim,n+m,M);
PP = zeros(dim,dim,n+m,M);

for k = 1:M
    r = zeros(n+m,1);
    x = zeros(dim,n+m);
    P = zeros(dim,dim,n+m);
    for i = 1:m
        if bestAssign(k,i) <= n
            r(bestAssign(k,i)) = 1-eps;
            x(:,bestAssign(k,i)) = xupd(:,bestAssign(k,i),i+1);
            P(:,:,bestAssign(k,i)) = Pupd(:,:,bestAssign(k,i),i+1);
        else
            r(i+n) = rnew(bestAssign(k,i)-n);
            x(:,i+n) = xnew(:,bestAssign(k,i)-n);
            P(:,:,i+n) = Pnew(:,:,bestAssign(k,i)-n);
        end
    end
    
    for i = 1:n
        if r(i) == 0
            r(i) = rupd(i,1);
            x(:,i) = xupd(:,i,1);
            P(:,:,i) = Pupd(:,:,i,1);
        end
    end
    
    rr(:,k) = r;
    xx(:,:,k) = x;
    PP(:,:,:,k) = P;
end

end

            