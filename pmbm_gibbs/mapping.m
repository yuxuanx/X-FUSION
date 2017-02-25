function [ w,r,x,P ] = mapping( assignments,wupd,rupd,xupd,Pupd,...
    w_new,r_new,x_new,P_new,valid_idx,idx_out,threshold )

wnew = w_new(valid_idx);
rnew = r_new(valid_idx);
xnew = x_new(:,valid_idx);
Pnew = P_new(:,:,valid_idx);
M = size(assignments,1);
r = cell(M,1);
x = cell(M,1);
P = cell(M,1);
w = ones(M,1)*prod(w_new(idx_out));
m = length(valid_idx);
n = size(wupd,1);
N = n+m;

for i = 1:M
    rr = zeros(N,1);
    xx = zeros(4,N);
    PP = zeros(4,4,N);
    for j = 1:m
        if assignments(i,j)>0
            rr(assignments(i,j)) = 1;
            xx(:,assignments(i,j)) = xupd(:,assignments(i,j),j+1);
            PP(:,:,assignments(i,j)) = Pupd(:,:,assignments(i,j),j+1);
            w(i) = w(i)*wupd(assignments(i,j),j+1);
        else
            rr(n+j) = rnew(j);
            xx(:,n+j) = xnew(:,j);
            PP(:,:,n+j) = Pnew(:,:,j);
            w(i) = w(i)*wnew(j);
        end
    end
    for k = 1:n
        if rr(k)==0
            rr(k) = rupd(k,1);
            xx(:,k) = xupd(:,k,1);
            PP(:,:,k) = Pupd(:,:,k,1);
            w(i) = w(i)*wupd(k,1);
        end
    end
    rc = [rr;r_new(idx_out)];
    xc = [xx x_new(:,idx_out)];
    Pc = cat(3,PP,P_new(:,:,idx_out));
    idx = rc > threshold;
    r{i} = rc(idx);
    x{i} = xc(:,idx);
    P{i} = Pc(:,:,idx);
end


end

