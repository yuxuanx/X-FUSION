function [r,x,P] = make_assign(bestAssign,rupd,xupd,Pupd,rnew,xnew,Pnew)
%Make measurement to track assignment

n = size(rupd,1);
m = length(rnew);
r = zeros(n+m,1);
x = zeros(4,n+m);
P = zeros(4,4,n+m);

for i=1:n
    if bestAssign(i) ~= 0
        r(i) = 1;
        x(:,i) = xupd(:,i,bestAssign(i)+1);
        P(:,:,i) = Pupd(:,:,i,bestAssign(i)+1);
    else
        r(i) = rupd(i,1);
        x(:,i) = xupd(:,i,1);
        P(:,:,i) = Pupd(:,:,i,1);
    end
end

for i=n+1:n+m
    if ~any(bestAssign(:)==i-n)
        r(i) = rnew(i-n);
        x(:,i) = xnew(:,i-n);
        P(:,:,i) = Pnew(:,:,i-n);
    end
end

end

