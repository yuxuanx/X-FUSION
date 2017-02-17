function [ C,r_hat,x_hat,P_hat ] = Mstep( qmiss,rupd,xupd,Pupd )

[H,J] = size(qmiss);
r_hat = zeros(J,1);
x_hat = zeros(4,J);
P_hat = zeros(4,4,J);
for j = 1:J
    r_hat(j) = rupd(j,:)*qmiss(:,j);
    for h = 1:H
        x_hat(:,j) = x_hat(:,j) + qmiss(h,j)*rupd(j,h)*xupd(:,j,h);
    end
    x_hat(:,j) = x_hat(:,j)/r_hat(j);
    for h = 1:H
        v = xupd(:,j,h) - x_hat(:,j);
        P_hat(:,:,j) = P_hat(:,:,j) + qmiss(h,j)*rupd(j,h)*(Pupd(:,:,j,h) + v*v');
    end
    P_hat(:,:,j) = P_hat(:,:,j)/r_hat(j);
end

C = zeros(H,J);
r_hat(r_hat==1) = 1-eps;
r_hat(r_hat==0) = eps;
for h = 1:H
    for j = 1:J
        v = xupd(:,j,h) - x_hat(:,j);
        inv_P = inv(P_hat(:,:,j));
        C(h,j) = -(1-rupd(j,h))*log(1-r_hat(j)) - rupd(j,h)*log(r_hat(j)) + ...
            rupd(j,h)/2*(trace(inv_P*Pupd(:,:,j,h)) + v'*inv_P*v + log(det(2*pi*P_hat(:,:,j))));
    end
end


end

