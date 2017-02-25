function [ C,r_hat,x_hat,P_hat ] = cost(phi,h_r,h_x,h_p)

dim = 4;
[H,N] = size(phi);
x_hat = zeros(dim,N);
P_hat = zeros(dim,dim,N);
v = zeros(dim,H,N);

phi = phi + eps;
r_hat = (h_r*phi)';
r_hat(r_hat>1-eps) = 1-eps;

for j = 1:N
    x_hat(:,j) = sum(repmat(h_r.*phi(:,j)',dim,1).*h_x,2)/r_hat(j);
    for h = 1:H
        v(:,h,j) = h_x(:,h) - x_hat(:,j);
        P_hat(:,:,j) = P_hat(:,:,j) + phi(h,j)*h_r(h)*(h_p(:,:,h)+v(:,h,j)*v(:,h,j)');
    end
end
P_hat = abs(P_hat)./reshape(kron(r_hat,ones(dim))',dim,dim,N);

C = zeros(H,N);
for j = 1:N
    for h = 1:H
        temp = trace(P_hat(:,:,j)\h_p(:,:,h)) + v(:,h,j)'/P_hat(:,:,j)*v(:,h,j) + log(abs(det(2*pi*P_hat(:,:,j))));
        C(h,j) = -(1-h_r(h))*log(1-r_hat(j)) - h_r(h)*log(r_hat(j)) + h_r(h)/2*temp;
    end
end

end

