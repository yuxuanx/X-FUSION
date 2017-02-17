function [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model)
%PREDICT: PREDICT MULTI-BERNOULLI AND POISSON COMPONENTS

% Get multi-Bernoulli prediction parameters from model
F = model.F;
Q = model.Q;
Ps = model.Ps;

% Get birth parameters from model
lambdab = model.lambdab;
xb = model.xb;
Pb = model.Pb;

% Implement prediction algorithm
n = length(r);
% Predict existing tracks
r = Ps*r;
x = F*x;
for i = 1:n
    P(:,:,i) = F*P(:,:,i)*F' + Q;
end

% Predict existing PPP intensity
lambdau = Ps*lambdau;
xu = F*xu;
nu = size(xu,2);
for k = 1:nu
    Pu(:,:,k) = F*Pu(:,:,k)*F' + Q;
end

% Incorporate birth intensity into PPP
lambdau = [lambdau;lambdab];
xu = [xu xb];
Pu = cat(3,Pu,Pb);

% Truncate low weight components
ss = lambdau > model.H_threshold;
lambdau = lambdau(ss);
xu = xu(:,ss);
Pu = Pu(:,:,ss);





