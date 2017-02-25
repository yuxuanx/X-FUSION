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

% Interpret length of inputs
len = length(r);

for l = 1:len
    n = length(r{l});
    % Implement prediction algorithm
    
    % Predict existing tracks
    r{l} = Ps*r{l};
    x{l} = F*x{l};
    for i = 1:n
        P{l}(:,:,i) = F*P{l}(:,:,i)*F' + Q;
    end
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



