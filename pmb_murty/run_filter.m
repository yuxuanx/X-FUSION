function [ est ] = run_filter( model,meas,IF_recycle )

est.X= cell(meas.K,1);

% State dimension
dim = model.x_dim;
% Multi-Bernoulli representation
n = 0;
r = zeros(0,1);
x = zeros(dim,n);
P = zeros(dim,dim,n);

% Unknown target PPP parameters
lambdau = model.lambdab;
xu = model.xb;
Pu = model.Pb;

% Loop through time
for t = 1:meas.K
    % Predict
    [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
    
    % Update
    [lambdau,xu,Pu,r,x,P,est.X{t}] = updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model,IF_recycle);
    
end


end

