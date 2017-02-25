function est = run_filter(model,meas,IF_recycle)

est.X= cell(meas.K,1);

% Multi-Bernoulli representation
n = 0;
r = cell(0,1);
r{1} = zeros(n,1);
x = cell(0,1);
x{1} = zeros(4,n);
P = cell(0,1);
P{1} = zeros(4,4,n);
w_update = 1;

% Unknown target PPP parameters
lambdau = model.lambdab;
xu = model.xb;
Pu = model.Pb;

for t = 1:meas.K
    
    % Predict
    [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
    
    % Update
    [lambdau,xu,Pu,r,x,P,est.X{t},w_update] = ...
        updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model,w_update,IF_recycle);
    
end


end

