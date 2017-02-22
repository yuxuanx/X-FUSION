clc;clear
dbstop if error
% Generate model
model= gen_model(0.75,10);
% Recycling indicator
IF_recycle = true;
% Monte Carlo simulations
numTrial = 200;
K = 100;
% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(K,4,numTrial);

parfor trial = 1:numTrial
    
    % Generate ground truth
    truth= gen_truth(model);
    
    % Generate measurements
    meas=  gen_meas(model,truth);
    
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
    
    for t = 1:K
        
        % Predict
        [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model);
        
        % Update
        [lambdau,xu,Pu,r,x,P,x_est,w_update] = ...
            updating(lambdau,xu,Pu,r,x,P,meas.Z{t},model,w_update,IF_recycle);
        
        % GOSPA
        [gospa_vals(t,:,trial)] = ...
            gospa_dist(get_comps(truth.X{t},[1 3]),get_comps(x_est,[1 3])...
            ,gospa_c,gospa_p,gospa_alpha);
    end
end

averGospa = mean(gospa_vals,3);
mean(averGospa)


