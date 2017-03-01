clc;clear
dbstop if error
% Generate model
model= gen_model(0.75,10);
load('truth_10_75');
load('meas_10_75');
IF_recycle = false;
% Monte Carlo simulations
numTrial = 50;
K = 100;
% GOSPA parameters
gospa_p= 1;
gospa_c= 100;
gospa_alpha= 2;
gospa_vals= zeros(K,4,numTrial);

time = zeros(numTrial,1);

parfor trial = 1:numTrial
    
%     truth= gen_truth(model);
%     meas=  gen_meas(model,truth);
    tic
    est = run_filter(model,meas{trial},IF_recycle);
    time(trial) = toc;
    
    for t = 1:K
        % Performance evaluation using GOSPA metric
        [gospa_vals(t,:,trial)] = gospa_dist(get_comps(truth{trial}.X{t},[1 3]),...
            get_comps(est.X{t},[1 3]),gospa_c,gospa_p,gospa_alpha);
    end
    
end

averGospa = mean(gospa_vals,3);
mean(averGospa)



