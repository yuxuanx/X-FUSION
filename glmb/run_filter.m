function est = run_filter(model,meas)

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);

%filter parameters
filter.H_bth= 10;                    %requested number of birth components/hypotheses
filter.H_sur= 3000;                 %requested number of surviving components/hypotheses
filter.H_upd= 3000;                 %requested number of updated components/hypotheses
filter.H_max= 3000;                 %cap on number of posterior components/hypotheses
filter.hyp_threshold= 1e-5;        %pruning threshold for components/hypotheses

filter.P_G= 0.999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value

est.filter= filter;

%=== Filtering

%initial prior
glmb_update.tt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
glmb_update.w= 1;               %vector of GLMB component/hypothesis weights
glmb_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
glmb_update.n= 0;               %vector of GLMB component/hypothesis cardinalities
glmb_update.cdn= 1;             %cardinality distribution of GLMB (vector of cardinality distribution probabilities)

%recursive filtering
for k=1:meas.K
    
    %prediction and update
%     glmb_predict= predict(glmb_update,model,filter);
%     glmb_update= updating(glmb_predict,model,filter,meas,k);
    glmb_update = joint_prediction_update(glmb_update,model,filter,meas,k);
    
    %pruning and truncation
    glmb_update= prune(glmb_update,filter);
    glmb_update= cap(glmb_update,filter);
    
    %state estimation and display diagnostics
    [est.X{k},est.N(k)]= extract_estimates(glmb_update,model);
    
end
        
end


