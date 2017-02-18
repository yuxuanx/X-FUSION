function est = run_filter(model,meas)


%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.T_max= 100;                  %maximum number of tracks
filter.track_threshold= 1e-4;       %threshold to prune tracks

filter.H_bth= 10;                    %requested number of birth components/hypotheses (for LMB to GLMB casting before update)
filter.H_sur= 50;                  %requested number of surviving components/hypotheses (for LMB to GLMB casting before update)
filter.H_upd= 100;                  %requested number of updated components/hypotheses (for GLMB update)
filter.H_max= 100;                  %cap on number of posterior components/hypotheses (not used yet)
filter.hyp_threshold= 1e-4;        %pruning threshold for components/hypotheses (not used yet)

filter.L_max= 10;                   %limit on number of Gaussians in each track
filter.elim_threshold= 1e-4;        %pruning threshold for Gaussians in each track
filter.merge_threshold= 1;          %merging threshold for Gaussians in each track

filter.P_G= 0.999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value

est.filter= filter;

%=== Filtering

%initial prior
tt_lmb_update= cell(0,1);      %track table for LMB (cell array of structs for individual tracks)

%recursive filtering
for k=1:meas.K
    
    %prediction
    [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update,model,k);
    
    %update
    glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter);
    glmb_update= updating(glmb_predict,model,filter,meas,k);
    
    %pruning and truncation
    tt_lmb_update= glmb2lmb(glmb_update);
    
    %pruning, truncation and track cleanup
    tt_lmb_update= clean_lmb(tt_lmb_update,filter);
    
    %state estimation
    [est.X{k},est.N(k),est.L{k}]= extract_estimates(tt_lmb_update,model);
    
end
end
