function [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update,model,k)
%---generate birth tracks
tt_lmb_birth= cell(length(model.r_birth),1);                                           %initialize cell array
for tabbidx=1:length(model.r_birth)
    tt_lmb_birth{tabbidx}.r= model.r_birth(tabbidx);                                   %birth prob for birth track
    tt_lmb_birth{tabbidx}.m= model.m_birth{tabbidx};                                   %means of Gaussians for birth track
    tt_lmb_birth{tabbidx}.P= model.P_birth{tabbidx};                                   %covs of Gaussians for birth track
    tt_lmb_birth{tabbidx}.w= model.w_birth{tabbidx}(:);                                %weights of Gaussians for birth track
    tt_lmb_birth{tabbidx}.l= [k;tabbidx];                                              %track label
end

%---generate surviving tracks
tt_lmb_survive= cell(length(tt_lmb_update),1);                                                                              %initialize cell array
for tabsidx=1:length(tt_lmb_update)
    tt_lmb_survive{tabsidx}.r= model.P_S*tt_lmb_update{tabsidx}.r;                                                          %predicted existence probability for surviving track
    [mtemp_predict,Ptemp_predict]= kalman_predict_multiple(model,tt_lmb_update{tabsidx}.m,tt_lmb_update{tabsidx}.P);        %kalman prediction
    tt_lmb_survive{tabsidx}.m= mtemp_predict;                                                                               %means of Gaussians for surviving track
    tt_lmb_survive{tabsidx}.P= Ptemp_predict;                                                                               %covs of Gaussians for predicted track
    tt_lmb_survive{tabsidx}.w= tt_lmb_update{tabsidx}.w;                                                                    %weights of Gaussians for predicted track
    tt_lmb_survive{tabsidx}.l= tt_lmb_update{tabsidx}.l;                                                                    %track label
end
end

