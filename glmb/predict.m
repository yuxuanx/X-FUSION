function glmb_predict= predict(glmb_update,model,filter)
%---generate birth hypotheses/components
%create birth tracks
tt_birth= cell(length(model.r_birth),1);           %initialize cell array
for tabbidx=1:length(model.r_birth)
    tt_birth{tabbidx}.m= model.m_birth{tabbidx};   %means of Gaussians for birth track
    tt_birth{tabbidx}.P= model.P_birth{tabbidx};   %covs of Gaussians for birth track
%     tt_birth{tabbidx}.w= model.w_birth{tabbidx}(:);%weights of Gaussians for birth track
%     tt_birth{tabbidx}.l= [k;tabbidx];              %track label
end
glmb_birth.tt= tt_birth;                           %copy track table back to GLMB struct

%calculate best birth hypotheses/components
costv= model.r_birth./(1-model.r_birth);           %cost vector
neglogcostv= -log(costv);                          %negative log cost
[bpaths,nlcost]= kshortestwrap_pred(neglogcostv,filter.H_bth);%k-shortest path to calculate k-best births hypotheses/components

%generate corrresponding birth hypotheses/components
for hidx=1:length(nlcost)
    birth_hypcmp_tmp= bpaths{hidx}(:);
    glmb_birth.w(hidx)= sum(log(1-model.r_birth))-nlcost(hidx);%hypothesis/component weight
    glmb_birth.I{hidx}= birth_hypcmp_tmp;                      %hypothesis/component tracks (via indices to track table)
    glmb_birth.n(hidx)= length(birth_hypcmp_tmp);              %hypothesis/component cardinality
end
glmb_birth.w= exp(glmb_birth.w-logsumexp(glmb_birth.w));       %normalize weights

%extract cardinality distribution
for card=0:max(glmb_birth.n)
    glmb_birth.cdn(card+1)= sum(glmb_birth.w(glmb_birth.n==card));%extract probability of n targets
end

%---generate survival hypotheses/components
%create surviving tracks - via time prediction (single target CK)
tt_survive= cell(length(glmb_update.tt),1);                                                                                 %initialize cell array
for tabsidx=1:length(glmb_update.tt)
    [mtemp_predict,Ptemp_predict]= kalman_predict_multiple(model,glmb_update.tt{tabsidx}.m,glmb_update.tt{tabsidx}.P);      %kalman prediction for GM
    tt_survive{tabsidx}.m= mtemp_predict;                                                                                   %means of Gaussians for surviving track
    tt_survive{tabsidx}.P= Ptemp_predict;                                                                                   %covs of Gaussians for surviving track
%     tt_survive{tabsidx}.w= glmb_update.tt{tabsidx}.w;                                                                       %weights of Gaussians for surviving track
%     tt_survive{tabsidx}.l= glmb_update.tt{tabsidx}.l;                                                                       %track label
end
glmb_survive.tt= tt_survive;                                                                                                %copy track table back to GLMB struct

%loop over posterior components/hypotheses
runidx= 1;                                                                                      %counter and index variable for new hypotheses/components
for pidx=1:length(glmb_update.w)
    if glmb_update.n(pidx)==0 %no target means no deaths
        glmb_survive.w(runidx)= log(glmb_update.w(pidx));       %hypothesis/component weight
        glmb_survive.I{runidx}= glmb_update.I{pidx};            %hypothesis/component tracks (via indices to track table)
        glmb_survive.n(runidx)= glmb_update.n(pidx);            %hypothesis component cardinality
        runidx= runidx+1;
    else %perform prediction for survivals
        %calculate best surviving hypotheses/components
        costv= model.P_S/model.Q_S*ones(glmb_update.n(pidx),1);                                                                     %cost vector
        neglogcostv= -log(costv);                                                                                                   %negative log cost
        [spaths,nlcost]= kshortestwrap_pred(neglogcostv,ceil(filter.H_sur*glmb_update.w(pidx)));    %k-shortest path to calculate k-best surviving hypotheses/components
        %generate corrresponding surviving hypotheses/components
        for hidx=1:length(nlcost)
            survive_hypcmp_tmp= spaths{hidx}(:);
            glmb_survive.w(runidx)= glmb_update.n(pidx)*log(model.Q_S)+log(glmb_update.w(pidx))-nlcost(hidx);                       %hypothesis/component weight
            glmb_survive.I{runidx}= glmb_update.I{pidx}(survive_hypcmp_tmp);                                                        %hypothesis/component tracks (via indices to track table)
            glmb_survive.n(runidx)= length(survive_hypcmp_tmp);                                                                     %hypothesis/component cardinality
            runidx= runidx+1;
        end
    end
end
glmb_survive.w= exp(glmb_survive.w-logsumexp(glmb_survive.w));                                                                      %normalize weights

%extract cardinality distribution
for card=0:max(glmb_survive.n)
    glmb_survive.cdn(card+1)= sum(glmb_survive.w(glmb_survive.n==card));                                                            %extract probability of n targets
end

%---generate predicted hypotheses/components (by convolution of birth and survive GLMBs)
%perform convolution - just multiplication
glmb_predict.tt= cat(1,glmb_birth.tt,glmb_survive.tt);                                                                              %concatenate track table
for bidx= 1:length(glmb_birth.w)
    for sidx= 1:length(glmb_survive.w)
        hidx= (bidx-1)*length(glmb_survive.w)+sidx;
        glmb_predict.w(hidx)= glmb_birth.w(bidx)*glmb_survive.w(sidx);                                                              %hypothesis/component weight
        glmb_predict.I{hidx}= [glmb_birth.I{bidx}; length(glmb_birth.tt)+glmb_survive.I{sidx}];                                     %hypothesis/component tracks (via indices to track table)
        glmb_predict.n(hidx)= glmb_birth.n(bidx)+glmb_survive.n(sidx);                                                              %hypothesis/component cardinality
    end
end

%remove duplicate entries and clean track table
glmb_predict= clean_predict(glmb_predict);

glmb_predict.w= glmb_predict.w/sum(glmb_predict.w);                                                                                 %normalize weights

%extract cardinality distribution
for card=0:max(glmb_predict.n)
    glmb_predict.cdn(card+1)= sum(glmb_predict.w(glmb_predict.n==card));                                                            %extract probability of n targets
end
end

