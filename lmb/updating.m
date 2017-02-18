function glmb_update= updating(glmb_predict,model,filter,meas,k)
%gating by tracks
m_tracks= [];
P_tracks= [];
for tabidx=1:length(glmb_predict.tt)
    m_tracks= cat(2,m_tracks,glmb_predict.tt{tabidx}.m);
    P_tracks= cat(3,P_tracks,glmb_predict.tt{tabidx}.P);
end
meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_tracks,P_tracks);
%create updated tracks (single target Bayes update)
m= size(meas.Z{k},2);
tt_update= cell((1+m)*length(glmb_predict.tt),1);
%missed detection tracks (legacy tracks)
for tabidx= 1:length(glmb_predict.tt)
    tt_update{tabidx}= glmb_predict.tt{tabidx};
end

%measurement updated tracks (all pairs)
allcostm= zeros(length(glmb_predict.tt),m);
for emm= 1:m
    for tabidx= 1:length(glmb_predict.tt)
        stoidx= length(glmb_predict.tt)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted_tracks*j + i)
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k}(:,emm),model,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P);
        w_temp= qz_temp.*glmb_predict.tt{tabidx}.w+eps;
        tt_update{stoidx}.m= m_temp;
        tt_update{stoidx}.P= P_temp;
        tt_update{stoidx}.w= w_temp/sum(w_temp);
        tt_update{stoidx}.l = glmb_predict.tt{tabidx}.l;
        allcostm(tabidx,emm)= sum(w_temp);
    end
end
glmb_update.tt= tt_update;

%component updates
if m==0 %no measurements means all missed detections
    glmb_update.w= -model.lambda_c+glmb_predict.n*log(model.Q_D)+log(glmb_predict.w);
    glmb_update.I= glmb_predict.I;
    glmb_update.n= glmb_predict.n;
else %loop over predicted components/hypotheses
    runidx= 1;
    for pidx=1:length(glmb_predict.w)
        if glmb_predict.n(pidx)==0 %no target means all clutter
            glmb_update.w(runidx)= -model.lambda_c+m*log(model.lambda_c*model.pdf_c)+log(glmb_predict.w(pidx));             %weight of hypothesis/component
            glmb_update.I{runidx}= glmb_predict.I{pidx};                                                                    %tracks in hypothesis/component
            glmb_update.n(runidx)= glmb_predict.n(pidx);                                                                    %cardinality of hypothesis/component
            runidx= runidx+1;
        else %otherwise perform update for component
            %calculate best updated hypotheses/components
            costm= model.P_D/model.Q_D*allcostm(glmb_predict.I{pidx},:)/(model.lambda_c*model.pdf_c);                                           %cost matrix
            neglogcostm= -log(costm);                                                                                                           %negative log cost
            [uasses,nlcost]= mbestwrap_updt_custom(neglogcostm,ceil(filter.H_upd*glmb_predict.w(pidx)));    	%murty's algo to calculate m-best assignment hypotheses/components
            
            %generate corrresponding surviving hypotheses/components
            for hidx=1:length(nlcost)
                update_hypcmp_tmp= uasses(hidx,:)';
                glmb_update.w(runidx)= -model.lambda_c+m*log(model.lambda_c*model.pdf_c)+glmb_predict.n(pidx)*log(model.Q_D)+log(glmb_predict.w(pidx))-nlcost(hidx);    %weight of hypothesis/component
                glmb_update.I{runidx}= length(glmb_predict.tt).*update_hypcmp_tmp+glmb_predict.I{pidx};                                                                 %tracks in hypothesis/component
                glmb_update.n(runidx)= glmb_predict.n(pidx);                                                                                                            %cardinality of hypothesis/component
                runidx= runidx+1;
            end
        end
    end
end
glmb_update.w= exp(glmb_update.w-logsumexp(glmb_update.w));                                                                                                             %normalize weights

%extract cardinality distribution
for card=0:max(glmb_update.n)
    glmb_update.cdn(card+1)= sum(glmb_update.w(glmb_update.n==card));                                                                                                   %extract probability of n targets
end
end

