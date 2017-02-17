function [ glmb_update ] = joint_prediction_update( glmb_update,model,filter,meas,k )

%get birth and survive table length
len_birth = length(model.r_birth);
len_existing = length(glmb_update.tt);

%create birth tracks
tt_birth= cell(len_birth,1); %initialize cell array
for tabbidx=1:len_birth
    tt_birth{tabbidx}.m= model.m_birth{tabbidx}; %means of Gaussians for birth track
    tt_birth{tabbidx}.P= model.P_birth{tabbidx}; %covs of Gaussians for birth track
end

%create surviving tracks - via time prediction (single target CK)
tt_survive= cell(len_existing,1); %initialize cell array
for tabsidx=1:len_existing
    %kalman prediction
    [tt_survive{tabsidx}.m,tt_survive{tabsidx}.P]= ...
        kalman_predict_multiple(model,glmb_update.tt{tabsidx}.m,glmb_update.tt{tabsidx}.P);
end

glmb_predict.tt= cat(1,tt_birth,tt_survive);

%gating by tracks
m_tracks= [];
P_tracks= [];
for tabidx=1:length(glmb_predict.tt)
    m_tracks= cat(2,m_tracks,glmb_predict.tt{tabidx}.m);
    P_tracks= cat(3,P_tracks,glmb_predict.tt{tabidx}.P);
end
meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_tracks,P_tracks);
m= size(meas.Z{k},2); %number of measurements

tt_birth_update= cell((1+m)*len_birth,1); %initialize cell array
for tabidx= 1:len_birth
    tt_birth_update{tabidx}= tt_birth{tabidx}; %same track table for missed detection
end

tt_exist_update= cell((1+m)*len_existing,1);
for tabidx= 1:len_existing
    tt_exist_update{tabidx}= tt_survive{tabidx};
end

%calculate the cost matrix for detected tracks (both birth and preexisted)
cost_updated_birth = zeros(len_birth,m);
cost_updated_exist = zeros(len_existing,m);

for emm= 1:m
    for tabidx= 1:len_birth
        stoidx= len_birth*emm + tabidx;
        %kalman update for this track and this measurement
        [qz_temp,tt_birth_update{stoidx}.m,tt_birth_update{stoidx}.P] = ...
            kalman_update_multiple(meas.Z{k}(:,emm),model,tt_birth{tabidx}.m,tt_birth{tabidx}.P); 
        %predictive likelihood
        cost_updated_birth(tabidx,emm)= qz_temp*model.r_birth(tabidx)*model.P_D/(model.lambda_c*model.pdf_c); 
    end
    for tabidx= 1:len_existing
        stoidx= len_existing*emm + tabidx;
        [qz_temp,tt_exist_update{stoidx}.m,tt_exist_update{stoidx}.P] = ...
            kalman_update_multiple(meas.Z{k}(:,emm),model,tt_survive{tabidx}.m,tt_survive{tabidx}.P);
        cost_updated_exist(tabidx,emm)= qz_temp*model.P_S*model.P_D/(model.lambda_c*model.pdf_c);
    end
end
neglogcost_updated_birth = -log(cost_updated_birth);
neglogcost_updated_exist = -log(cost_updated_exist);

%calculate the cost matrix for undetected birth tracks
neglogcost_birth_undetected = -log(diag(model.r_birth*model.Q_D));
%calculate the cost matrix for unconfirmed birth tracks
neglogcost_birth_die = -log(diag(1-model.r_birth));

glmb_update.tt= cat(1,tt_birth_update,tt_exist_update);

%ranked assignments
%cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
idx = 1;
for h = 1:length(glmb_update.w)
    labels = glmb_update.I{h}; %labels for current hypothesis
    num_exist = glmb_update.n(h);
    %calculate the cost matrix for undetected preexisted tracks
    neglogcost_exist_undetected = -log(diag(model.P_S*model.Q_D*ones(num_exist,1)));
    %calculate the cost matrix for non-survival tracks
    neglogcost_exist_nonsurvival = -log(diag(model.Q_S*ones(num_exist,1)));

    neglogcost_existing = [neglogcost_updated_exist(labels,:) neglogcost_exist_undetected Inf(num_exist,len_birth)...
    neglogcost_exist_nonsurvival Inf(num_exist,len_birth)];
    neglogcost_birth = [neglogcost_updated_birth Inf(len_birth,num_exist) neglogcost_birth_undetected...
    Inf(len_birth,num_exist) neglogcost_birth_die];
    neglogcost_total = [neglogcost_existing;neglogcost_birth];
    
    %Murty
    [uasses,nCost] = mbestwrap_updt(neglogcost_total,ceil(filter.H_upd*glmb_update.w(h)));
    
    for i = 1:length(nCost)
        update.I{idx} = [];
        for j = 1:num_exist % existed tracks
            if uasses(i,j) == m+j % missed detection
                update.I{idx} = [update.I{idx};length(tt_birth_update)+labels(j)];
            elseif uasses(i,j) <= m % updated by measurement
                update.I{idx} = [update.I{idx};length(tt_birth_update)+len_existing*uasses(i,j)+labels(j)];
            end
        end
        for j = num_exist+1:num_exist+len_birth % new born tracks
            if uasses(i,j) == m+j
                update.I{idx} = [update.I{idx};j-num_exist];  
            elseif uasses(i,j) <= m
                update.I{idx} = [update.I{idx};len_birth*uasses(i,j)+j-num_exist];
            end
        end
        update.n(idx) = length(update.I{idx});
        update.w(idx) = exp(-nCost(i))*glmb_update.w(h);
        idx = idx + 1;
    end
end

glmb_update.w = update.w;
glmb_update.I = update.I;
glmb_update.n = update.n;

%remove duplicate entries and clean track table
glmb_update= clean_predict(glmb_update);
glmb_update= clean_update(glmb_update);

glmb_update.w = glmb_update.w/sum(glmb_update.w);
%extract cardinality distribution
for card=0:max(glmb_update.n)
    glmb_update.cdn(card+1)= sum(glmb_update.w(glmb_update.n==card));
end


end

