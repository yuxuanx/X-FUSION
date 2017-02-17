function [lambdau,xu,Pu,r_update,x_update,p_update,x_est,w_update] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model,w_update,IF_recycle)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Extract parameters from model
Pd = model.Pd;
H_max = model.H_max;
H_threshold = model.H_threshold;

% Only PPP update at the first time step
if isempty(r{1})
    [~,rr,xx,pp] = ppp_update(lambdau,xu,Pu,z,model);
    lambdau = lambdau*(1-Pd);
    
    if IF_recycle
        idx = rr > model.threshold;
        r_update{1} = rr(idx);
        x_update{1} = xx(:,idx);
        p_update{1} = pp(:,:,idx);
        x_est = x_update{1}(:,r_update{1}>0.5);
        
        idx = (rr < model.threshold) & (rr > H_threshold);
        lambdau = [lambdau;rr(idx)];
        xu = [xu xx(:,idx)];
        Pu = cat(3,Pu,pp(:,:,idx));
    else
        idx = rr > model.H_threshold;
        r_update{1} = rr(idx);
        x_update{1} = xx(:,idx);
        p_update{1} = pp(:,:,idx);
        x_est = x_update{1}(:,r_update{1}>0.5);
    end
    
    return;
end

% Extract number of global hypothesis
len = length(r);

% Allocate memory for existing tracks
r_update = cell(0,1);
x_update = cell(0,1);
p_update = cell(0,1);

% Update unknown targets
[w_new,r_new,x_new,P_new] = ppp_update(lambdau,xu,Pu,z,model);
lambdau = (1-Pd)*lambdau;

% Loop through the global hypothesis
w = cell(len,1);
for l = 1:len
    
    % Gating
    [valid_idx,idx_out] = gate_meas_gms(z,model,x{l},P{l});
    
    % Update existing tracks
    [wupd,rupd,xupd,Pupd] = mbm_update(z(:,valid_idx),r{l},x{l},P{l},model);
    
    % Calculate cost function
    cost = -log(wupd(:,2:end)./repmat(w_new(valid_idx)',length(r{l}),1));
    
    % Find M-best assignment
    Mt = ceil(H_max*w_update(l));
    [bestAssign, nCost] = mbestwrap_updt_custom(cost,Mt,wupd(:,1));
    
    % Update single target hypothesis
    [rMurty,xMurty,pMurty] = hypo_update(bestAssign,rupd,xupd,Pupd,...
        r_new,x_new,P_new,valid_idx,idx_out,nCost,H_threshold);
    
    w{l} = exp(-nCost)'*w_update(l);
    r_update = cat(1,r_update,rMurty);
    x_update = cat(1,x_update,xMurty);
    p_update = cat(1,p_update,pMurty);
end

% Normalisation
w_update = normalize(cell2mat(w));

% Pruning
[w_update,r_update,x_update,p_update] = pruning(w_update,r_update,...
    x_update,p_update,H_threshold);

% Capping
[w_update,r_update,x_update,p_update] = capping(w_update,r_update,...
    x_update,p_update,H_max);

% Recycling
if IF_recycle
    [lambdau,xu,Pu,r_update,x_update,p_update] = recycling_v2(w_update,...
        r_update,x_update,p_update,lambdau,xu,Pu,model.threshold);
end

% Best state extraction
x_est = state_extract( w_update, r_update, x_update);