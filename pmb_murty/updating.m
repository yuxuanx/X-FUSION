function [lambdau,xu,Pu,r_update,x_update,p_update,x_est] = ...
    updating(lambdau,xu,Pu,r,x,P,z,model,IF_recycle)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS

% Only PPP update at the first time step
if isempty(r)
    [~,rr,xx,PP] = ppp_update(lambdau,xu,Pu,z,model);
    lambdau = (1-model.Pd)*lambdau;
    if IF_recycle
        idx = rr > model.threshold;
        r_update = rr(idx);
        x_update = xx(:,idx);
        p_update = PP(:,:,idx);
        x_est = x_update(:,r_update>0.5);
        idx = (rr > model.H_threshold) & (rr < model.threshold);
        lambdau = [lambdau;rr(idx)];
        xu = [xu xx(:,idx)];
        Pu = cat(3,Pu,PP(:,:,idx));
    else
        idx = rr > model.H_threshold;
        r_update = rr(idx);
        x_update = xx(:,idx);
        p_update = PP(:,:,idx);
        x_est = x_update(:,r_update>0.5);
    end
    return;
end

% Gating
[z_gate,idx_out] = gate_meas_gms(z,model,x,P);

% Update unknown targets
[wnew,rnew,xnew,Pnew] = ppp_update(lambdau,xu,Pu,z_gate,model);
[~,rout,xout,Pout] = ppp_update(lambdau,xu,Pu,z(:,idx_out),model);
lambdau = (1-model.Pd)*lambdau;

% Update existing tracks
[wupd,rupd,xupd,Pupd] = mbm_update(z_gate,r,x,P,model);

% Calculate cost function
cost = -log(wupd(:,2:end)./repmat(wupd(:,1),1,size(z_gate,2)));

% Find M-best assignment
[bestAssign, nCost] = mbestwrap_updt_custom(cost,model.M,wnew);

% Update single target hypothesis
[r_update,x_update,p_update,lambdau,xu,Pu] = hypo_update(bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,model,nCost,lambdau,xu,Pu,IF_recycle);

% Best state extraction
x_est = state_extract(r_update,x_update);
