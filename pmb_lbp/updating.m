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
%         x_est = x_update(:,r_update>0.5);
        x_est = state_extract(r_update,x_update);
        idx = (rr > model.H_threshold) & (rr < model.threshold);
        lambdau = [lambdau;rr(idx)];
        xu = [xu xx(:,idx)];
        Pu = cat(3,Pu,PP(:,:,idx));
    else
        idx = rr > model.H_threshold;
        r_update = rr(idx);
        x_update = xx(:,idx);
        p_update = PP(:,:,idx);
%         x_est = x_update(:,r_update>0.5);
        x_est = state_extract(r_update,x_update);
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

% Loopy belief propagation
[pupd,pnew] = lbp(wupd,wnew);
ph = sum(pupd,1)';
qmiss = pupd';

% Generate cost function
[ C,r_hat,x_hat,P_hat ] = Mstep( qmiss,rupd,xupd,Pupd );

% Solve linear transport problem
[Cmin,qmiss] = LP_transport(C,ph);
temp = Cmin;
maxIterations = 1e2;
iteration = 0;
while(1)
    [ C,r_temp,x_temp,P_temp ] = Mstep( qmiss,rupd,xupd,Pupd );
    [Cmin,qmiss] = LP_transport(C,ph);
    iteration = iteration + 1;
    if (temp - Cmin < 1e-4 && temp >= Cmin) || (iteration > maxIterations)
        r_hat = r_temp;
        x_hat = x_temp;
        P_hat = P_temp;
        break;
    else
        temp = Cmin;
    end
end
rr = [r_hat;pnew.*rnew;rout];
xx = [x_hat xnew xout];
PP = cat(3,cat(3,P_hat,Pnew),Pout);

% Pruning and recycling
if IF_recycle
    idx = rr > model.threshold;
    r_update = rr(idx);
    x_update = xx(:,idx);
    p_update = PP(:,:,idx);
    idx = (rr > model.H_threshold) & (rr < model.threshold);
    lambdau = [lambdau;rr(idx)];
    xu = [xu xx(:,idx)];
    Pu = cat(3,Pu,PP(:,:,idx));
else
    idx = rr > model.H_threshold;
    r_update = rr(idx);
    x_update = xx(:,idx);
    p_update = PP(:,:,idx);
end

% Best states estimation
x_est = state_extract(r_update,x_update);
% x_est = state_extract2(r_update,x_update,0.5);

