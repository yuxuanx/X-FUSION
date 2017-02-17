function [ r_update,x_update,p_update,lambdau,xu,Pu] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    rnew,xnew,Pnew,rout,xout,Pout,model,nCost,lambdau,xu,Pu,IF_recycle )
%Update single target hypothesis according to the assignment

n = size(rupd,1);
m = length(rnew);
% Making assignments
[rr,xx,PP] = make_assign(bestAssign,rupd,xupd,Pupd,rnew,...
    xnew,Pnew,n,m,model);
w = exp(-nCost-logsumexp(-nCost))';
if isempty(w)
    rr = cat(1,rr,rout);
    xx = cat(2,xx,xout);
    PP = cat(3,PP,Pout);
else
    % Generate q(h,j) r_h,x_h,P_h
    [pn,ph,phi,h_r,h_x,h_p] = hypo_all(w,rr,xx,PP,m+n,model);
    
    % Generate cost function
    [C,r_hat,x_hat,P_hat] = cost(phi,h_r,h_x,h_p,model);
    
    % Solve linear transport problem
    [Cmin,phi] = LP_transport(C,pn,ph);
    temp = Cmin;
    maxIterations = 1e2;
    iteration = 0;
    while(1)
        [C,r_temp,x_temp,P_temp] = cost(phi,h_r,h_x,h_p,model);
        [Cmin,phi] = LP_transport(C,pn,ph);
        iteration = iteration + 1;
        if (temp - Cmin < model.H_threshold && temp >= Cmin) || (iteration > maxIterations)
            r_hat = r_temp;
            x_hat = x_temp;
            P_hat = P_temp;
            break;
        else
            temp = Cmin;
        end
    end
    
    rr = cat(1,r_hat,rout);
    xx = cat(2,x_hat,xout);
    PP = cat(3,P_hat,Pout);
end

% Prune low track weights
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


end

