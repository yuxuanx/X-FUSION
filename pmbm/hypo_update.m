function [rMurty,xMurty,pMurty] = hypo_update( bestAssign,rupd,xupd,Pupd,...
    r_new,x_new,P_new,valid_idx,idx_out,nCost,threshold)
%Update single target hypothesis according to the assignment

M = length(nCost);
rMurty = cell(M,1);
xMurty = cell(M,1);
pMurty = cell(M,1);

for i = 1:M
    [rr,xx,PP] = make_assign(bestAssign(i,:),rupd,xupd,Pupd,...
        r_new(valid_idx),x_new(:,valid_idx),P_new(:,:,valid_idx));
    
    rr = cat(1,rr,r_new(idx_out));
    xx = cat(2,xx,x_new(:,idx_out));
    PP = cat(3,PP,P_new(:,:,idx_out));
    
    % Prune low weight tracks
    idx = rr >= threshold;
    rMurty{i} = rr(idx);
    xMurty{i} = xx(:,idx);
    pMurty{i} = PP(:,:,idx);
    
end

end

