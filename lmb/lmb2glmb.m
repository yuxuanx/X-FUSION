function glmb= lmb2glmb(tt_lmb,H_req)
%express LMB in GLMB structure
rvect= get_rvals(tt_lmb);                                   %vector of existence probabilities
costv= rvect./(1-rvect);                                    %cost vector
neglogcostv= -log(costv);                                   %negative log cost
[paths,nlcost]= kshortestwrap_pred(neglogcostv,H_req);      %k-shortest path to calculate k-best surviving hypotheses/components
glmb.tt= tt_lmb;
for hidx=1:length(nlcost)
    glmb.w(hidx)= sum(log(1-rvect))-nlcost(hidx);           %weight of hypothesis/component
    glmb.I{hidx}= paths{hidx}(:);                           %tracks in hypothesis/component
    glmb.n(hidx)= length(paths{hidx}(:));                   %cardinality of hypothesis/component
end
glmb.w= exp(glmb.w-logsumexp(glmb.w));                      %normalize weights

%extract cardinality distribution
for card=0:max(glmb.n)
    glmb.cdn(card+1)= sum(glmb.w(glmb.n==card));            %extract probability of n targets
end
end

