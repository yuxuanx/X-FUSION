function glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter)
%express birth and surviving LMBs in GLMB structure
glmb_birth= lmb2glmb(tt_lmb_birth,filter.H_bth);
glmb_survive= lmb2glmb(tt_lmb_survive,filter.H_sur);

%generate predicted hypotheses/components (by convolution of birth and survive GLMBs)
glmb_predict.tt= cat(1,glmb_birth.tt,glmb_survive.tt);
for bidx= 1:length(glmb_birth.w)
    for sidx= 1:length(glmb_survive.w)
        hidx= (bidx-1)*length(glmb_survive.w)+sidx;
        glmb_predict.w(hidx)= glmb_birth.w(bidx)*glmb_survive.w(sidx);
        glmb_predict.I{hidx}= [glmb_birth.I{bidx}; length(glmb_birth.tt)+glmb_survive.I{sidx}];
        glmb_predict.n(hidx)= glmb_birth.n(bidx)+glmb_survive.n(sidx);
    end
end
glmb_predict.w= glmb_predict.w/sum(glmb_predict.w);

%extract cardinality distribution
for card=0:max(glmb_predict.n)
    glmb_predict.cdn(card+1)= sum(glmb_predict.w(glmb_predict.n==card));
end
end
