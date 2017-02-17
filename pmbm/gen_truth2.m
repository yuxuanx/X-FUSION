function [ truth ] = gen_truth2( model )

truth.K = 101;
truth.X = cell(truth.K,1);
truth.N = zeros(truth.K,1);
truth.total_tracks= 0;

midpoint = (truth.K+1)/2; 
numfb = midpoint-1;

Pmid = 1e-6*eye(model.x_dim); % Initial mid point covariance
numTruth = 6; % number of targets
birthtime = zeros(numTruth,1);

x = chol(Pmid)'*randn(size(model.F,1),numTruth);
xf = x; 
xb = x;
truth.X{midpoint} = x;
truth.N(midpoint) = size(x,2);

for t = 1:numfb
  % Run forward and backward simulation process
  xf = model.F*xf + model.Qc*randn(2,size(x,2));
  xb = model.F\(xb + model.Qc*randn(2,size(x,2)));
  truth.X{midpoint-t} = xb(:,midpoint-t>birthtime);
  truth.N(midpoint-t) = sum(midpoint-t>birthtime);
  truth.X{midpoint+t} = xf; % note that all targets exist after midpoint
  truth.N(midpoint+t) = size(xf,2);
end

truth.total_tracks = numTruth;

end


