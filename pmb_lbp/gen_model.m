function model= gen_model(pd,lambda)

% basic parameters
model.x_dim= 4;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector

% dynamical model parameters (CV model)
model.T= 1;                                     %sampling period
model.A0= [ 1 model.T; 0 1 ];                         %transition matrix                     
model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];
model.B0= [ (model.T^2)/2; model.T ];
model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
model.sigma_v = 5;
model.Q= (model.sigma_v)^2* model.B*model.B';   %process noise covariance

% survival/death parameters
model.Ps= 0.99;

% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ 10; 10 ]); 
model.R= model.D*model.D';              %observation noise covariance

% detection parameters
model.Pd= pd;   %probability of detection in measurements

% clutter parameters
model.lfai= lambda;                             %poisson average rate of uniform clutter (per scan)
model.range_c= [ -1000 1000; -1000 1000 ];      %uniform clutter region
model.lambda_fa= model.lfai/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

model.threshold = 1e-1; % Threshold for pruning low weights track
model.M = 100; % number of best assignments
model.H_threshold = 1e-4; % Pruning threshold

% Initialise new target parameter structure
birthNum = 4;
model.xb = zeros(4,birthNum);
model.xb(:,1) = [ 0; 0; 0; 0 ];
model.xb(:,2) = [ 400; 0; -600; 0 ];
model.xb(:,3) = [ -800; 0; -200; 0 ]; 
model.xb(:,4) = [ -200; 0; 800; 0 ]; 
model.Pb = zeros(4,4,birthNum);
for i = 1:birthNum
    model.Pb(:,:,i) = diag([ 10; 10; 10; 10 ])*diag([ 10; 10; 10; 10 ])';
end

model.lambdab = 0.1*ones(birthNum,1);

% Gating parameters
P_G= 0.999;                    %gate size in percentage
model.gamma= chi2inv(P_G,model.z_dim);   %inv chi^2 dn gamma value

