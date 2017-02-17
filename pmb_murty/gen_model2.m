function model= gen_model2(pd,lambda)

% basic parameters
model.x_dim= 4;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector

% dynamical model parameters (CV model)
model.T= 1;                                     %sampling period
model.A0= [ 1 model.T; 0 1 ];                         %transition matrix                     
model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];
model.B0= [ (model.T^2)/2; model.T ];
model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
model.sigma_v = 0.2;
model.Q= (model.sigma_v)^2* model.B*model.B';   %process noise covariance
model.Qc = model.sigma_v*model.B;

% survival/death parameters
model.Ps= .99;

% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ 1; 1 ]); 
model.R= model.D*model.D';              %observation noise covariance

% detection parameters
model.Pd= pd;   %probability of detection in measurements

% clutter parameters
model.lfai= lambda;                             %poisson average rate of uniform clutter (per scan)
model.range_c= [ -100 100; -100 100 ];      %uniform clutter region
model.lambda_fa= model.lfai/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

model.M = 100;
model.threshold = 1e-1; % Threshold for pruning low weights track
model.H_max = 100; % capping threshold
model.H_threshold = 1e-4; % Pruning threshold

% Initialise new target parameter structure

model.xb = zeros(4,1);
model.Pb = diag([100 1 100 1].^2);
model.lambdab = 0.1;
model.lambdau = 10;

% Gating parameters
P_G= 0.999;                    %gate size in percentage
model.gamma= chi2inv(P_G,model.z_dim);   %inv chi^2 dn gamma value

