function [ wnew, rnew, xnew, Pnew ] = ppp_update( lambdau, xu, Pu, z, model )
%Create a new track for each measurement by updating PPP with measurement
% Interpret sizes from inputs
Pd = model.Pd;
H = model.H;
R = model.R;
lambda_fa = model.lambda_fa;

[stateDimensions,nu] = size(xu);
[measDimensions,m] = size(z);

% Allocate memory for new tracks
wnew = zeros(m,1);
rnew = zeros(m,1);
xnew = zeros(stateDimensions,m);
Pnew = zeros(stateDimensions,stateDimensions,m);

% Allocate temporary working for new tracks
Sk = zeros(measDimensions,measDimensions,nu);
Kk = zeros(stateDimensions,measDimensions,nu);
Pk = zeros(stateDimensions,stateDimensions,nu);
ck = zeros(nu,1);
yk = zeros(stateDimensions,nu);

% Create a new track for each measurement by updating PPP with measurement
for k = 1:nu
    Sk(:,:,k) = H*Pu(:,:,k)*H' + R;
    Kk(:,:,k) = Pu(:,:,k)*H'/Sk(:,:,k);
    Pk(:,:,k) = Pu(:,:,k) - Kk(:,:,k)*Sk(:,:,k)*Kk(:,:,k)';
end
C = zeros(m,1);
for j = 1:m
    for k = 1:nu
        v = z(:,j) - H*xu(:,k);
        try
            ck(k) = lambdau(k)*Pd*mvnpdf(v,0,Sk(:,:,k));
        catch
            ck(k) = eps;
        end
        yk(:,k) = xu(:,k) + Kk(:,:,k)*v;
    end
    
    C(j) = sum(ck);
    wnew(j) = C(j) + lambda_fa;
    rnew(j) = C(j)/wnew(j);
    ck = ck/C(j);
    xnew(:,j) = yk*ck;
    
    for k = 1:nu
        v = xnew(:,j) - yk(:,k);
        Pnew(:,:,j) = Pnew(:,:,j) + ck(k)*(Pk(:,:,k) + v*v');
    end
    
end

end

