function [ wupd, rupd, xupd, Pupd ] = mbm_update( z, r, x, P, model )

% Interpret sizes from inputs
Pd = model.Pd;
H = model.H;
R = model.R;
stateDimensions = 4;
m = size(z,2);
n = length(r);

%Update existing tracks
wupd = zeros(n,m+1);
rupd = zeros(n,m+1);
xupd = zeros(stateDimensions,n,m+1);
Pupd = zeros(stateDimensions,stateDimensions,n,m+1);


% Update existing tracks
% Create missed detection hypothesis
wupd(:,1) = 1 - r + r*(1-Pd);
rupd(:,1) = r*(1-Pd)./wupd(:,1);
xupd(:,:,1) = x;
Pupd(:,:,:,1) = P;
for i = 1:n
    % Create hypotheses with measurement updates
    S = H*P(:,:,i)*H' + R;
    K = P(:,:,i)*H'/S;
    Pplus = P(:,:,i) - K*S*K';
    
    rupd(:,2:end) = 1;
    for j = 1:m
        v = z(:,j) - H*x(:,i);
        try
            wupd(i,j+1) = r(i)*Pd*mvnpdf(v,0,S);
        catch
            wupd(i,j+1) = eps;
        end
        xupd(:,i,j+1) = x(:,i) + K*v;
        Pupd(:,:,i,j+1) = Pplus;
    end
end

end

