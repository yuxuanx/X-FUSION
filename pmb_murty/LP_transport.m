function [ Cmin,q ] = LP_transport( C,pn,ph )

[H,N] = size(C);
c = reshape(C',H*N,1);
beq = [ph;pn];
Aeq1 = kron(eye(H),ones(1,N));
Aeq2 = repmat(eye(N),1,H);
Aeq = [Aeq1;Aeq2];
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
[x,Cmin] = linprog(abs(c),[],[],Aeq,beq,zeros(1,H*N),ones(1,H*N),options);
q = reshape(x,N,H)';

end

