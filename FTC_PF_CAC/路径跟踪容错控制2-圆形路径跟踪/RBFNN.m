function Hx_hat = RBFNN( miu1,tau, miu1_hat, ts, B )
%RBFNN estimate and cancel the unknown fault and model uncertainties
% Input: miu1 = [u,r]'
%          tau = [tu,tr]'
% Weight: W
% estimate error: Ef
% control eror: E
% sampleTime: ts
persistent W nh no
if isempty(W)
    nh = 7;
    no = 2;
    W = zeros(nh,no);  
end
miu1_tilde = miu1-miu1_hat;
Z = [miu1', tau', miu1_hat']';
Lamda = 20000;
sigma = 0;
C = [-1.5 1 0.5 0 0.5 1 1.5
     -0.6 -0.4 -0.2 0 0.2 0.4 0.6
     -6 -4  -2 0 2 4 6
     -2 -1 -1.5 0 1 1.5 2 
     -1.5 1 0.5 0 0.5 1 1.5
     -0.6 -0.4 -0.2 0 0.2 0.4 0.6];

b = 1.5;
PhiZ = zeros(nh,1);
for i=1:nh
    PhiZ(i)=exp(-norm(Z-C(:,i))^2/(2*b^2));
end
Hx_hat = W'*PhiZ;
Wdot = -Lamda*(PhiZ*(B*miu1_tilde)'+sigma*W);
% Wdot = -Lamda*(PhiZ*miu1_tilde'*B+sigma*W);
W = W+Wdot*ts;

end

