function [ y1,Vchat_dot,y2 ] = observer( p,U,X,ts )
%OBSERVER [ Vchat,p_tilde ] = observer( p,U,X,ts )
%   恒定海流观测器 返回海流观测值
persistent phat Vchat K1 K2
if isempty(phat)
    phat = p; Vchat = p-phat;
    K1 = diag([2 2]);
    K2 = diag([1 1]);
end

mn = [cos(X) sin(X)]';
p_tilde = p-phat;
phat_dot = U*mn+Vchat+K1*p_tilde;
Vchat_dot = K2*p_tilde;
phat = euler2(phat_dot,phat,ts);
Vchat = euler2(Vchat_dot,Vchat,ts);

y1 = Vchat;
y2 = phat;

end

