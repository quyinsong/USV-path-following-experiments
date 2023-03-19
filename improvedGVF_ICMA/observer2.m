function [ y1,vchat_dot,y2, phat_dot ] = observer2( p,u,psi,ts )
%OBSERVER [ Vchat,p_tilde ] = observer( p,U,X,ts )
%   海流观测器 返回海流观测值
persistent phat vchat K1 K2
if isempty(phat)
    phat = p; vchat = p-phat;
    K1 = diag([8 8]);
    K2 = diag([10 10]);
end

mn = [cos(psi) sin(psi)]';
p_tilde = p-phat;
phat_dot = u*mn+vchat+K1*p_tilde;
vchat_dot = K2*p_tilde;
phat = euler2(phat_dot,phat,ts);
vchat = euler2(vchat_dot,vchat,ts);

y1 = vchat;
y2 = phat;

end

