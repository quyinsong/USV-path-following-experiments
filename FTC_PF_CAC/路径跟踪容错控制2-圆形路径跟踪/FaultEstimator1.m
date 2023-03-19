function [y, Hx_hat, D_hat ] = FaultEstimator1( miu1,tau,B,ts, Fx )
%FAULTESTIMATOR 此处显示有关此函数的摘要
%   此处显示详细说明
persistent z Q Lambda1 Lambda2
if isempty(z)
    z = [0 0]';
    Q = [0 0]';
    Lambda1 = diag([10 10]);
    Lambda2 = diag([0.8 0.8]);

end
Gamma = 0;

Hx_hat = B\Lambda1*(Q-miu1);
D_hat = -Lambda2*(z-miu1);
Q_dot = -Gamma*(Q-miu1)+Fx+tau-B*Hx_hat+D_hat;
z_dot = Fx+tau-B*Hx_hat+D_hat;

x_dot = [Q_dot' z_dot']';
x = [Q' z']';
x = euler2(x_dot,x,ts);

Q = x(1:2);
z = x(3:4);
y = Q;

end

