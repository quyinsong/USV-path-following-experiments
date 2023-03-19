function [y, Hx_hat ] = FaultEstimator2( miu1,tau,B,ts, Fx )
%FAULTESTIMATOR 此处显示有关此函数的摘要
%   此处显示详细说明
persistent miu1_hat Q_hat L G
if isempty(miu1_hat)
    miu1_hat = [0 0]';
    Q_hat = [0 0]';
    L = diag([5 5]);
    G = diag([15 15]);
end

Hx_hat = -B\(Q_hat+G*miu1_hat);
miu1_hat_dot = Fx+tau-B*Hx_hat+L*(miu1-miu1_hat);
Q_hat_dot = -G*Q_hat-G*(Fx+tau+G*miu1_hat);


x_dot = [miu1_hat_dot' Q_hat_dot']';
x = [miu1_hat' Q_hat']';
x = euler2(x_dot,x,ts);

miu1_hat = x(1:2);
Q_hat = x(3:4);

y = miu1_hat;

end

