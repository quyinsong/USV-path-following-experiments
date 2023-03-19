function [D,D_tilde] = FaultDetection( D, Fx, miu1, tau, ts, L)
%FAULTDETECTION 此处显示有关此函数的摘要
%   fault-detection
D_dot = Fx+tau-L*(D-miu1);
D = euler2(D_dot,D,ts);
D_tilde = D-miu1;
end

