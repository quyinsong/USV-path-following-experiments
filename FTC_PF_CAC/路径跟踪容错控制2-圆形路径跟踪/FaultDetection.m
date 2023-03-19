function [D,D_tilde] = FaultDetection( D, Fx, miu1, tau, ts, L)
%FAULTDETECTION �˴���ʾ�йش˺�����ժҪ
%   fault-detection
D_dot = Fx+tau-L*(D-miu1);
D = euler2(D_dot,D,ts);
D_tilde = D-miu1;
end

