function [y1, y2 ,tuc] = controller1( Ud,U,tuc1,Delta_tau_u, ts )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

m11 = 50.5;

% ESO
persistent Uhat fUhat kU1 kU2 kU lambda k_l Udf TUd
if isempty(Uhat)
    kU = 1;
    Uhat = U; fUhat = 0; kU1 = 5; kU2 = 300;
    lambda = 0; k_l = 2;
    Udf = Ud; TUd = 0.1;
end

Udf_dot = -(Udf-Ud)/TUd; Udf = euler2(Udf_dot,Udf,ts);

Uhat_dot = fUhat+kU1*(U-Uhat)+tuc1/m11; Uhat = euler2(Uhat_dot,Uhat,ts);
fUhat_dot = kU2*(U-Uhat); fUhat = euler2(fUhat_dot,fUhat,ts);
% 辅助动态系统
lambda_dot = (-k_l*lambda-Delta_tau_u)/m11;
lambda = euler2(lambda_dot,lambda,ts);

tuc = -k_l*lambda+(-fUhat-kU*(U-Ud-lambda)+Udf_dot)*m11;
y1 = Uhat; y2 = fUhat; 

end

