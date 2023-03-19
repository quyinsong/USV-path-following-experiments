function [y1, y2 ,tuc, lambdau] = ESOcontroller1( ud,u,tuc1,Delta_tau_u, ts )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

m11 = 50.5;

% ESO
persistent uhat fuhat ku1 ku2 ku lambda k_l udf Tud
if isempty(uhat)
    ku = 5;
    uhat = u; fuhat = 0; ku1 = 5; ku2 = 20;
    lambda = 0; k_l = 20;
    udf = ud; Tud = 0.1;
end

udf_dot = -(udf-ud)/Tud; udf = euler2(udf_dot,udf,ts);

uhat_dot = fuhat+ku1*(u-uhat)+tuc1/m11; uhat = euler2(uhat_dot,uhat,ts);
fuhat_dot = ku2*(u-uhat); fuhat = euler2(fuhat_dot,fuhat,ts);
% 辅助动态系统
lambda_dot = (-k_l*lambda-Delta_tau_u)/m11;
lambda = euler2(lambda_dot,lambda,ts);

tuc = -k_l*lambda+(-fuhat-ku*(u-ud-lambda)+udf_dot)*m11;
y1 = uhat; y2 = fuhat; 
lambdau = lambda;

end

