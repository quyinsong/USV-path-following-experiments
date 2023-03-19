function [y1, y2, trc, lambdar] = ESOcontroller2( rd,r,tur1,Delta_tau_r,ts )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

m33 = 17.21;

% ESO
persistent rhat frhat kr1 kr2 kr rdf Trdf lambda k_l
if isempty(rhat)
    kr = 2;
    rhat = 0; frhat = 0; kr1 = 5; kr2 = 20;
    rdf = rd; Trdf = 0.1;
    lambda = 0; k_l = 5;
end

rdf_dot = -(rdf-rd)/Trdf; rdf = euler2(rdf_dot,rdf,ts);

rhat_dot = frhat+kr1*(r-rhat)+tur1/m33; rhat = euler2(rhat_dot,rhat,ts);
frhat_dot = kr2*(r-rhat); frhat = euler2(frhat_dot,frhat,ts);

% 辅助动态系统
lambda_dot = (-k_l*lambda-Delta_tau_r)/m33;
lambda = euler2(lambda_dot,lambda,ts);

trc = -k_l*lambda+(-frhat+rdf_dot-kr*(r-rd-lambda))*m33;

% trc = (-fr-kw*(w-wd))*m33;

y1 = rhat; y2 = frhat;
lambdar = lambda;
end

