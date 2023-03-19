function [y1, y2, trc] = controller2( x,wd,w,tur1,Delta_tau_w,ts )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

m11 = 50.5; m22 = 84.36; m33 = 17.21;
d11 = 151.57; d22 = 132.5; d33 = 34.56;
u = x(1); v = x(2); r = x(3);
fr = -d33*r/m33+(m11-m22)*u*v/m33;

% ESO
persistent what fwhat kw1 kw2 kw wdf Twdf lambda k_l
if isempty(what)
    kw = 5;
    what = 0; fwhat = 0; kw1 = 15; kw2 = 30;
    wdf = wd; Twdf = 0.1;
    lambda = 0; k_l = 15;
end

wdf_dot = -(wdf-wd)/Twdf; wdf = euler2(wdf_dot,wdf,ts);

what_dot = fwhat+kw1*(w-what)+tur1/m33; what = euler2(what_dot,what,ts);
fwhat_dot = kw2*(w-what); fwhat = euler2(fwhat_dot,fwhat,ts);

% 辅助动态系统
lambda_dot = (-k_l*lambda-Delta_tau_w)/m33;
lambda = euler2(lambda_dot,lambda,ts);

trc = -k_l*lambda+(-fwhat+wdf_dot-kw*(w-wdf-lambda))*m33;

% trc = (-fr-kw*(w-wd))*m33;

y1 = what; y2 = fwhat;

end

