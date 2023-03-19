function [Fhat, y_lambda, Tc ] = NNcontroller( wd,w,T,Tc,F,ts )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

m11 = 50.5; dp = 0.26; m33 = 17.21;
B = [1/m11 1/m11; dp/m33 -dp/m33];
DeltaT = Tc-T;
u = w(1); r = w(2);
X = [u r u*r]';  % ����������

% ��ʼ��
persistent What Kw kw A Gamma c b m n p lambda wdf Twf
if isempty(What)
    Kw = diag([4 2]);
    kw = 0.005; Gamma = 2.5;
    m = 9; n = length(X); p = 2; % mΪ��������nΪ������,pΪ�����
    What = zeros(m,p);
    c = zeros(n,m);
    b = 2*ones(1,m);
    lambda = B*DeltaT; A = diag([1 0.5]);
    wdf = wd; Twf = 0.1;
end

% һ��΢�ָ�����
wdf_dot = -(wdf-wd)/Twf;
wdf = euler2(wdf_dot,wdf,ts);

% ������̬ϵͳ
lambda_dot = -A*lambda-B*DeltaT;
lambda = euler2(lambda_dot,lambda,ts);

% ������ƽ�F
PhiX = zeros(m,1);
for i = 1:m
    PhiX(i) = exp(-norm(X-c(:,i))^2/b(i)^2);
end
Fhat = What'*PhiX;

% ������
ew = w-wd-lambda;
Tc = B\(wdf_dot-A*lambda-Fhat-Kw*ew);

% ��������
What_dot = Gamma*(PhiX*ew'-kw*What);
What = euler2(What_dot,What,ts);

% �������
y_lambda = lambda;
end

