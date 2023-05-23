% calculate S+ and S-
% 2022 5 14
% ���CKF���ݻ������ S+ S-
clc
clear all
close all
n =5;
E = eye(n);
% ��ʽ1
% k=1;
% for i=1:n-1
%     A = E(:,i+1:n);
%     S(:,k:k-1+n-i) = sqrt(0.5)*(repmat(E(:,i),[1,n-i])+A);
%     S_(:,k:k-1+n-i) = sqrt(0.5)*repmat(E(:,i),[1,n-i])-A;
%     k=k+n-i;
% end
% ��ʽ2
i=1;
for k = 2:n
    A = E(:,1:k-1);
    S(:,i:k-1+i-1) = sqrt(0.5)*(repmat(E(:,k),[1,k-1])+A);
    S_(:,i:k-1+i-1) = sqrt(0.5)*(repmat(E(:,k),[1,k-1])-A);
    i = k-1+i;
end
disp(S);
disp(S_);