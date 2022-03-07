% Author: Quyinsong
% Data: 5th March 2022
% test the USV function

% ²âÊÔUSV03ÖÐ
clc
clear all
close all

m11 = 215; m22 = 265; m33 = 80;
Xu = 70; Xuu = 100; Yv = 100; Yvv = 200; Nr = 50; Nrr = 100;
%% initial
ts=0.01;                 % sample time
tfinal =600;              % simulation final time
Ns =tfinal/ts;           % step number of simulation

tao=[0 0 0]'; 
d=[0 0 0]';
x=[0 0 0 70 0 90*pi/180]';
u_e_sum = 0;
%% simulation start
disp('Simulation ... ');
for k=1:1:Ns
    time(k)=(k-1)*ts;
    
%     if x(6)*180/pi>=360
%         x(6)=x(6)-2*pi;
%     end
%     if x(6)*180/pi<=-360
%         x(6)=x(6)+2*pi;
%     end
%     if k*ts == 50
%         x_turn = x(4);
%         y_turn = x(5);
%     end    
%     if k*ts>=50
%        tao(3) = 5;
%     end
    u = x(1); v = x(2); r = x(3);
    ud = 4;
    u_e = u-ud;
    u_e_sum = u_e+u_e_sum;
    Kpu = 5; Kiu = 0;
    tao(1) = -m22*v*r+Xu*u+Xuu*abs(u)*u-m11*(Kpu*u_e+Kiu*u_e_sum);
    % udot = m22*v*r/m11-Xu*u/m11-Xuu*abs(u)*u/m11+Fu/m11+d(1);
    % USV03
    x=USV03(x,tao,d);
    % store time series
    xout(k,:)=x';
    Ttao(k,:)=tao';
end
%% plot
u=xout(:,1);
v=xout(:,2);
r=xout(:,3);
N=xout(:,4);
E=xout(:,5);
psai=xout(:,6);

disp('plot ...');
for k=1:1:Ns
    pos =[N(k) E(k)]';
    if k==1
        modelplot(pos,psai(k));
    end
    if rem(k,1000)==0
        modelplot(pos,psai(k));
    end   
end
plot(E,N,'r','linewidth',1)
% rectangle('position',[y_turn,x_turn,2,2],'curvature',[1 1]);
hold off;
figure(2);
plot(time,psai*180/pi,'r','linewidth',2);
xlabel('time/s');ylabel('psai/deg');
figure(3);
plot(E,N,'r','linewidth',2)
xlabel('E');ylabel('N');
figure(4);
plot(time,u,'r','linewidth',2)
xlabel('time/s');ylabel('u (m/s)');
figure(5);
plot(time,Ttao(:,1),'r',time,Ttao(:,2),'k',time,Ttao(:,3),'b','linewidth',2)
xlabel('time/s');ylabel('T (m/s)');


    