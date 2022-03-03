% Author: Quyinsong
% Data: 3rd March 2022
% test the USV function

% 在测试USV02中，加入海流干扰因素
clc
clear all
close all
%% initial
ts=0.01;                 % sample time
tfinal =30;              % simulation final time
Ns =tfinal/ts;           % step number of simulation
Vc = 1; betac = 120*pi/180; % current speed and orientation in inertia coordinate
current=[Vc betac]';       % current

tao=[40 0 0]'; 
d=[0 0 0]';
x=[0 0 0 40 0 90*pi/180]';
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
    
%     if k*ts>=15
%        tao = [20 0 1]';
%     end
   
    % USV02
    x=USV02(x,tao,current,d);
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
    if rem(k,100)==0
        modelplot(pos,psai(k));
    end   
end
plot(E,N,'r','linewidth',1)
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


    