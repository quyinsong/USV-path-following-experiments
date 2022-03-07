% Author: Quyinsong
% Data: 5th March 2022
% test the USV function

% 在测试USV02中，加入海流干扰因素
clc
clear all
close all
%% USV parameters
% -------------------惯性参数-----------------------------
global m; global Xudot; global Nvdot; global Iz; global Yvdot; global Nrdot; global xg; global Yrdot;
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
global m11; global m22; global m23; global m32; global m33; global m0;             
m11 = 25.8; m22 = 33.8; m23 = 1.0948; m32 = m23;  m33 = 2.76; m0 = 92.0894; % M 参数
% m11 = m-Xudot; m22 = m-Yvdot; m23 = m*xg-Yrdot;
% m32 = m*xg-Nvdot; m33 = Iz-Nrdot; m0 = m22*m33-m23*m32;  
% -------------------水动力参数---------------------------
global Xu; global Xuu; global Yv; global Yr; global Yvv; global Yrv; global Yvr; 
global Yrr; global Nv; global Nr; global Nvv; global Nrv; global Nvr; global Nrr; 
Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
                     Yvv=-36.47287;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;    
%% initial
ts=0.01;                 % sample time
tfinal =30;              % simulation final time
Ns =tfinal/ts;           % step number of simulation
Vc = 0; betac = 120*pi/180; % current speed and orientation in inertia coordinate
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
    
    if k*ts>=5
       tao = [40 0 5]';
    end
   
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


    