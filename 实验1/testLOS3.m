% line of sight guidence used in straight line path following
% �����testLOS1���÷��������������ٶȿ����������
% �����testLOS02���÷��濼���˺����������ص�Ӱ�죬����ILOS�Ƶ��ʲ���LOS��
% �ο����ף� 1.��Integral LOS Control for Path Following of Underactuated Marine Surface Vessels in the Presence of Constant Ocean Currents��
% 2. MODE LING, IDENTIFICATION, AND ADAPTIVE MANEUVERING OF CYBERSffiP 11: A COMPLETE DESIGN WITH EXPERIMENTS
% control law use PID
% USV���˸����µ�ģ��Ϊ USV02
clc;
clear ;
close all;
%% USV parameters 
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
% ------------------------------------------------------
Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
                     Yvv=-36.47287;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;               
% ----------------------------------------------------


%% initial
% generate point sets------------------------------------
% point_database =[0 0; 40 40; 80 40; 90 20; 90 10; 80 0; -20 -30; -20 60; 80 60; 80 90; -40 90; ...
%                  -40 -45; 0 -45; 0 0]'; % ���ⲻ����·��

% point_database =[0 0; 40 40; 80 40; 90 20; 90 10; 80 0; -20 -30; 40 -30]';
% point_database =[5 0; 5 20; 25 20; 25 40; 5 40; 5 50; ]';  % ��״·��

point_database =[5 0; 60 60; 5 15]'; % ���㷵��

% kk = 0:2*pi/4:2*pi;
% xx = 20*cos(kk)+40;
% yy = 20*sin(kk)+40;
% point_database = [xx;yy]; % Բ��·��

% kk = 0:2*pi/20:2*pi;
% xx = 5*kk;
% yy = 10*sin(2*kk)+20;
% point_database = [xx;yy];  % ����·��

%--------------------------------------------------------- 

pointer = 1; % ������ָ�룬��ʼֵΪ1��ָ���һ��������
Pk = point_database(:,pointer);
Pk1 = point_database(:,pointer+1);
afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));
%----------------------------------------------
ts =0.01; % ����ʱ��
tfinal=60; % �������ʱ��
Ns=tfinal/ts; % ���沽��
x=[0 0 0 0 0 pi/2]'; % USV��ʼ״̬
ek_1=2; % ����Ǹ���ǰһʱ������ʼ��
psaid_1 = 0.1; psaid_2 = 0.05; % ���������ǰ��ʱ�̳�ʼ��
yint = 0;
% �����ڴ�洢ʱ������
xout=zeros(Ns,6); % USV״̬
YE=zeros(Ns,1); % �������
Ek=zeros(Ns,1); % ��������  
Ttao=zeros(Ns,3); % ������������

%% simulation
disp('Simulation ...');
for k=1:1:Ns
    time(k)=k*ts;
    
    % ILOS law
    deta=3;
    Rk = sqrt((x(4)-Pk1(1))^2+(x(5)-Pk1(2))^2);
    if pointer ~= length(point_database(1,:))-1
        if Rk<= 3 % LOS�������л��뾶Ϊ5m
            pointer = pointer+1;
            Pk = point_database(:,pointer);
            Pk1 = point_database(:,pointer+1);
        end
    end
    afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));
    % ��ֹLOS·���л�ʱ�����ϴ��ͻ��------------
%     if pointer == 1
%        Pk0 = point_database(:,pointer);
%        Pk01 = point_database(:,pointer+1);
%     else
%        Pk0 = point_database(:,pointer-1);
%        Pk01 = point_database(:,pointer);
%     end
%     afak0 = atan2(Pk01(2)-Pk0(2),Pk01(1)-Pk0(1));
%     if abs(afak-afak0)*180/pi >= 180
%         if afak > 0
%             afak = afak - 2*pi;
%         elseif afak < 0
%             afak = afak + 2*pi;
%         end
%     end
    %--------------------------------------------------

    ye=-(x(4)-Pk(1))*sin(afak)+(x(5)-Pk(2))*cos(afak);
    beta=atan2(x(2),x(1));
    sigma = 2;
    psai_LOS = atan2(-(ye+sigma*yint)/deta,1);  % ILOS�Ƶ���
    yint_d = deta*ye/((ye+sigma*yint)^2+deta^2);
    yint = euler2(yint_d,yint,ts);
    psaid=afak+psai_LOS-beta;
    % yaw control law
    u = x(1);v=x(2);r= x(3);
    ek=x(6)-psaid;
    Kdr = 2; 
    Kpr = 10;
    % ----------------------------------------------------
    m11 = m-Xudot; 
    m22 = m-Yvdot;
    m23 = m*xg-Yrdot;
    m32 = m*xg-Nvdot;
    m33 = Iz-Nrdot;
    % -----------------------------------------------------
    c13 = -m23*r-m22*v;  % �����޸ģ�ԭ������ʽΪ��c13=-m*(xg*r+v), c23 = m*u;
    c23 = m11*u;
    c31 = -c13; c32 = -c23;
    % -----------------------------------------------------
    d11=-Xu-Xuu*abs(u);
    d22=-Yv-Yvv*abs(v)-Yrv*abs(r);
    d23=-Yr-Yvr*abs(v)-Yrr*abs(r);
    d32=-Nv-Nvv*abs(v)-Nrv*abs(r);
    d33=-Nr-Nvr*abs(v)-Nrr*abs(r);
    m0 = m22*m33-m23*m32;
    %------------------------------------------------------
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
    
    psaidd = (psaid-2*psaid_1+psaid_2)/ts^2;
    psaid_2=psaid_1; psaid_1 = psaid;
    tpid=(-Kpr*ek-Kdr*(ek-ek_1)/ts-fr+psaidd)*m0/m22; % �����������  
    tr = tpid;
    ek_1=ek;
    % surge control law
    fu = (-c13*r-d11*u)/m11;
    ud = 3;
    eu = u-ud;
    Kpu = 5;
    tu = m11*(-fu-Kpu*eu); % �����ٶȿ�����
    tao=[tu 0 tr]';
    % USV
    Vc = 0.2; betac = 30*pi/180;
    d = [0 0 0]'; % ������
    x=USV02(x,tao,[Vc,betac]',d);
    % store time series
    xout(k,:)=x';% USV״̬
    YE(k)=ye; % �������
    Ek(k)=ek; % ��������  
    Ttao(k,:)=tao'; % ������������
end

%% plot

u=xout(:,1);
v=xout(:,2);
r=xout(:,3);
N=xout(:,4);
E=xout(:,5);
psai=xout(:,6);

disp('Plot ...');
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
plot(point_database(2,:),point_database(1,:),'b',E,N,'r--','linewidth',1)
hold off
figure
plot(point_database(2,:),point_database(1,:),'b-',E,N,'r-','linewidth',2)
xlabel('E');ylabel('N');
figure
plot(time,psai*180/pi,'r','linewidth',2);
xlabel('time/s');ylabel('psai/deg');
figure
plot(time,u,'r','linewidth',2)
xlabel('time/s');ylabel('u (m/s)');
figure
plot(time,Ttao(:,1),'r',time,Ttao(:,3),'b','linewidth',2)
legend('surge force','yaw torch');
figure
plot(time,YE,'r','linewidth',2)
xlabel('time/s');ylabel('YE (m)');

