% line of sight guidence used in straight line path following
% 相比于testLOS1，该仿真增加了纵向速度控制器的设计
% control law use PID
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
%                  -40 -45; 0 -45; 0 0]';

% point_database =[0 0; 40 40; 80 40; 90 20; 90 10; 80 0; -20 -30; 40 -30]';

kk = 0:2*pi/15:2*pi;
xx = 20*cos(kk);
yy = 20*sin(kk);
point_database = [xx;yy];

% kk = 0:2*pi/100:2*pi;
% xx = 20*kk;
% yy = 10*sin(2*kk);
% point_database = [xx;yy];

%--------------------------------------------------------- 

pointer = 1; % 航迹点指针，初始值为1，指向第一个航迹点
Pk = point_database(:,pointer);
Pk1 = point_database(:,pointer+1);
afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));

ts =0.01; % 采样时间
tfinal=50; % 仿真结束时间
Ns=tfinal/ts; % 仿真步数
x=[0 0 0 18 0 pi/2]'; % USV初始状态
ek_1=2; % 艏向角跟踪前一时刻误差初始化

psaid_1 = 0.1; psaid_2 = 0.05; % 期望艏向角前两时刻初始化

% 申请内存存储时间序列
xout=zeros(Ns,6); % USV状态
YE=zeros(Ns,1); % 横向误差
Ek=zeros(Ns,1); % 艏向角误差  
Ttao=zeros(Ns,3); % 控制力和力矩

%% simulation
disp('Simulation ...');
for k=1:1:Ns
    time(k)=k*ts;
    
    % LOS law
    deta=3;
    Rk = sqrt((x(4)-Pk1(1))^2+(x(5)-Pk1(2))^2);
    if pointer ~= length(point_database(1,:))-1
        if Rk<=5  % LOS航迹点切换半径为5m
            pointer = pointer+1;
            Pk = point_database(:,pointer);
            Pk1 = point_database(:,pointer+1);
        end
    end
    afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));
    % 防止LO路径切换时产生较大的突变------------
    if pointer == 1
       Pk0 = point_database(:,pointer);
       Pk01 = point_database(:,pointer+1);
    else
       Pk0 = point_database(:,pointer-1);
       Pk01 = point_database(:,pointer);
    end
    afak0 = atan2(Pk01(2)-Pk0(2),Pk01(1)-Pk0(1));
    if abs(afak-afak0)*180/pi >= 180
        if afak > 0
            afak = afak - 2*pi;
        elseif afak < 0
            afak = afak + 2*pi;
        end
    end
    %--------------------------------------------------

    ye=-(x(4)-Pk(1))*sin(afak)+(x(5)-Pk(2))*cos(afak);
    beta=atan2(x(2),x(1));
    psaid=afak+atan2(-ye/deta,1)-beta;
    % yaw control law
    u = x(1);v=x(2);r= x(3);
    ek=x(6)-psaid;
    Kdr = 6; 
    Kpr =4;
    % ----------------------------------------------------
    m11 = m-Xudot; 
    m22 = m-Yvdot;
    m23 = m*xg-Yrdot;
    m32 = m*xg-Nvdot;
    m33 = Iz-Nrdot;
    % -----------------------------------------------------
    c13 = -m23*r-m22*v;  % 做了修改，原来的形式为：c13=-m*(xg*r+v), c23 = m*u;
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
    tpid=(-Kpr*ek-Kdr*(ek-ek_1)/ts-fr*m0/m22+psaidd); % 艏向控制力矩  
    tr = tpid;
    ek_1=ek;
    % surge control law
    fu = (-c13*r-d11*u)/m11;
    ud = 3;
    eu = u-ud;
    Kpu = 4;
    tu = m11*(-fu-Kpu*eu); % 纵向速度控制器
    tao=[tu 0 tr]';
    % USV
    d = [0 0 0]'; % 外界干扰
    xdot=USV01(x,tao,[0,0]',d);
    % state update
    x=euler2(xdot,x,ts);
    % store time series
    xout(k,:)=x';% USV状态
    YE(k)=ye; % 横向误差
    Ek(k)=ek; % 艏向角误差  
    Ttao(k,:)=tao'; % 控制力和力矩
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

