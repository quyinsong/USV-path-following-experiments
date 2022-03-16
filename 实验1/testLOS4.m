% line of sight guidence used in straight line path following
% 相比于testLOS1，该仿真增加了纵向速度控制器的设计
% 相比于testLOS02，该仿真考虑了海流干扰因素的影响，采用ILOS制导率产生LOS角
% 参考文献： 1.《Integral LOS Control for Path Following of Underactuated Marine Surface Vessels in the Presence of Constant Ocean Currents》
% 2. MODE LING, IDENTIFICATION, AND ADAPTIVE MANEUVERING OF CYBERSffiP 11: A COMPLETE DESIGN WITH EXPERIMENTS
% control law use PID
% USV海浪干扰下的模型为 USV02
% 相较于testLOS3，采用反步法设计控制器，并利用TD来计算一阶导数，并且考虑执行器饱和
% 设计辅助动态系统来补偿输入饱和，参考文献： 滑膜变结构控制 MATLAB仿真　　作者：刘金琨
clc;
clear ;
close all;
%% USV parameters 
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
% ------------------------------------------------------
Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
Xuuu=-5.86643;       Yvv=-36.47287;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;               
% ----------------------------------------------------


%% initial
% generate point sets------------------------------------
% point_database =[0 0; 40 40; 80 40; 90 20; 90 10; 80 0; -20 -30; -20 60; 80 60; 80 90; -40 90; ...
%                  -40 -45; 0 -45; 0 0]'; % 任意不规则路径

% point_database =[0 0; 40 40; 80 40; 90 20; 90 10; 80 0; -20 -30; 40 -30]'; % 半心形路径
% point_database =[5 5; 25 5; 35 5+10*sqrt(3); 35 25+10*sqrt(3); 25 25+20*sqrt(3); 5 25+20*sqrt(3)]'; % 正五边形路径
% point_database =[5 0; 5 20; 25 20; 25 40; 5 40; 5 50; ]';  % 梳状路径
% point_database =[5 0; 60 60; 30 5]'; % 定点返航
% point_database =[5 0; 60 60]'; % 直线1
point_database =[10 0; 10 100]'; % 直线2
% kk = 0:2*pi/4:2*pi;
% xx = 20*cos(kk)+40;
% yy = 20*sin(kk)+40;
% point_database = [xx;yy]; % 圆形路径

% kk = 0:2*pi/20:2*pi;
% xx = 5*kk;
% yy = 10*sin(2*kk)+20;
% point_database = [xx;yy];  % 正弦路径

%--------------------------------------------------------- 

pointer = 1; % 航迹点指针，初始值为1，指向第一个航迹点
Pk = point_database(:,pointer);
Pk1 = point_database(:,pointer+1);
afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));
%----------------------------------------------
ts =0.01; % 采样时间
tfinal=150; % 仿真结束时间
Ns=tfinal/ts; % 仿真步数
% USV
Vc = 0.8; betac = 0*pi/180;
d = [0 0 0]'; % 外界干扰
x=[0 0 0 0 0 0]'; % USV初始状态
Tu_max = 50; Tu_min = 0; Tr_max =10;
yint = 0;  % 运动学控制器漂角自适应初始值
% TD 初始化
R = 1000; h=0.01;
x1 = [0 0]'; x2 = [0 0]';
x3 = [0 0]'; x4 = [0 0]';
% auxiliary dynamic system
lu = [0 0]'; lr = [0 0]'; 
cr1 = 100; cr2 = 100; cu1 = 100; cu2 = 100; br = 0.5; bu = 0.5;
% ctr parameters
k1 = 8; k2 = 10; k3 = 4;
% 申请内存存储时间序列
tout = zeros(Ns,6); % 时间序列
xout=zeros(Ns,6); % USV状态
YE=zeros(Ns,1); % 横向误差
Ek=zeros(Ns,1); % 艏向角误差  
Tlout=zeros(Ns,3); % 限幅后控制力和力矩
Tout=zeros(Ns,3); % 限幅前控制力和力矩
D_t = zeros(Ns,2); % 饱和补偿

%% simulation
disp('Simulation ...');
for k=1:1:Ns
    % ILOS law
    deta=3;
    Rk = sqrt((x(4)-Pk1(1))^2+(x(5)-Pk1(2))^2);
    if pointer ~= length(point_database(1,:))-1
        if Rk<= 5 % LOS航迹点切换半径为5m
            pointer = pointer+1;
            Pk = point_database(:,pointer);
            Pk1 = point_database(:,pointer+1);
        end
    end
    afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));
    % 防止LOS路径切换时产生较大的突变------------
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
    sigma = 0.3;
    psai_LOS = atan2(-ye-deta*yint,deta);  % ILOS制导率
    yint_d = sigma*ye/((ye+sigma*yint)^2+deta^2);
    yint = euler2(yint_d,yint,ts);
    psaid=afak+psai_LOS;
    % TD1
    x1_d = [x1(2);
            fhan(x1(1)-psaid,x1(2),R,h) ];
    x1 = euler2(x1_d,x1,ts);  % 跟踪期望艏向角及其一阶导
    % yaw control law
    u = x(1);v=x(2);r= x(3); psi = x(6);
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
    d11=-Xu-Xuu*abs(u)-Xuuu*u^2;
    d22=-Yv-Yvv*abs(v)-Yrv*abs(r);
    d23=-Yr-Yvr*abs(v)-Yrr*abs(r);
    d32=-Nv-Nvv*abs(v)-Nrv*abs(r);
    d33=-Nr-Nvr*abs(v)-Nrr*abs(r);
    m0 = m22*m33-m23*m32;
    %------------------------------------------------------
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
    z1 = psi-psaid-lr(1);
    afa1 = x1(2)-k1*z1+x3(2);
    z2 = r-afa1;
    % TD2
    x2_d = [x2(2);
            fhan(x2(1)-afa1,x2(2),R,h) ];
    x2 = euler2(x2_d,x2,ts);  % 跟踪afa1及其一阶导
    % yaw control law
    tr=(-fr+x2(2)-k2*z2-z1)*m0/m22; % 艏向控制力矩  
    % surge control law
    ud = 2;
    fu = (-c13*r-d11*u)/m11;
    z3 = u-ud-lu(1);
    tu = m11*(-fu-k3*z3+x4(2)); % 纵向速度控制器
    % input limited
    if tu > Tu_max
        tul = Tu_max*sign(tu); 
    elseif tu < Tu_min
        tul = 0;
    else
        tul = tu;
    end
    if abs(tr) > Tr_max
        trl = Tr_max*sign(tr); 
    else 
        trl = tr;
    end
    % auxiliary dynamic system
    D_r = trl - tr;  D_u = tul - tu;
    lr_d = [-cr1*lr(1)+lr(2); -cr2*lr(2)+br*D_r];
    lr = euler2(lr_d,lr,ts);
    lu_d = [-cu1*lu(1)+lu(2); -cu2*lu(2)+bu*D_u];
    lu = euler2(lu_d,lu,ts);
    % TD3
    x3_d = [x3(2);
            fhan(x3(1)-lr(1),x3(2),R,h) ];
    x3 = euler2(x3_d,x3,ts);  % 跟踪lr(1)及其一阶导
    % TD4
    x4_d = [x4(2);
            fhan(x4(1)-lr(1),x4(2),R,h) ];
    x4 = euler2(x4_d,x4,ts);  % 跟踪lu(1)及其一阶导   
    % USV
    tao=[tul 0 trl]';
    x=USV02(x,tao,[Vc,betac]',d);
    % store time series
    tout(k)=k*ts;
    xout(k,:)=x';% USV状态
    YE(k)=ye; % 横向误差
    Ek(k)=z1; % 艏向角误差  
    Tlout(k,:)=tao'; % 限幅后的控制力和力矩
    Tout(k,:)=[tu 0 tr]'; % 限幅前的控制力和力矩
    D_t(k,:) = [D_u D_r];
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
    if rem(k,300)==0
        modelplot(pos,psai(k));
    end   
end
plot(E,N,'r','linewidth',1)
plot(point_database(2,:),point_database(1,:),'b',E,N,'r--','linewidth',1);
hold off
figure()
plot(point_database(2,:),point_database(1,:),'b-',E,N,'r-','linewidth',2);
xlabel('E');ylabel('N');
figure()
plot(tout,psai*180/pi,'r','linewidth',2);
xlabel('time/s');ylabel('psai/deg');
figure()
plot(tout,u,'r','linewidth',2);
xlabel('time/s');ylabel('u (m/s)');
figure()
plot(tout,Tlout(:,1),'r-',tout,Tlout(:,3),'b-','linewidth',2);
legend('surge force','yaw torch');
figure()
plot(tout,YE,'r','linewidth',2);
xlabel('time/s');ylabel('YE (m)');
figure()
subplot(211)
plot(tout,Tlout(:,1),'r-',tout,Tout(:,1),'b--',tout,D_t(:,1),'g--','linewidth',2)
legend('tul','tu','D_u');
subplot(212)
plot(tout,Tlout(:,3),'r-',tout,Tout(:,3),'b--',tout,D_t(:,2),'g--','linewidth',2)
legend('trl','tr','D_r');

