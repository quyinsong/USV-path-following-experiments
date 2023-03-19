% test fault-tolerant control control experiments
% date: 7/6/2022
% Author: Quyinsong
clc
clear all
close all

ts = 0.02;
tfinal = 100;
Ns = tfinal/ts;
%% With AFTC
tf = 30;
x = [0 0 0 8 2 30*pi/180]'; w = 0;
Tmax = 113; Tmin = -55;
% fault detection
D = [0;0];
L = diag([40 40]);
D_up = 0.005;
Tc1 = [0 0]';
tauw = 1;
% adaptive RBFNN fault observer 
miu1_hat = [0 0]'; Hx_hat = [0 0]'; Hx = [0 0]'; Q_hat = [0 0]'; z = 0;
% Auxliary system
A = diag([50 50]);
lamda = [0 0]';

m11 = 50.05; m33 = 17.21; Fx = [0 0]';
B = [1/m11 1/m11; 0.26/m33 -0.26/m33];
% LOS
pc = 3; 

k1 = 3; delta = 2;  T1 = 0.1; T2 = 0.1; T3 = 0.2;
k2 = 0.8;  x1 = 0; x2 = 0;  
ud_min = 0.6; ud_max = 0.6;
x3 = [0;0]; x3_d = [0;0];
% control law
K = diag([1 1]);

for k = 1:Ns
    t = (k-1)*ts;
    % diturbances
    switch tauw
        case 1
            tau_w = [1.2+0.8*sin(t)*cos(0.5*t) 1.5+0.5*sin(t)*cos(0.5*t) 0.8+0.5*sin(t)*cos(0.5*t)]';
        case 0
            tau_w = [0 0 0]';
    end
    
    % state
    xn = x(4); yn = x(5); psi = x(6); u = x(1); v = x(2); r = x(3); U = sqrt(u^2+v^2); beta = atan2(v,u);
    % LOS guidance
    [xd, yd, Xsf_d, psif, kc, up , xe, ye, w] = LOS( pc, k1, delta, xn, yn, psi, U, beta, w, ts );
    % desired surge speed
    ud = ud_max;
    % backstepping control law (virtual control law)
    [ x1, x1_dot] = TD1( beta-Xsf_d,x1,T1,ts );
    Xsf = psi+beta-psif;
    Xsf_tilde = ssa(Xsf-Xsf_d); chi_tilde = Xsf_tilde;
    alpha_r = -x1_dot+kc*up-k2*Xsf_tilde;
    [ x2, x2_dot] = TD1( alpha_r,x2,T2,ts );
    e_r = r - x2; r_tilde = e_r;
    e_u = u-ud; u_tilde = e_u;
    v_tilde = v;
    
    % fault detection and diagnosis
    miu1 = [u,r]';
%     Fx = [1.6855*v*r-3.1483*u
%           -1.9936*u*v-2.0081*r];
    
    [D,D_tilde] = FaultDetection( D, Fx, miu1, B*Tc1, ts, L);
    % decision mechanism
    Hx_tilde = Hx-Hx_hat;
    
    if norm(D_tilde) >= D_up
        flag = 1;
    else
        flag = 0;
    end
    %  fault observer       
    Ed_dot = [0,x2_dot]';
    E = [e_u,e_r]'-lamda;
    switch flag
        case 1
            [miu1_hat, Hx_hat ] = FaultEstimator2( miu1,B*Tc1, B ,ts, Fx );
        case 0
            Hx_hat = [0 0]';
    end
    
    if flag == 1 && norm(Hx_tilde) <= 6
        tau = -Fx-K*E+Ed_dot-A*lamda-0.05*sign(E)+B*Hx_hat;
    else
        tau = -Fx-K*E+Ed_dot-A*lamda-0.05*sign(E);
    end

    Tc = B\tau;
    Tc0 = Tc;
    % Ö´ÐÐÆ÷±¥ºÍ
    for i = 1:2
        if Tc(i)>=Tmax
            Tc(i) = Tmax;
        elseif Tc(i)<=Tmin
            Tc(i)=Tmin;
        end
    end
    Tc1 = Tc;
    % Auxliary system
    lamda_dot = -A*lamda-B*(Tc0-Tc1);
    lamda = euler2(lamda_dot,lamda,ts);
    % ¹ßÐÔ»·½Ú
    for i=1:2
        [x3(i),x3_d(i)] = TD1(Tc1(i),x3(i),T3,ts);
        Tc1(i) = x3(i);
    end
    % USV
    [x, Tf, Hx, Fx] = USV( x, Tc1, tau_w, ts, t, tf );
    % estimate value
    xout_FTC(k,:) = [t x' Tc1' Tf' xd yd xe ye D_tilde' norm(D_tilde) flag miu1_hat' Hx_hat' Hx' norm(Hx_tilde)...
        u_tilde v_tilde chi_tilde r_tilde];
end
u_FTC = xout_FTC(:,2);
v_FTC = xout_FTC(:,3);
r_FTC = xout_FTC(:,4);
xn_FTC = xout_FTC(:,5);
yn_FTC = xout_FTC(:,6);
psi_FTC = xout_FTC(:,7);
Tc1_FTC = xout_FTC(:,8:9);
Tf_FTC = xout_FTC(:,10:11);
xd_FTC = xout_FTC(:,12);
yd_FTC = xout_FTC(:,13);
xe_FTC = xout_FTC(:,14);
ye_FTC = xout_FTC(:,15);
D_tilde_FTC = xout_FTC(:,16:17);
norm_Dtilde_FTC = xout_FTC(:,18);
flag = xout_FTC(:,19);
miu1_hat_FTC = xout_FTC(:,20:21);
Hx_hat_FTC = xout_FTC(:,22:23);
Hx_FTC = xout_FTC(:,24:25);
norm_Hxtilde_FTC = xout_FTC(:,26);
u_tilde_FTC = xout_FTC(:,27);
v_tilde_FTC = xout_FTC(:,28);
chi_tilde_FTC = xout_FTC(:,29);
r_tilde_FTC = xout_FTC(:,30);
Hx_tilde_FTC = Hx_FTC-Hx_hat_FTC;


%% Without AFTC
% initial
x = [0 0 0 8 2 30*pi/180]'; w = 0;
Tmax = 113; Tmin = -55;

% Auxliary system
A = diag([50 50]);
lamda = [0 0]';

m11 = 50.05; m33 = 17.21; Fx = [0 0]';
B = [1/m11 1/m11; 0.26/m33 -0.26/m33];
% LOS
tau = [0 0]';
k1 = 3; delta = 2;  T1 = 0.1; T2 = 0.1; T3 = 0.2;
k2 = 0.8;  x1 = 0; x2 = 0; 

x3 = [0;0]; x3_d = [0;0];
% control law
K = diag([1 1]);

for k = 1:Ns
    t = (k-1)*ts;
    % diturbances
    switch tauw
        case 1
            tau_w = [1.2+0.8*sin(t)*cos(0.5*t) 1.5+0.5*sin(t)*cos(0.5*t) 0.8+0.5*sin(t)*cos(0.5*t)]';
        case 0
            tau_w = [0 0 0]';
    end
    
    % state
    xn = x(4); yn = x(5); psi = x(6); u = x(1); v = x(2); r = x(3); U = sqrt(u^2+v^2); beta = atan2(v,u);
    % LOS guidance
    [xd, yd, Xsf_d, psif, kc, up , xe, ye, w] = LOS( pc, k1, delta, xn, yn, psi, U, beta, w, ts );
    % desired surge speed
    ud = ud_max;
    % backstepping control law (virtual control law)
    [ x1, x1_dot] = TD1( beta-Xsf_d,x1,T1,ts );
    Xsf = psi+beta-psif;
    Xsf_tilde = ssa(Xsf-Xsf_d); chi_tilde = Xsf_tilde;
    alpha_r = -x1_dot+kc*up-k2*Xsf_tilde;
    [ x2, x2_dot] = TD1( alpha_r,x2,T2,ts );
    e_r = r - x2; r_tilde = e_r;
    e_u = u-ud; u_tilde = e_u;
    v_tilde = v;
    % fault detection and diagnosis
    miu1 = [u,r]';
%     Fx = [1.6855*v*r-3.1483*u
%           -1.9936*u*v-2.0081*r];
    Ed_dot = [0,x2_dot]';
    E = [e_u,e_r]'-lamda;
    tau = -Fx-K*E+Ed_dot-A*lamda-0.15*sign(E);
    
    Tc = B\tau;
    Tc0 = Tc;
    % Ö´ÐÐÆ÷±¥ºÍ
    for i = 1:2
        if Tc(i)>=Tmax
            Tc(i) = Tmax;
        elseif Tc(i)<=Tmin
            Tc(i)=Tmin;
        end
    end
    Tc1 = Tc;
    % Auxliary system
    lamda_dot = -A*lamda-B*(Tc0-Tc1);
    lamda = euler2(lamda_dot,lamda,ts);
    % ¹ßÐÔ»·½Ú
    for i=1:2
        [x3(i),x3_d(i)] = TD1(Tc1(i),x3(i),T3,ts);
        Tc1(i) = x3(i);
    end
    % USV
    [x, Tf, Hx, Fx] = USV( x, Tc1, tau_w, ts, t, tf );
    % estimate value
    xout(k,:) = [t x' Tc1' Tf' xd yd xe ye u_tilde v_tilde chi_tilde r_tilde];
end
t = xout(:,1);
u = xout(:,2);
v = xout(:,3);
r = xout(:,4);
xn = xout(:,5);
yn = xout(:,6);
psi = xout(:,7);
Tc1 = xout(:,8:9);
Tf = xout(:,10:11);
xd = xout(:,12);
yd = xout(:,13);
xe = xout(:,14);
ye = xout(:,15);
u_tilde = xout(:,16);
v_tilde = xout(:,17);
chi_tilde = xout(:,18);
r_tilde = xout(:,19);


%% PLOTS
% Initial
fontsize = 12.5; fontname = 'TimesNewRoman'; linewid = 1.5;
figure(1);hold on
xrange=[-5 5 25]; yrange = [-5 5 25];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

plot(yd_FTC,xd_FTC,'k-','linewid',linewid);
plot(yn,xn,'b--','linewid',linewid);
plot(yn_FTC,xn_FTC,'r--','linewid',linewid);
h1=legend('desired path','real path-without AFTC','real path-with AFTC');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
y_tick = 16.2205; x_tick = 2.4914;
text(x_tick+0.5,y_tick,'\leftarrow Faults occur','color','r','FontSize',12,'linewid',2.5);
plot(x_tick,y_tick,'*','MarkerSize',15,'color','r')

for k=1:1:Ns
    pos =[xn_FTC(k) yn_FTC(k)]';
    if k==1
        modelplot(pos,psi_FTC(k),xrange,yrange);
    end
    if rem(k,500)==0
        modelplot(pos,psi_FTC(k),xrange,yrange);
    end   
end

for k=1:1:Ns
    pos =[xn(k) yn(k)]';
    if k==1
        modelplot(pos,psi(k),xrange,yrange);
    end
    if rem(k,500)==0
        modelplot1(pos,psi(k),xrange,yrange);
    end   
end



set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);

xlabel('y (m)','FontSize',fontsize,'FontName',fontname);
ylabel('x (m)','FontSize',fontsize,'FontName',fontname);

% axes('Position',[0.2,0.6,0.3,0.3]); % Éú³É×ÓÍ¼   
% plot(yd_FTC,xd_FTC,'k-','linewid',linewid); hold on
% plot(yn,xn,'b--','linewid',linewid);
% plot(yn_FTC,xn_FTC,'r--','linewid',linewid);
% xlim([20,26]); % ÉèÖÃ×ø±êÖá·¶Î§  
% ylim([20,26]); hold off

grid
box on
hold off;

figure(2)
plot(t,u,'b-',t,u_FTC,'r-','linewid',linewid);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylabel('surge speed (m/s)','FontSize',fontsize,'FontName',fontname);
h1=legend('u-without AFTC','u-with AFTC');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
figure(3)
plot(t,xe,'g-',t,ye,'k-',t,xe_FTC,'r-',t,ye_FTC,'b-','linewid',linewid);
h1=legend('xe-without AFTC','ye-without AFTC','xe-with AFTC','ye-with AFTC');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylabel('tracking errors (m)','FontSize',fontsize,'FontName',fontname);
figure(4)
Threshold = D_up*ones(length(t),1);
plot(t,D_tilde_FTC(:,1),'r-',t,D_tilde_FTC(:,2),'b-',t,norm_Dtilde_FTC,'k-',t,Threshold','m--','linewid',linewid)
h1=legend('$$\tilde{D}_1$$','$$\tilde{D}_2$$','$$\|\tilde{D}\|$$','$$Threshold$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylabel('Faults detection residual','FontSize',fontsize,'FontName',fontname);
figure(5)
subplot(2,1,1),plot(t,Tc1_FTC(:,1),'r-',t,Tc1_FTC(:,2),'b-',t,Tf_FTC(:,1),'g--',t,Tf_FTC(:,2),'c--','linewid',linewid)
h1=legend('$$T_{1}$$','$$T_{2}$$','$$T_{f1}$$','$$T_{f2}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'Orientation','Horizontal');
h1=text(80,100,'with AFTC');
set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
ylabel('Thrust (N)','FontSize',fontsize,'FontName',fontname);
% xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylim([-60,120]); % ÉèÖÃ×ø±êÖá·¶Î§ 
subplot(2,1,2),plot(t,Tc1(:,1),'r-',t,Tc1(:,2),'b-',t,Tf(:,1),'g--',t,Tf(:,2),'c--','linewid',linewid)
h1=legend('$$T_{1}$$','$$T_{2}$$','$$T_{f1}$$','$$T_{f2}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'Orientation','Horizontal');
h1=text(75,100,'without AFTC');
set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
ylabel('Thrust (N)','FontSize',fontsize,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylim([-60,120]);% ÉèÖÃ×ø±êÖá·¶Î§ 
figure(6)
subplot(2,1,1),plot(t,miu1_hat_FTC(:,1),'r-',t,u_FTC,'r--','linewid',linewid); hold on
h1=legend('$$\hat{\xi}_1$$','$$\xi_1$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylabel('(m/s)','FontSize',fontsize,'FontName',fontname);
subplot(2,1,2),plot(t,miu1_hat_FTC(:,2),'b-',t,r_FTC,'b--','linewid',linewid); hold off
h1=legend('$$\hat{\xi}_2$$','$$\xi_2$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylabel('(rad/s)','FontSize',fontsize,'FontName',fontname);

figure(7)
plot(t,-Hx_hat_FTC(:,1),'r-',t,-Hx_hat_FTC(:,2),'b-','linewid',linewid); hold on
plot(t,-Hx_FTC(:,1),'r--',t,-Hx_FTC(:,2),'b--','linewid',linewid); hold off
h1 = legend('$$-\hat{H}_1$$','$$-\hat{H}_2$$','$$-H_1$$','$$-H_2$$');
set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
ylabel('Faults estimation (N)','FontSize',fontsize,'FontName',fontname);

figure(8)
Threshold = 6*ones(length(t),1);
plot(t,norm_Hxtilde_FTC,'r-',t,Threshold,'m--','linewid',linewid);
h1 = legend('$$\|\tilde{H}\|$$','$$Threshold$$');
set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$\|\tilde{H}\| (N)$$ ');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);

figure(9)

subplot(3,2,1),plot(t,xe,'b',t,xe_FTC,'r--','linewid',linewid)
h1 = legend('$$without AFTC$$','$$with AFTC$$');
set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$x_e\ $$(m)');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);
xlim([0,100]); % ÉèÖÃ×ø±êÖá·¶Î§ 
ylim([-2 0.5]);

subplot(3,2,2),plot(t,ye,'b',t,ye_FTC,'r--','linewid',linewid)
% h1 = legend('$$without AFTC$$','$$with AFTC$$');
% set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$y_e\ $$(m)');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);
xlim([0,100]); % ÉèÖÃ×ø±êÖá·¶Î§ 
ylim([-1 2]);

subplot(3,2,3),plot(t,chi_tilde,'b',t,chi_tilde_FTC,'r--','linewid',linewid)
% h1 = legend('$$without AFTC$$','$$with AFTC$$');
% set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$\tilde{\psi}\ $$(rad)');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);
xlim([0,100]); % ÉèÖÃ×ø±êÖá·¶Î§ 
ylim([-0.5 0.5]);

subplot(3,2,4),plot(t,r_tilde,'b',t,r_tilde_FTC,'r--','linewid',linewid)
% h1 = legend('$$without AFTC$$','$$with AFTC$$');
% set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$\tilde{r}\ $$(rad/s)');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);
xlim([0,100]); % ÉèÖÃ×ø±êÖá·¶Î§
ylim([-0.5 1]);

subplot(3,2,5),plot(t,u_tilde,'b',t,u_tilde_FTC,'r--','linewid',linewid)
% h1 = legend('$$without AFTC$$','$$with AFTC$$');
% set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$\tilde{u}\ $$(m/s)');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);
xlim([0,100]); % ÉèÖÃ×ø±êÖá·¶Î§ 
ylim([-0.8 0.2]);

subplot(3,2,6),plot(t,v_tilde,'b',t,v_tilde_FTC,'r--','linewid',linewid)
% h1 = legend('$$without AFTC$$','$$with AFTC$$');
% set(h1,'Interpreter','latex');set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
xlabel('Time (s)','FontSize',fontsize,'FontName',fontname);
h1=ylabel('$$v\ $$(m/s)');
set(h1,'Interpreter','latex');set(h1,'FontSize',12); set(h1,'FontName',fontname);
xlim([0,90]); % ÉèÖÃ×ø±êÖá·¶Î§ 
ylim([-0.05 0.05]);







