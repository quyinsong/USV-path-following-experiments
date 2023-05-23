% EKF UKF 三阶CKF 五阶CKF　蒙特卡洛对比仿真实验
% 分别调用LOSdemoEKF.slx  LOSdemoUKF.slx  LOSdemoCKF3.slx  LOSdemoCKF5.slx 文件
% 进行Mc次蒙特卡洛仿真
clc
close all
% EKF and UKF and CKF3  
Mc = 1;
%% EKF
disp('EKF')
sumx1 = 0;
sumy1 = 0;
sumU1 = 0;
sumX1 = 0;
sumWx1 = 0;
tic;  % 记录开始时间
for kc = 1:Mc
    disp(kc);
    sim LOSdemoEKF;
    sumx1 = sumx1+(x-xhat).^2;
    sumy1 = sumy1+(y-yhat).^2;
    sumU1 = sumU1+(U-Uhat).^2;
    sumX1 = sumX1+(X-Xhat).^2;
    sumWx1 = sumWx1+(Wx-Wxhat).^2; 
end
x1 = x; y1 = y; U1 = U; X1 = X; Wx1 = Wx;
TEKF = toc;   % 记录结束时间
tekf = TEKF/Mc;  % 计算Mc次蒙特卡洛仿真平均运行时间
N = Mc;
% 计算各个状态RMSE
    RMSE_x1 = sqrt(sumx1/N);
    RMSE_y1 = sqrt(sumy1/N);
    RMSE_U1 = sqrt(sumU1/N);
    RMSE_X1 = sqrt(sumX1/N);
    RMSE_Wx1 = sqrt(sumWx1/N);
    mean_RMSE_x1 = mean(RMSE_x1);
    mean_RMSE_y1 = mean(RMSE_y1);
    mean_RMSE_U1 = mean(RMSE_U1);
    mean_RMSE_X1 = mean(RMSE_X1);
    mean_RMSE_Wx1 = mean(RMSE_Wx1);
    var_RMSE_x1 = mean((RMSE_x1-mean_RMSE_x1).^2);
    var_RMSE_y1 = mean((RMSE_y1-mean_RMSE_y1).^2);
    var_RMSE_U1 = mean((RMSE_U1-mean_RMSE_U1).^2);
    var_RMSE_X1 = mean((RMSE_X1-mean_RMSE_X1).^2);
    var_RMSE_Wx1 = mean((RMSE_Wx1-mean_RMSE_Wx1).^2);
%% UKF
disp('UKF')
sumx2 = 0;
sumy2 = 0;
sumU2 = 0;
sumX2 = 0;
sumWx2 = 0;
tic;
for kc = 1:Mc
    disp(kc);
    sim LOSdemoUKF;
    sumx2 = sumx2+(x-xhat).^2;
    sumy2 = sumy2+(y-yhat).^2;
    sumU2 = sumU2+(U-Uhat).^2;
    sumX2 = sumX2+(X-Xhat).^2;
    sumWx2 = sumWx2+(Wx-Wxhat).^2; 
end
x2 = x; y2 = y;U2 = U; X2 = X; Wx2 = Wx;
TUKF = toc;
tukf = TUKF/Mc;
N = Mc;
    RMSE_x2 = sqrt(sumx2/N);
    RMSE_y2 = sqrt(sumy2/N);
    RMSE_U2 = sqrt(sumU2/N);
    RMSE_X2 = sqrt(sumX2/N);
    RMSE_Wx2 = sqrt(sumWx2/N);
    mean_RMSE_x2 = mean(RMSE_x2);
    mean_RMSE_y2 = mean(RMSE_y2);
    mean_RMSE_U2 = mean(RMSE_U2);
    mean_RMSE_X2 = mean(RMSE_X2);
    mean_RMSE_Wx2 = mean(RMSE_Wx2);
    var_RMSE_x2 = mean((RMSE_x2-mean_RMSE_x2).^2);
    var_RMSE_y2 = mean((RMSE_y2-mean_RMSE_y2).^2);
    var_RMSE_U2 = mean((RMSE_U2-mean_RMSE_U2).^2);
    var_RMSE_X2 = mean((RMSE_X2-mean_RMSE_X2).^2);
    var_RMSE_Wx2 = mean((RMSE_Wx2-mean_RMSE_Wx2).^2);
%% 三阶CKF 
disp('CKF3')
sumx3 = 0;
sumy3 = 0;
sumU3 = 0;
sumX3 = 0;
sumWx3 = 0;
tic;
for kc = 1:Mc
    disp(kc);
    sim LOSdemoCKF3;
    sumx3 = sumx3+(x-xhat).^2;
    sumy3 = sumy3+(y-yhat).^2;
    sumU3 = sumU3+(U-Uhat).^2;
    sumX3 = sumX3+(X-Xhat).^2;
    sumWx3 = sumWx3+(Wx-Wxhat).^2; 
end
x3 = x; y3 = y;U3 = U; X3 = X; Wx3 = Wx;
TCKF3 = toc;
tckf3 = TCKF3/Mc;
N = Mc;
    RMSE_x3 = sqrt(sumx3/N);
    RMSE_y3 = sqrt(sumy3/N);
    RMSE_U3 = sqrt(sumU3/N);
    RMSE_X3 = sqrt(sumX3/N);
    RMSE_Wx3 = sqrt(sumWx3/N);
    mean_RMSE_x3 = mean(RMSE_x3);
    mean_RMSE_y3 = mean(RMSE_y3);
    mean_RMSE_U3 = mean(RMSE_U3);
    mean_RMSE_X3 = mean(RMSE_X3);
    mean_RMSE_Wx3 = mean(RMSE_Wx3);
    var_RMSE_x3 = mean((RMSE_x3-mean_RMSE_x3).^2);
    var_RMSE_y3 = mean((RMSE_y3-mean_RMSE_y3).^2);
    var_RMSE_U3 = mean((RMSE_U3-mean_RMSE_U3).^2);
    var_RMSE_X3 = mean((RMSE_X3-mean_RMSE_X3).^2);
    var_RMSE_Wx3 = mean((RMSE_Wx3-mean_RMSE_Wx3).^2); 
%% 五阶CKF 
disp('CKF5')
sumx4 = 0;
sumy4 = 0;
sumU4 = 0;
sumX4 = 0;
sumWx4 = 0;
tic;
for kc = 1:Mc
    disp(kc);
    sim LOSdemoCKF5;
    sumx4 = sumx4+(x-xhat).^2;
    sumy4 = sumy4+(y-yhat).^2;
    sumU4 = sumU4+(U-Uhat).^2;
    sumX4 = sumX4+(X-Xhat).^2;
    sumWx4 = sumWx4+(Wx-Wxhat).^2; 
end
x4 = x; y4 = y;U4 = U; X4 = X; Wx4 = Wx;
TCKF5 = toc;
tckf5 = TCKF5/Mc;
N = Mc;
    RMSE_x4 = sqrt(sumx4/N);
    RMSE_y4 = sqrt(sumy4/N);
    RMSE_U4 = sqrt(sumU4/N);
    RMSE_X4 = sqrt(sumX4/N);
    RMSE_Wx4 = sqrt(sumWx4/N);
    mean_RMSE_x4 = mean(RMSE_x4);
    mean_RMSE_y4 = mean(RMSE_y4);
    mean_RMSE_U4 = mean(RMSE_U4);
    mean_RMSE_X4 = mean(RMSE_X4);
    mean_RMSE_Wx4 = mean(RMSE_Wx4);
    var_RMSE_x4 = mean((RMSE_x4-mean_RMSE_x4).^2);
    var_RMSE_y4 = mean((RMSE_y4-mean_RMSE_y4).^2);
    var_RMSE_U4 = mean((RMSE_U4-mean_RMSE_U4).^2);
    var_RMSE_X4 = mean((RMSE_X4-mean_RMSE_X4).^2);
    var_RMSE_Wx4 = mean((RMSE_Wx4-mean_RMSE_Wx4).^2); 
%% LOSdemo0
sim LOSdemo0;  % 为了和上面滤波方法进行对照
x5 = x; y5 = y;U5 = U; X5 = X; Wx5 = Wx;

%% disp
disp('EKF mean RMSE and var RMSE : ');
disp(['mean_RMSE_x1 : ', num2str(mean_RMSE_x1)]);
disp(['mean_RMSE_y1 : ', num2str(mean_RMSE_y1)]);
disp(['mean_RMSE_U1 : ', num2str(mean_RMSE_U1)]);
disp(['mean_RMSE_X1 : ', num2str(mean_RMSE_X1)]);
disp(['mean_RMSE_Wx1 : ', num2str(mean_RMSE_Wx1)]);
disp(['var_RMSE_x1 : ', num2str(var_RMSE_x1)]);
disp(['var_RMSE_y1 : ', num2str(var_RMSE_y1)]);
disp(['var_RMSE_U1 : ', num2str(var_RMSE_U1)]);
disp(['var_RMSE_X1 : ', num2str(var_RMSE_X1)]);
disp(['var_RMSE_Wx1 : ', num2str(var_RMSE_Wx1)]);
disp('UKF mean RMSE and var RMSE : ');
disp(['mean_RMSE_x2 : ', num2str(mean_RMSE_x2)]);
disp(['mean_RMSE_y2 : ', num2str(mean_RMSE_y2)]);
disp(['mean_RMSE_U2 : ', num2str(mean_RMSE_U2)]);
disp(['mean_RMSE_X2 : ', num2str(mean_RMSE_X2)]);
disp(['mean_RMSE_Wx2 : ', num2str(mean_RMSE_Wx2)]);
disp(['var_RMSE_x2 : ', num2str(var_RMSE_x2)]);
disp(['var_RMSE_y2 : ', num2str(var_RMSE_y2)]);
disp(['var_RMSE_U2 : ', num2str(var_RMSE_U2)]);
disp(['var_RMSE_X2 : ', num2str(var_RMSE_X2)]);
disp(['var_RMSE_Wx2 : ', num2str(var_RMSE_Wx2)]);
disp('CKF3 mean RMSE and var RMSE : ');
disp(['mean_RMSE_x3 : ', num2str(mean_RMSE_x3)]);
disp(['mean_RMSE_y3 : ', num2str(mean_RMSE_y3)]);
disp(['mean_RMSE_U3 : ', num2str(mean_RMSE_U3)]);
disp(['mean_RMSE_X3 : ', num2str(mean_RMSE_X3)]);
disp(['mean_RMSE_Wx3 : ', num2str(mean_RMSE_Wx3)]);
disp(['var_RMSE_x3 : ', num2str(var_RMSE_x3)]);
disp(['var_RMSE_y3 : ', num2str(var_RMSE_y3)]);
disp(['var_RMSE_U3 : ', num2str(var_RMSE_U3)]);
disp(['var_RMSE_X3 : ', num2str(var_RMSE_X3)]);
disp(['var_RMSE_Wx3 : ', num2str(var_RMSE_Wx3)]);
disp('CKF5 mean RMSE and var RMSE : ');
disp(['mean_RMSE_x4 : ', num2str(mean_RMSE_x4)]);
disp(['mean_RMSE_y4 : ', num2str(mean_RMSE_y4)]);
disp(['mean_RMSE_U4 : ', num2str(mean_RMSE_U4)]);
disp(['mean_RMSE_X4 : ', num2str(mean_RMSE_X4)]);
disp(['mean_RMSE_Wx4 : ', num2str(mean_RMSE_Wx4)]);
disp(['var_RMSE_x4 : ', num2str(var_RMSE_x4)]);
disp(['var_RMSE_y4 : ', num2str(var_RMSE_y4)]);
disp(['var_RMSE_U4 : ', num2str(var_RMSE_U4)]);
disp(['var_RMSE_X4 : ', num2str(var_RMSE_X4)]);
disp(['var_RMSE_Wx4 : ', num2str(var_RMSE_Wx4)]);
disp('mean run time : ');
disp(['EKF mean time: ',num2str(tekf)]);
disp(['UKF mean time: ',num2str(tukf)]);
disp(['CKF3 mean time: ',num2str(tckf3)]);
disp(['CKF5 mean time: ',num2str(tckf5)]);
%% plot
disp('Plot ...');
for k=1:1:length(x)
    pos =[x(k) y(k)]';
    if k==1
        modelplot(pos,psi(k));
    end
    if rem(k,1000)==0
        modelplot(pos,psi(k));
    end   
end
plot(y,x,'r-','linewidth',1);
point_database =[300 300; 1500 300; 2900 1900; 2900 3900; 1500 5500; 300 5500]';
% plot(point_database(2,:),point_database(1,:),'b-','linewidth',1)

for k = 1:length(point_database(1,:))
    plot(point_database(2,k),point_database(1,k),'ro','linewidth',4);
end
title('Map-XY');
xlabel('E(m)');ylabel('N(m)');
hold off

figure()
plot(t,RMSE_x1,'r-',t,RMSE_x2,'b-',t,RMSE_x3,'g-',t,RMSE_x4,'k-','linewidth',1.5);
legend('RMSE(EKF x)','RMSE(UKF x)','RMSE(CKF3 x)','RMSE(CKF5 x)');
title('North position RMSE'); 
xlabel('t(s)');ylabel('x-RMSE(m)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([0,6]);
% axes('Position',[0.2,0.6,0.3,0.3]); % 生成子图   
% plot(t,RMSE_x1,'r-',t,RMSE_x2,'b--',t,RMSE_x3,'g--',t,RMSE_x4,'k--','linewidth',1.5);% 绘制局部曲线图                        
% xlim([250,350]); % 设置坐标轴范围  
% ylim([0.8,0.82]);

figure()
plot(t,RMSE_y1,'r-',t,RMSE_y2,'b-',t,RMSE_y3,'g-',t,RMSE_y4,'k-','linewidth',1.5);
legend('RMSE(EKF y)','RMSE(UKF y)','RMSE(CKF3 y)','RMSE(CKF5 y)');
title('East position RMSE');
xlabel('t(s)');ylabel('y-RMSE(m)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([0,5]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图   
% plot(t,RMSE_y1,'r-',t,RMSE_y2,'b--',t,RMSE_y3,'g--',t,RMSE_y4,'k--','linewidth',1); % 绘制局部曲线图                        
% xlim([250,350]); % 设置坐标轴范围  
% ylim([0,0.008]);

figure()
plot(t,RMSE_U1,'r-',t,RMSE_U2,'b-',t,RMSE_U3,'g-',t,RMSE_U4,'k-','linewidth',1.5);
legend('RMSE(EKF U)','RMSE(UKF U)','RMSE(CKF3 U)','RMSE(CKF5 U)');
title('course speed RMSE');
xlabel('t(s)');ylabel('U-RMSE(m/s)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([0,0.1]);
% axes('Position',[0.5,0.25,0.3,0.3]); % 生成子图   
% plot(t,RMSE_U1,'r-',t,RMSE_U2,'b-',t,RMSE_U3,'g-',t,RMSE_U4,'k-','linewidth',1.5); % 绘制局部曲线图                        
% xlim([400,500]); % 设置坐标轴范围  
% ylim([0,0.001]);

figure()
plot(t,RMSE_X1,'r-',t,RMSE_X2,'b-',t,RMSE_X3,'g-',t,RMSE_X4,'k-','linewidth',1.5);
legend('RMSE(EKF X)','RMSE(UKF X)','RMSE(CKF3 X)','RMSE(CKF5 X)');
title('course angle RMSE');
xlabel('t(s)');ylabel('X-RMSE(rad)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([0,0.03]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图   
% plot(t,RMSE_X1,'r-',t,RMSE_X2,'b--',t,RMSE_X3,'g--',t,RMSE_X4,'k-','linewidth',1); % 绘制局部曲线图                        
% xlim([250,350]); % 设置坐标轴范围  
% ylim([0.025,0.035]);

figure()
plot(t,RMSE_Wx1,'r-',t,RMSE_Wx2,'b-',t,RMSE_Wx3,'g-',t,RMSE_Wx4,'k-','linewidth',1.5);
legend('RMSE(EKF WX)','RMSE(UKF Wx)','RMSE(CKF3 Wx)','RMSE(CKF5 Wx)');
title('course anguler rate RMSE');
xlabel('t(s)');ylabel('Wx-RMSE(rad/s)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([0,0.04]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图   
% plot(t,RMSE_Wx1,'r-',t,RMSE_Wx2,'b-',t,RMSE_Wx3,'g--',t,RMSE_Wx4,'k-','linewidth',1); % 绘制局部曲线图                        
% xlim([250,350]); % 设置坐标轴范围  
% ylim([0.005,0.01]);

figure()
plot(y5,x5,'c-','linewidth',1.5);hold on
plot(y1,x1,'r-','linewidth',1.5);
plot(y2,x2,'g-','linewidth',1.5);
plot(y3,x3,'b-','linewidth',1.5);
plot(y4,x4,'k-','linewidth',1.5);
legend('real value','EKF','UKF','CKF3','CKF5')
xlabel('E(m)'),ylabel('N(m)');
point_database =[300 300; 1500 300; 2900 1900; 2900 3900; 1500 5500; 300 5500]';

% plot(point_database(2,:),point_database(1,:),'b-','linewidth',1.5)
for k = 1:length(point_database(1,:))
    plot(point_database(2,k),point_database(1,k),'ro','linewidth',3);
end
xlim([0,5600]);
ylim([-100,3000]);

axes('Position',[0.3,0.2,0.3,0.3]); % 生成子图   
plot(y5,x5,'c-','linewidth',1.5);hold on
plot(y1,x1,'r-','linewidth',1.5);
plot(y2,x2,'g-','linewidth',1.5);
plot(y3,x3,'b--','linewidth',1.5);
plot(y4,x4,'k-','linewidth',1.5);
xlim([2280,2330]); % 设置坐标轴范围  
ylim([2880,2910]);hold off
hold off

% course angle estimate 
figure()
plot(t,X5,'c-',t,X1,'r-',t,X2,'g-',t,X3,'b-',t,X4,'k-','linewidth',1.5);
legend('real value','EKF','UKF','CKF3','CKF5');
title('course angle estimate');
xlabel('t(s)');ylabel('course angle(rad)');
% axes('Position',[0.5,0.2,0.2,0.2]); % 生成子图   
% plot(t,X5,'c-',t,X1,'r-',t,X2,'g-',t,X3,'b-',t,X4,'k-','linewidth',1.5);
% xlim([700,750]); % 设置坐标轴范围  
% ylim([1.56,1.58]);
% course anduler rate estimate value
figure()
plot(t,Wx5,'c-',t,Wx1,'r-',t,Wx2,'g-',t,Wx3,'b-',t,Wx4,'k-','linewidth',1.5);
legend('real value','EKF','UKF','CKF3','CKF5');
title('course anguler rate estimate');
xlabel('t(s)');ylabel('course anguler rate(rad/s)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([-0.03,0.02]);
axes('Position',[0.4,0.2,0.2,0.2]); % 生成子图   
plot(t,Wx5,'c-',t,Wx1,'r-',t,Wx2,'g-',t,Wx3,'b-',t,Wx4,'k-','linewidth',1.5);
xlim([650,700]); % 设置坐标轴范围  
ylim([-0.00005,0.00005]);

% disp

disp('祝老师身体健康，学术长青！')
