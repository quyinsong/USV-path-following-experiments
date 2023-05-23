% KF 蒙特卡洛仿真实验
% 调用LOSdemoKF.slx文件
clc
close all
% KF
Mc = 1;
disp('KF')
sumksi2 = 0;
sumpsi = 0;
sumr = 0;
tic;
for kc = 1:Mc
    disp(kc);
    sim LOSdemoKF;
    sumksi2 = sumksi2+((ksi2-ksi2hat)).^2;
    sumpsi = sumpsi+((psi-psihat)*180/pi).^2;
    sumr = sumr+(r-rhat).^2;
end
TKF = toc;
tkf = TKF/Mc;
    RMSE_ksi2 = sqrt(sumksi2/Mc);
    RMSE_psi = sqrt(sumpsi/Mc);
    RMSE_r = sqrt(sumr/Mc);
    mean_RMSE_ksi2 = mean(RMSE_ksi2);
    mean_RMSE_psi = mean(RMSE_sumpsi);
    mean_RMSE_r = mean(RMSE_sumr);

% disp
disp('KF mean RMSE and var RMSE : ');
disp(['mean_RMSE_ksi2 : ', num2str(mean_RMSE_ksi2)]);
disp(['mean_RMSE_psi : ', num2str(mean_RMSE_psi)]);
disp(['mean_RMSE_r : ', num2str(mean_RMSE_r)]);



% plot

figure()
plot(t,psim*180/pi,'g-','linewidth',0.1);hold on
plot(t,psihat*180/pi,'r-','linewidth',1.5);
plot(t,psi*180/pi,'b--','linewidth',2);hold off
legend('psi measure','psihat','psi ');
title('艏向角估计对比');
xlabel('时间(s)');ylabel('艏向角(deg)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([-15,120]);
axes('Position',[0.5,0.17,0.2,0.2]); % 生成子图   
plot(t,psim*180/pi,'g-','linewidth',0.1); hold on% 绘制局部曲线图   
plot(t,psihat*180/pi,'r-',t,psi*180/pi,'b--','linewidth',2.5); hold off% 绘制局部曲线图   
xlim([400,500]); % 设置坐标轴范围  
ylim([20,50]);

figure
plot(t,rhat,'r-','linewidth',1.5);hold on
plot(t,r,'b--','linewidth',2);hold off
legend('rhat','r');
title('转艏率估计对比');
xlabel('时间(s)');ylabel('艏向角速率(rad/s)');
xlim([0,1200]); % 设置坐标轴范围  
ylim([-0.02,0.04]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图   
% plot(t,r,'b',t,rhat,'r','linewidth',0.5); % 绘制局部曲线图                        
% xlim([300,400]); % 设置坐标轴范围  
% ylim([-1e-2,1e-2]);

figure
plot(t,ksi2hat,'r-','linewidth',0.5);hold on 
plot(t,ksi2,'b--','linewidth',0.5); hold off
legend('ksi2hat','ksi2');
title('海浪估计对比');
xlabel('时间(s)');ylabel('艏向角速率(rad/s)');
xlim([0,100]); % 设置坐标轴范围  
ylim([-0.02,0.04]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图   
% plot(t,r,'b',t,rhat,'r','linewidth',0.5); % 绘制局部曲线图                        
% xlim([300,400]); % 设置坐标轴范围  
% ylim([-1e-2,1e-2]);

figure()
plot(t,RMSE_ksi2,'r-','linewidth',1);
title('wave motion RMSE');
xlim([1,1200]); % 设置坐标轴范围  
ylim([0,0.008]);
xlabel('t(s)');ylabel('wave motion RMSE(rad^2)');
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图   
% plot(t,RMSE_ksi2,'r-','linewidth',1);
% xlim([1,5]); % 设置坐标轴范围  
% ylim([0,0.008]);

figure()
plot(t,RMSE_psi,'r-','linewidth',0.5);
title('heading anguler rate RMSE');
xlabel('t(s)');ylabel('heading anguler rate(deg)');
xlim([1,1200]); % 设置坐标轴范围  
ylim([0,10]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图  
% plot(t,RMSE_psi,'r-','linewidth',0.5);

figure()
plot(t,RMSE_r,'r-','linewidth',0.5);
title('heading anguler rate RMSE');
xlabel('t(s)');ylabel('heading anguler rate(deg/s)');
xlim([1,1200]); % 设置坐标轴范围  
ylim([0,0.015]);
% axes('Position',[0.2,0.5,0.3,0.3]); % 生成子图  
% plot(t,RMSE_r,'r-','linewidth',0.5);
% xlim([1,5]); % 设置坐标轴范围  
% ylim([0,0.008]);

figure
plot(t,delta*180/pi,'r','linewidth',0.5);
title('舵角');
xlabel('时间(s)');ylabel('舵角(deg)');


% disp

disp('祝老师身体健康，学术长青！')


