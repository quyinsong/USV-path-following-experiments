% KF ���ؿ������ʵ��
% ����LOSdemoKF.slx�ļ�
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
title('����ǹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('�����(deg)');
xlim([0,1200]); % ���������᷶Χ  
ylim([-15,120]);
axes('Position',[0.5,0.17,0.2,0.2]); % ������ͼ   
plot(t,psim*180/pi,'g-','linewidth',0.1); hold on% ���ƾֲ�����ͼ   
plot(t,psihat*180/pi,'r-',t,psi*180/pi,'b--','linewidth',2.5); hold off% ���ƾֲ�����ͼ   
xlim([400,500]); % ���������᷶Χ  
ylim([20,50]);

figure
plot(t,rhat,'r-','linewidth',1.5);hold on
plot(t,r,'b--','linewidth',2);hold off
legend('rhat','r');
title('ת���ʹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('���������(rad/s)');
xlim([0,1200]); % ���������᷶Χ  
ylim([-0.02,0.04]);
% axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ   
% plot(t,r,'b',t,rhat,'r','linewidth',0.5); % ���ƾֲ�����ͼ                        
% xlim([300,400]); % ���������᷶Χ  
% ylim([-1e-2,1e-2]);

figure
plot(t,ksi2hat,'r-','linewidth',0.5);hold on 
plot(t,ksi2,'b--','linewidth',0.5); hold off
legend('ksi2hat','ksi2');
title('���˹��ƶԱ�');
xlabel('ʱ��(s)');ylabel('���������(rad/s)');
xlim([0,100]); % ���������᷶Χ  
ylim([-0.02,0.04]);
% axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ   
% plot(t,r,'b',t,rhat,'r','linewidth',0.5); % ���ƾֲ�����ͼ                        
% xlim([300,400]); % ���������᷶Χ  
% ylim([-1e-2,1e-2]);

figure()
plot(t,RMSE_ksi2,'r-','linewidth',1);
title('wave motion RMSE');
xlim([1,1200]); % ���������᷶Χ  
ylim([0,0.008]);
xlabel('t(s)');ylabel('wave motion RMSE(rad^2)');
% axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ   
% plot(t,RMSE_ksi2,'r-','linewidth',1);
% xlim([1,5]); % ���������᷶Χ  
% ylim([0,0.008]);

figure()
plot(t,RMSE_psi,'r-','linewidth',0.5);
title('heading anguler rate RMSE');
xlabel('t(s)');ylabel('heading anguler rate(deg)');
xlim([1,1200]); % ���������᷶Χ  
ylim([0,10]);
% axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ  
% plot(t,RMSE_psi,'r-','linewidth',0.5);

figure()
plot(t,RMSE_r,'r-','linewidth',0.5);
title('heading anguler rate RMSE');
xlabel('t(s)');ylabel('heading anguler rate(deg/s)');
xlim([1,1200]); % ���������᷶Χ  
ylim([0,0.015]);
% axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ  
% plot(t,RMSE_r,'r-','linewidth',0.5);
% xlim([1,5]); % ���������᷶Χ  
% ylim([0,0.008]);

figure
plot(t,delta*180/pi,'r','linewidth',0.5);
title('���');
xlabel('ʱ��(s)');ylabel('���(deg)');


% disp

disp('ף��ʦ���彡����ѧ�����࣡')


