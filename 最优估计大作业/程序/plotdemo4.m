clc
close all
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
plot(y,x,'r-','linewidth',1)
point_database =[300 300; 2900 2500; 2900 4800]';
plot(point_database(2,:),point_database(1,:),'b-','linewidth',1)
title('Map-XY');
hold off
figure
plot(t,delta*180/pi,'r','linewidth',0.5);
title('���');
xlabel('ʱ��(s)');ylabel('���(deg)');
figure
plot(t,psim*180/pi,'g',t,psi*180/pi,'b',t,psihat*180/pi,'r','linewidth',0.5);
legend('psi measure','psi','psi estimate');
title('����ǹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('�����(deg)');
axes('Position',[0.2,0.17,0.3,0.3]); % ������ͼ   
plot(t,psim*180/pi,'g',t,psi*180/pi,'b',t,psihat*180/pi,'r','linewidth',0.5);% ���ƾֲ�����ͼ                        
xlim([400,500]); % ���������᷶Χ  
ylim([20,50]);
figure
plot(t,r,'b',t,rhat,'r','linewidth',0.5);
legend('r','rhat');
title('ת���ʹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('���������(rad/s)');
axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ   
plot(t,r,'b',t,rhat,'r','linewidth',0.5); % ���ƾֲ�����ͼ                        
xlim([300,400]); % ���������᷶Χ  
ylim([-1e-2,1e-2]);
figure
plot(t,ksi2,'b',t,ksi2hat,'r','linewidth',0.5);
legend('ksi2','ksi2hat');
title('��Ƶ�˶����ƶԱ�');
xlabel('ʱ��(s)');ylabel('��Ƶ�˶�(rad/s)');
axes('Position',[0.2,0.5,0.3,0.3]); % ������ͼ   
plot(t,ksi2,'b',t,ksi2hat,'r','linewidth',0.5); % ���ƾֲ�����ͼ                        
xlim([400,500]); % ���������᷶Χ  
ylim([-0.05,0.05]);

