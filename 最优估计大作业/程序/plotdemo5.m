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
plot(y,x,'r-','linewidth',1);
point_database =[300 300; 1500 300; 2900 1900; 2900 3900; 1500 5500; 300 5500]';
% plot(point_database(2,:),point_database(1,:),'b-','linewidth',1)

for k = 1:length(point_database(1,:))
    plot(point_database(2,k),point_database(1,k),'ro','linewidth',4);
end
title('Map-XY');
xlabel('E(m)');ylabel('N(m)');
hold off
figure
plot(t,x,'r',t,xm,'g',t,xhat,'b','linewidth',0.5);
legend('x','xm','xhat');
title('����λ�ù��ƶԱ�');
xlabel('ʱ��(s)');ylabel('����λ��(m)');
figure
plot(t,y,'r',t,ym,'g',t,yhat,'b','linewidth',0.5);
title('����λ�ù��ƶԱ�');
xlabel('ʱ��(s)');ylabel('����λ��(m)');
legend('y','ym','yhat');
figure
plot(t,U,'r',t,Uhat,'b','linewidth',0.5);
title('�����ٶȹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('�����ٶ�(m/s)');
legend('U','Uhat');
figure
plot(t,delta*180/pi,'r','linewidth',0.5);
title('���');
xlabel('ʱ��(s)');ylabel('���(deg)');
figure
plot(t,X*180/pi,'r',t,Xhat*180/pi,'b','linewidth',0.5);
legend('X measure','X estimate');
title('����ǹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('�����(deg)');
legend('X','Xhat');
figure
plot(t,Wx,'r',t,Wxhat,'b','linewidth',0.5);
legend('Wx','Wxhat');
title('����Ǳ仯�ʹ��ƶԱ�');
xlabel('ʱ��(s)');ylabel('���������(rad/s)');
legend('Wx','Wxhat');

