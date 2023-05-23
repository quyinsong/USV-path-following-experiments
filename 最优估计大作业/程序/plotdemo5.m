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
title('北向位置估计对比');
xlabel('时间(s)');ylabel('北向位置(m)');
figure
plot(t,y,'r',t,ym,'g',t,yhat,'b','linewidth',0.5);
title('东向位置估计对比');
xlabel('时间(s)');ylabel('东向位置(m)');
legend('y','ym','yhat');
figure
plot(t,U,'r',t,Uhat,'b','linewidth',0.5);
title('航行速度估计对比');
xlabel('时间(s)');ylabel('航行速度(m/s)');
legend('U','Uhat');
figure
plot(t,delta*180/pi,'r','linewidth',0.5);
title('舵角');
xlabel('时间(s)');ylabel('舵角(deg)');
figure
plot(t,X*180/pi,'r',t,Xhat*180/pi,'b','linewidth',0.5);
legend('X measure','X estimate');
title('航向角估计对比');
xlabel('时间(s)');ylabel('航向角(deg)');
legend('X','Xhat');
figure
plot(t,Wx,'r',t,Wxhat,'b','linewidth',0.5);
legend('Wx','Wxhat');
title('航向角变化率估计对比');
xlabel('时间(s)');ylabel('航向角速率(rad/s)');
legend('Wx','Wxhat');

