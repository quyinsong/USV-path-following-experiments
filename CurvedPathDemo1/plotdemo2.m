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
point_database =PathPoint';
plot(point_database(2,:),point_database(1,:),'b-','linewidth',1)
hold off
figure
plot(t,psi*180/pi,'r','linewidth',2);
xlabel('time/s');ylabel('psai/deg');
figure
plot(t,u,'r','linewidth',2)
xlabel('time/s');ylabel('u (m/s)');
figure
plot(t,Tu,'r',t,Tr,'b','linewidth',2)
legend('surge force','yaw torch');
