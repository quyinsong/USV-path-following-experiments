function  modelplot( pos, psai )
% MODELPLOT modelplot(pos,psai) :plot the geometrical model of ship in real time
%
% Input:  pos; ship real time position in {n}, pos=[x y]'; 
%         psai: real time yaw angle
% Output: none
%
% Auther: Quyinsong
% Data:   12nd Jan 2022

%% check of input dimensions
if nargin~=2 
    error('the number of input must be 2');end
if length(pos)~=2 
    error('the number of pos must be 2');end
if length(psai)~=1 
    error('the number of psai must be 1');end
%% ship geometrical shape is characterized by five point:
% xb1,xb2,xb3,xb4,xb5, their coodinates are described in {b} 
% ship  length: 4.5m   width: 1.5m
L = 1.65; W = 0.75;
xb1=[0.35 -0.375]';xb2=[1.3 0]';xb3=[0.35 0.375]';xb4=[-0.35 0.375]';xb5=[-0.35 -0.375]';
%% rotation matrix from {b} to {n}
Rb_n=[cos(psai) -sin(psai);
      sin(psai) cos(psai)]; 
%% trasfer the five point in {b} to {n}
xn1=Rb_n*xb1+pos;
xn2=Rb_n*xb2+pos;
xn3=Rb_n*xb3+pos;
xn4=Rb_n*xb4+pos;
xn5=Rb_n*xb5+pos;
%% plot
figure(1);hold on
N=pos(1); E=pos(2);
plot(E,N,'r-','linewidth',2);  % plot the navagation line
plot([xn1(2) xn2(2)],[xn1(1) xn2(1)],'g-',[xn2(2) xn3(2)],[xn2(1) xn3(1)],'g-',...
[xn3(2) xn4(2)],[xn3(1) xn4(1)],'g-',[xn4(2) xn5(2)],[xn4(1) xn5(1)],'g-',...
[xn5(2) xn1(2)],[xn5(1) xn1(1)],'g-','linewidth',0.5); % plot ship model
XX = [xn1(1) xn2(1) xn3(1) xn4(1) xn5(1) ]; YY = [xn1(2) xn2(2) xn3(2) xn4(2) xn5(2) ];
fill(YY,XX,'g');
set(gca,'xTick',-20:10:100);
set(gca,'yTick',-20:10:100);
axis([-20 100,-20 100]);
grid on;
xlabel('E/¶«');ylabel('N/±±');
end

