function  plotCircle( R, pos, color, linewid )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

alpha=0:pi/20:2*pi;%角度[0,2*pi]
x=R*cos(alpha)+pos(1);
y=R*sin(alpha)+pos(2);
plot(x,y,color,'linewidth',linewid)
% fill(x,y,color);%填充
 
end

