function  plotCircle( R, pos, color, linewid )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

alpha=0:pi/20:2*pi;%�Ƕ�[0,2*pi]
x=R*cos(alpha)+pos(1);
y=R*sin(alpha)+pos(2);
plot(x,y,color,'linewidth',linewid)
% fill(x,y,color);%���
 
end

