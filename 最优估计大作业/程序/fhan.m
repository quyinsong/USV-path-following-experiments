function fh = fhan( x1,x2,r,h )
%FHAN 
% Date : 11st March 2022
% 跟踪微分器（TD2）实验
% 参考文献：自抗扰控制技术  作者：韩京清
if r<=0
    error('deta must be > 0!');
end
d = r*h^2; a0 = h*x2;
y = x1+a0;
a1 = sqrt(d*(d*8*abs(y)));
a2 = a0+sign(y)*(a1-d)/2;
sy = (sign(y+d)-sign(y-d))/2;
a = (a0+y-a2)*sy+a2;
sa = (sign(a+d)-sign(a-d))/2;
fh = -r*(a/d-sign(a))*sa-r*sign(a);

end

