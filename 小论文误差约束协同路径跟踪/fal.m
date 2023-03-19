function y = fal( x, a, d )
%FAL 此处显示有关此函数的摘要
%   此处显示详细说明

if (abs(x)<=d)
   y = x/d^(1-a); 
else
   y = abs(x)^a*sign(x);
end

end

