function y = fal( x, a, d )
%FAL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

if (abs(x)<=d)
   y = x/d^(1-a); 
else
   y = abs(x)^a*sign(x);
end

end

