function  miu1_hat  = DynamicSystem_( miu1_hat, miu1, Hx_hat, Q_hat,tau ,ts, B ,Fx)
%DYNAMICSYSTEM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
Gamma = diag([0 0]);
miu1_hat_dot = Gamma*(miu1-miu1_hat)+Fx+tau-B*Hx_hat+Q_hat;
miu1_hat = miu1_hat+miu1_hat_dot*ts;
end

