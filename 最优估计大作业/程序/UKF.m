function [sys,x0,str,ts,simStateCompliance] = UKF(t,x,u,flag,Q,R)
% Date: 2022.5.1
% Author: QuYinsong
switch flag

  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  case 1
    sys=mdlDerivatives(t,x,u);

  case 2
    sys=mdlUpdate(t,x,u,Q,R);

  case 3
    sys=mdlOutputs(t,x,u);

  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  case 9
    sys=mdlTerminate(t,x,u);

  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%

DiscStatesNum = 30;
OutputsNum = 30;
InputsNum = 2;
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = DiscStatesNum;
sizes.NumOutputs     = OutputsNum;
sizes.NumInputs      = InputsNum;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
Xk0 = [105 305 6 0.1 0.01]';
Pk0 = diag([100 100 1 1 1e-2]);
Pk0 = reshape(Pk0,[25,1]);
x0  = [Xk0;Pk0];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

sys = [];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u,Q,R)
h = 0.1;
% Input
xm = u(1);
ym = u(2);
Yk = [xm ym]';
% Discrect State
xx = x(1);
yy = x(2);
U = x(3);
X = x(4);
Wx = x(5);
Xk = [xx;yy;U;X;Wx];
Pk = x(6:30);
Pk = reshape(Pk,[5,5]);
% system parameters
afa1 = 1e-5; afa2 = 0.2;
ns = 5; % states number
% noise
Qk = Q; % process noise
Rk = R; % measurement noise

%% predictor
% get sigma points
k1 = 3-ns; afa3 = 0.5; lamda = afa3^2*(ns+k1)-ns; beta = 2;
Wim0 = lamda/(ns+lamda);
Wim1 = 1/(2*(ns+lamda));
Wic0 = lamda/(ns+lamda)+1+beta-afa3^2;
Wic1 = 1/(2*(ns+lamda));
Wim = [Wim0 repmat(Wim1,[1 2*ns])];
Wic = [Wic0 repmat(Wic1,[1 2*ns])];
Sk1k1=utchol(Pk);                                %%%产生方根矩阵
CPtArray = sqrt(ns+lamda)*[eye(ns) -eye(ns)];        %%%发生器
Xki1 = Xk+Sk1k1*CPtArray;

Xki11 = [Xk Xki1];
% substitute sigma points in states equation
Xki1_ = [Xki11(1,:)+h*Xki11(3,:).*cos(Xki11(4,:));
        Xki11(2,:)+h*Xki11(3,:).*sin(Xki11(4,:));
        (1-h*afa1)*Xki11(3,:);
        Xki11(4,:)+h*Xki11(5,:);
        (1-h*afa2)*Xki11(5,:)]; %  5 row 10 colum
% get  predictive states
Xk1 = Xki1_*Wim';
% get predictive convariance matrix
Ek = [0 0;0 0;h 0;0 0; 0 h];
Hk = [1 0 0 0 0;
      0 1 0 0 0];
Sum1 = 0;
for i = 1:2*ns+1
    Sum1 = Sum1+(Xki1_(:,i)-Xk1)*(Xki1_(:,i)-Xk1)'*Wic(i);
end
Pk1 = Sum1+Ek*Qk*Ek';
%% Corrector
% get sigma points
Skk1 = utchol(Pk1);

Xki2_ = Xk1+Skk1*CPtArray;

Xki22_ = [Xk1 Xki2_];
Yki_ = Hk*Xki22_;
Yk1 = Yki_*Wim';
% get measurement convariance matrix
Sum2 = 0; Sum3 = 0;
for i = 1:2*ns+1
    Sum2 = Sum2+(Yki_(:,i)-Yk1)*(Yki_(:,i)-Yk1)'*Wic(i);
end
Pzz = Sum2+Rk;
% get Xk1 and Yk1 convariance matrix
for i = 1:2*ns+1
    Sum3 = Sum3+(Xki22_(:,i)-Xk1)*(Yki_(:,i)-Yk1)'*Wic(i);
end
Pxz = Sum3;
% Kalman Gain
Kk = Pxz*inv(Pzz);
% corrector
Xk = Xk1+Kk*(Yk-Yk1);
Pk = Pk1-Kk*Pzz*Kk';
Pk = reshape(Pk,[25,1]);
sys = [Xk;Pk];
% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
Xk = [x(1);x(2);x(3);x(4);x(5)];
Pk = x(6:30);
sys = [Xk;Pk];

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 0.01;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
