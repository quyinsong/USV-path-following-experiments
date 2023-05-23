function [sys,x0,str,ts,simStateCompliance] = CKF5(t,x,u,flag,Q,R)
% Date: 2022.5.1
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
% measurement matrix and process noise matrix  
Ek = [0 0;0 0;h 0;0 0; 0 h];
Hk = [1 0 0 0 0;
      0 1 0 0 0];
%% predictor
Sk1k1=utchol(Pk);                                %%%产生方根矩阵

xk1k1 = Xk;
E=[eye(ns) -eye(ns)];
% S = [0.7071    0.7071    0.7071    0.7071         0         0         0         0         0    0;
%      0.7071         0         0         0    0.7071    0.7071    0.7071         0         0    0;
%      0    0.7071         0         0    0.7071         0         0    0.7071    0.7071         0;
%      0         0    0.7071         0         0    0.7071         0    0.7071         0    0.7071;
%      0         0         0    0.7071         0         0    0.7071         0    0.7071    0.7071];
% S_ = [0.7071    0.7071    0.7071    0.7071         0         0         0         0         0    0;
%      -1.0000         0         0         0    0.7071    0.7071    0.7071         0         0    0;
%       0   -1.0000         0         0   -1.0000         0         0    0.7071    0.7071         0;
%       0         0   -1.0000         0         0   -1.0000         0   -1.0000         0    0.7071;
%       0         0         0   -1.0000         0         0   -1.0000         0   -1.0000   -1.0000];

S = [0.7071    0.7071         0    0.7071         0         0    0.7071         0         0         0;
     0.7071         0    0.7071         0    0.7071         0         0    0.7071         0         0;
          0    0.7071    0.7071         0         0    0.7071         0         0    0.7071         0;
          0         0         0    0.7071    0.7071    0.7071         0         0         0    0.7071;
          0         0         0         0         0         0    0.7071    0.7071    0.7071    0.7071];
S_ = [-0.7071   -0.7071         0   -0.7071         0         0   -0.7071         0         0         0;
       0.7071         0   -0.7071         0   -0.7071         0         0   -0.7071         0         0;
            0    0.7071    0.7071         0         0   -0.7071         0         0   -0.7071         0;
            0         0         0    0.7071    0.7071    0.7071         0         0         0   -0.7071;
            0         0         0         0         0         0    0.7071    0.7071    0.7071    0.7071];
CPtArrayE = sqrt(ns+2)*E;        %%%发生器 5x10
CPtArrayS = sqrt(ns+2)*[S -S];        %%%发生器 5x20
CPtArrayS_ = sqrt(ns+2)*[S_ -S_];        %%%发生器 5x20

XEk1k1=repmat(xk1k1,1,2*ns) + Sk1k1*CPtArrayE;      %%%产生容积点
XSk1k1=repmat(xk1k1,1,4*ns) + Sk1k1*CPtArrayS;      %%%产生容积点
XS_k1k1=repmat(xk1k1,1,4*ns) + Sk1k1*CPtArrayS_;      %%%产生容积点
% 总容积点数为: 2nn+1=2*5*5+1=51
Xk1k1 = [xk1k1 XEk1k1 XSk1k1 XS_k1k1];  % 5x51
%%%传播容积点  
Xkk1 = [Xk1k1(1,:)+h*Xk1k1(3,:).*cos(Xk1k1(4,:));
        Xk1k1(2,:)+h*Xk1k1(3,:).*sin(Xk1k1(4,:));
        (1-h*afa1)*Xk1k1(3,:);
        Xk1k1(4,:)+h*Xk1k1(5,:);
        (1-h*afa2)*Xk1k1(5,:)];   % 5x51

% Xkk1 = [xk1k1(1,:)+h*Xk1k1(3,:).*cos(Xk1k1(4,:));
%         xk1k1(2,:)+h*Xk1k1(3,:).*sin(Xk1k1(4,:));
%         xk1k1(3,:)-h*afa1*Xk1k1(4,:);
%         xk1k1(4,:)+h*Xk1k1(5,:);
%         xk1k1(5,:)-h*afa2*Xk1k1(5,:)]; %%%传播容积点   

xkk1=2*Xkk1(:,1)/7-sum(Xkk1(:,2:11),2)/98+sum(Xkk1(:,12:31),2)/49+sum(Xkk1(:,32:51),2)/49;       

Pkk1=2*Xkk1(:,1)*Xkk1(:,1)'/7-Xkk1(:,2:11)*Xkk1(:,2:11)'/98+Xkk1(:,12:31)*Xkk1(:,12:31)'/49+...
     Xkk1(:,32:51)*Xkk1(:,32:51)'/49-xkk1*xkk1'+Ek*Qk*Ek';             


%% Corrector

Skk1=utchol(Pkk1);
Xi0 = xkk1;
XiE=  repmat(xkk1,1,2*ns) + Skk1*CPtArrayE;
XiS=  repmat(xkk1,1,4*ns) + Skk1*CPtArrayS;
XiS_=  repmat(xkk1,1,4*ns) + Skk1*CPtArrayS_;
Xi = [Xi0 XiE XiS XiS_];
Zi =Hk*Xi;

zkk1 = 2*Zi(:,1)/7-sum(Zi(:,2:11),2)/98+sum(Zi(:,12:31),2)/49+sum(Zi(:,32:51),2)/49;       

Pzzkk1 = 2*Zi(:,1)*Zi(:,1)'/7-Zi(:,2:11)*Zi(:,2:11)'/98+Zi(:,12:31)*Zi(:,12:31)'/49+...
         Zi(:,32:51)*Zi(:,32:51)'/49-zkk1*zkk1'+Rk;         % 5x5    
% Pzzkk1=Zi*Zi'/10-zkk1*zkk1'+Rk;             


Pxzkk1 = 2*Xi(:,1)*Zi(:,1)'/7-Xi(:,2:11)*Zi(:,2:11)'/98+Xi(:,12:31)*Zi(:,12:31)'/49+...
         Xi(:,32:51)*Zi(:,32:51)'/49-xkk1*zkk1';         % 5x2     
% Pxzkk1 = Xkk1*Zi'/10-xkk1*zkk1';         % 5x2   

Wk=Pxzkk1*inv(Pzzkk1);     % 5x2

Xk = xkk1 + Wk*(Yk - zkk1);

Pk=Pkk1-Wk*Pzzkk1*Wk';   

Pk = reshape(Pk,[25 1]);
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
