function [sys,x0,str,ts,simStateCompliance] = KF4(t,x,u,flag,Q,R)
% reference: Fossen 2011 P303
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

DiscStatesNum = 20;
OutputsNum = 20;
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
Xk0 = [0 0 0 0]';
Pk0 = diag([1 1 1 1]);
Pk0 = reshape(Pk0,[16,1]);
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
Yk = u(2);
Uk = -u(1);
% Discrect State
ksi = [x(1);x(2)];
psi = x(3);
r = x(4);
Xk = [ksi;psi;r];
Pk = x(5:20);
Pk = reshape(Pk,[4,4]);
% continue system matrix
w0 = 1.2; lamda = 0.1; K = 0.185; T = 107.3;
A = [0 1 0 0 ;
    -w0^2 -2*lamda*w0 0 0 ;
    0 0 0 1 ;
    0 0 0 -1/T;]; 
Ndelta = K/T;
B = [0;0;0;Ndelta];
E = [0 0;1 0;0 0;0 1];
% dicrect system matrix
I4 = diag(ones(1,4));
Faik = I4+h*A;
Bk = h*B;
Tk = h*E;
Hk = [0 1 1 0];
% noise
Qk = Q; % process noise
Rk = R; % measurement noise
% Kalman Gain
Kk = Pk*Hk'*inv(Hk*Pk*Hk'+Rk);
% Corrector
Pkhat = (I4-Kk*Hk)*Pk*(I4-Kk*Hk)'+Kk*Rk*Kk';
Xkhat = Xk+Kk*(Yk-Hk*Xk);
% predictor
Xk1 = Faik*Xkhat+Bk*Uk;
Pk1 = Faik*Pkhat*Faik'+Tk*Qk*Tk';
Pk1 = reshape(Pk1,[16,1]);
sys = [Xk1;Pk1];
% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
Xk = [x(1);x(2);x(3);x(4)];
Pk = x(5:20);
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
