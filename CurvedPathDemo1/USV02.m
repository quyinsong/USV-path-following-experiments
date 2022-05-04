function [sys,x0,str,ts,simStateCompliance] = USV02(t,x,u,flag)
switch flag
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
  case 1
    sys=mdlDerivatives(t,x,u);
  case {2,4,9}
    sys=[];
  case 3
    sys=mdlOutputs(t,x,u);
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
sizes = simsizes;

sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 7;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [0 0 0 20 5 0];

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

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
% state
% x1=x(1): surge velocity (m/s)
% x2=x(2): sway velocity (m/s)
% x3=x(3): yaw velocity (rad/s)
% x4=x(4): position in {n} x-direction (m)
% x5=x(5): position in {n} y-direction (m)
% x6=x(6): yaw angel in {n} (rad)
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5); x6 = x(6); r = x3; v = x2; 
U = [x1;x2;x3]; P = [x4;x5;x6];
Rbn=[ cos(x6) -sin(x6) 0;
      sin(x6) cos(x6)  0;
      0          0     1];
% disturbance
d = [u(1);u(2);u(3)];
Vc= u(4);
betac= u(5);
UC_I = [Vc*cos(betac) Vc*sin(betac) 0]'; % current velocity in {n} coordinate
UC_B = Rbn'*UC_I; % current velocity in {b} coordinate
ucb = UC_B(1);  vcb = UC_B(2); 
UC_B_dot = [r*vcb -r*ucb 0]'; 
ur = x1-ucb;  vr = x2-vcb;
Ur = [ur vr x3]';  % relative speed in {b}
% control input
Tu = u(6); Tr = u(7);
% inertial parameters
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
m11 = 25.8; m22 = 33.8; m23 = 1.0948; m32 = m23;  m33 = 2.76; m0 = 92.0894;
% hydrodynamic parameters
Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
Xuuu=-5.86643;       Yvv=-36.47287;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;    
% matrix
c13 = -m*(xg*r+v); c23 = m*x1;
c31 = -c13; c32 = -c23;            % C 参数
d11=-Xu-Xuu*abs(ur)-Xuuu*ur^2;
d22=-Yv-Yvv*abs(vr)-Yrv*abs(r);
d23=-Yr-Yvr*abs(vr)-Yrr*abs(r);
d32=-Nv-Nvv*abs(vr)-Nrv*abs(r);
d33=-Nr-Nvr*abs(vr)-Nrr*abs(r);    % D 参数
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33]; % M = MA+MRB
MA = [-Xudot   0        0;
       0     -Yvdot  -Yrdot;
       0     -Nvdot  -Nrdot];
CRB=[0     0     c13;
     0     0     c23;
    c31   c32     0 ];
CA = [0         0     Yvdot*vr;
      0         0    -Xudot*ur;
    -Yvdot*vr  Xudot*ur   0 ] ;
C = CA+CRB; % C = CA+CRB
D=  [d11  0     0;
     0    d22  d23;
     0    d32  d33];
% output
tao = [Tu;0;Tr];
Udot=M\(-C*U+MA*UC_B_dot+CA*UC_B-D*Ur+tao+d) ;
Pdot=Rbn*U+UC_I;
xdot=[Udot ;Pdot];
sys=xdot;
% end mdlDerivatives

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
sys=[x(1);x(2);x(3);x(4);x(5);x(6)];
% end mdlOutputs

