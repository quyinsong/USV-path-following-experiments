function [sys,x0,str,ts,simStateCompliance] = LOSguidance(t,x,u,flag,PointSet)
switch flag
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
  case 1
    sys=mdlDerivatives(t,x,u);
  case {2,4,9}
    sys=[];
  case 3
    sys=mdlOutputs(t,x,u,PointSet);
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

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

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
sys = [];
% end mdlDerivatives

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,PointSet)
% input
pointer = u(1);
x4 = u(2); x5 = u(3); 
deta = u(4);% look-ahead distance
R = u(5); % change radius
% PointSet =[0 0 2900 2900  2900 4900]';
WayPoint_database = reshape(PointSet,[2,length(PointSet)/2]); % WayPoint_database =[0 0 ;2900 2900;  2900 4900]';
% get the first piece path discribed by two points
Pk = WayPoint_database(:,pointer);
Pk1 = WayPoint_database(:,pointer+1);
% cal the angle between current path and N direction 
afak=atan2(Pk1(2)-Pk(2),Pk1(1)-Pk(1));
% cross-track error
ye=-(x4-Pk(1))*sin(afak)+(x5-Pk(2))*cos(afak);
% expected course angle
psid=afak+atan2(-ye,deta);
% change line
Rk = sqrt((x4-Pk1(1))^2+(x5-Pk1(2))^2);
if pointer ~= length(WayPoint_database(1,:))-1
   if Rk<= R % LOS change radius : 3m
       pointer = pointer+1;       
   end
end
% output
sys(1) = pointer;
sys(2) = psid;
% end mdlOutputs

