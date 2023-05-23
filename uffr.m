function [sys,x0,str,ts,simStateCompliance] = uffr(t,x,u,flag)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=[];

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=[];

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=[];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
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

sizes.NumContStates  = 7;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [2, 0.1, 0, 0.1, 0.1, 0.1, 0.1];

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
% 消防车物理参数
m = 200;
c = 0.05;
r = 0.1;
Jw = 0.005;
Jv = 10;
k = 6;
Fp = 50;
s = 1;
beta1 = 0.1*sin(t)*180/pi;
beta2 = 0.1*sin(t)*180/pi;
l = 0.3;

% 系统参数
a1 = 2*c/(m*r^2 + 2*Jw);
b1 = k*r/(m*r^2 + 2*Jw);
c1 = r*r/(m*r^2 + 2*Jw);

a2 = 2*c*l*l/(Jv*r^2 + 2*Jw*l^2);
b2 = r*l*k/(Jv*r^2 + 2*Jw*l^2);
c2 = r*r/(Jv*r^2 + 2*Jw*l^2);

g1 = 0.9;
g2 = 0.95;
B1 = b1*g1;
B2 = b2*g2;

Fdv = -Fp*cos(beta1+beta1);
Fdp = -s*Fp*sin(beta2);

% 自适应参数
ka1=1;
kto1=1;
ka2=1;
kto2=1;

% 控制参数
k1 = 10;
k2 = 10;
k3 = 10;
epsilon11 = 0.01;
epsilon31 = 0.01;
r2 = 1000;
% r2 = 100;

% 预定性能参数
% T1=5;
% T2=5;
% T1=3;
% T2=3;
T1=1;
T2=1;
ro10=3;
ro100=0.5;
ro20=0.2;
ro200=0.1;

detaud1=1;
detaba1=1;
detaud2=1;
detaba2=1;

% 期望跟踪信号
x1d = 2*sin(t);
x2d = 0.2*sin(t);
dx1d = 2*cos(t);
dx2d = 0.2*cos(t);
ddx1d = -2*sin(t);
ddx2d = -0.2*sin(t);

% 误差信号
z1 = x(1) - x1d;
z2 = x(2) - x2d;
dz2 = x(3) - dx2d;

% 期望性能型号
if t < T1
    ro1 = (ro10-ro100)*((T1-t)/T1)^3 + ro100;
else
    ro1 = ro100;
end

if t < T1  % 一阶导数
    dro1 = (-1/T1)*3*(ro10-ro100)*((T1-t)/T1)^2;
else
    dro1 = 0;
end

if t < T2
    ro2 = (ro20-ro200)*((T2-t)/T2)^3 + ro200;
else
    ro2 = ro200;
end

if t < T2  % 一阶导数
    dro2 = (-1/T2)*3*(ro20-ro200)*((T2-t)/T2)^2;
else
    dro2 = 0;
end

if t < T2  % 二阶导数
    ddro2 = (-1/T2)*(-1/T2)*6*(ro20-ro200)*((T2-t)/T2);
else
    ddro2 = 0;
end

% 转变的误差信号
s1 = z1/ro1;
s2 = z2/ro2;

z1hat = 0.5*log((s1+detaud1)/(detaba1-s1)) - 0.5*log(detaud1/detaba1);
z2hat = 0.5*log((s2+detaud2)/(detaba2-s2)) - 0.5*log(detaud2/detaba2);

R1 = (0.5/ro1)*((1/(s1+detaud1)) - (1/(s1-detaba1)));
R2 = (0.5/ro2)*((1/(s2+detaud2)) - (1/(s2-detaba2)));

dz2hat = R2*(x(3)- dx2d - dro2*z2/ro2);

ds2 = (dz2*ro2-z2*dro2)/(ro2*ro2);
dR2 = -0.25*(dro2/(ro2^2))*(1/(s2+detaud2) + 1/(detaba2-s2)) + (0.5/ro2)*((-1)/((s2+detaud2)^2)  + 1/((detaba2-s2)^2))*ds2;

% 控制信号
Uvbar = (k1/R1)*z1hat - x(4)*x(1) - dx1d - (z1*dro1)/ro1 + 0.5*r2*z1hat*R1;
Uv = -(z1hat*R1*x(5)^2*Uvbar^2)/((z1hat^2*R1^2*x(5)^2*Uvbar^2 + epsilon11^2)^0.5);

alfa2 = -(k2/R2)*z2hat + dx2d + z2*dro2/ro2;
dalfa2 = k2*dR2*z2hat/(R2^2) - (k2/R2)*dz2hat + ddx2d + (1/(ro2^2))*(ro2*(dz2*dro2 + z2*ddro2) - z2*dro2*dro2);

z3 = x(3) - alfa2;

Upbar = R2*z2hat - dalfa2 + k3*z3 - x(6)*x(2) + 0.5*r2*z3;
Up = -(z3*x(7)^2*Upbar^2)/((z3^2*x(7)^2*Upbar^2 + epsilon31^2)^0.5);

% 系统
sys(1) = -a1*x(1) + B1*Uv - c1*Fdv;
sys(2) = x(3);
sys(3) = -a2*x(2) + B2*Up - c2*Fdp;
%自适应a1,to1
sys(4) = -ka1*x(4) - z1hat*R1*x(1);  %a1hat
sys(5) = -kto1*x(5) + z1hat*R1*Uvbar;  %to1hat
%自适应a2,to2
sys(6) = -ka2*x(6) - z3*x(2);  %a2hat
sys(7) = -kto2*x(7) + z3*Upbar;  %to2hat

% end mdlDerivatives

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

% 消防车物理参数
m = 200;
c = 0.05;
r = 0.1;
Jw = 0.005;
Jv = 10;
k = 6;
Fp = 50;
s = 1;
beta1 = 0.1*sin(t)*180/pi;
beta2 = 0.1*sin(t)*180/pi;
l = 0.3;

% 系统参数
a1 = 2*c/(m*r^2 + 2*Jw);
b1 = k*r/(m*r^2 + 2*Jw);
c1 = r*r/(m*r^2 + 2*Jw);

a2 = 2*c*l*l/(Jv*r^2 + 2*Jw*l^2);
b2 = r*l*k/(Jv*r^2 + 2*Jw*l^2);
c2 = r*r/(Jv*r^2 + 2*Jw*l^2);

g1 = 0.9;
g2 = 0.95;
B1 = b1*g1;
B2 = b2*g2;

Fdv = -Fp*cos(beta1+beta1);
Fdp = -s*Fp*sin(beta2);

% 自适应参数
ka1=1;
kto1=1;
ka2=1;
kto2=1;

% 控制参数
k1 = 10;
k2 = 10;
k3 = 10;
epsilon11 = 0.01;
epsilon31 = 0.01;
r2 = 1000;
% r2 = 100;

% 预定性能参数
% T1=5;
% T2=5;
% T1=3;
% T2=3;
T1=1;
T2=1;
ro10=3;
ro100=0.5;
ro20=0.2;
ro200=0.1;

detaud1=1;
detaba1=1;
detaud2=1;
detaba2=1;

% 期望跟踪信号
x1d = 2*sin(t);
x2d = 0.2*sin(t);
dx1d = 2*cos(t);
dx2d = 0.2*cos(t);
ddx1d = -2*sin(t);
ddx2d = -0.2*sin(t);

% 误差信号
z1 = x(1) - x1d;
z2 = x(2) - x2d;
dz2 = x(3) - dx2d;

% 期望性能型号
if t < T1
    ro1 = (ro10-ro100)*((T1-t)/T1)^3 + ro100;
else
    ro1 = ro100;
end

if t < T1  % 一阶导数
    dro1 = (-1/T1)*3*(ro10-ro100)*((T1-t)/T1)^2;
else
    dro1 = 0;
end

if t < T2
    ro2 = (ro20-ro200)*((T2-t)/T2)^3 + ro200;
else
    ro2 = ro200;
end

if t < T2  % 一阶导数
    dro2 = (-1/T2)*3*(ro20-ro200)*((T2-t)/T2)^2;
else
    dro2 = 0;
end

if t < T2  % 二阶导数
    ddro2 = (-1/T2)*(-1/T2)*6*(ro20-ro200)*((T2-t)/T2);
else
    ddro2 = 0;
end

% 转变的误差信号
s1 = z1/ro1;
s2 = z2/ro2;

z1hat = 0.5*log((s1+detaud1)/(detaba1-s1)) - 0.5*log(detaud1/detaba1);
z2hat = 0.5*log((s2+detaud2)/(detaba2-s2)) - 0.5*log(detaud2/detaba2);

R1 = (0.5/ro1)*((1/(s1+detaud1)) - (1/(s1-detaba1)));
R2 = (0.5/ro2)*((1/(s2+detaud2)) - (1/(s2-detaba2)));

dz2hat = R2*(x(3)- dx2d - dro2*z2/ro2);

ds2 = (dz2*ro2-z2*dro2)/(ro2*ro2);
dR2 = -0.25*(dro2/(ro2^2))*(1/(s2+detaud2) + 1/(detaba2-s2)) + (0.5/ro2)*((-1)/((s2+detaud2)^2)  + 1/((detaba2-s2)^2))*ds2;

% 控制信号
Uvbar = (k1/R1)*z1hat - x(4)*x(1) - dx1d - (z1*dro1)/ro1 + 0.5*r2*z1hat*R1;
Uv = -(z1hat*R1*x(5)^2*Uvbar^2)/((z1hat^2*R1^2*x(5)^2*Uvbar^2 + epsilon11^2)^0.5);

alfa2 = -(k2/R2)*z2hat + dx2d + z2*dro2/ro2;
dalfa2 = k2*dR2*z2hat/(R2^2) - (k2/R2)*dz2hat + ddx2d + (1/(ro2^2))*(ro2*(dz2*dro2 + z2*ddro2) - z2*dro2*dro2);

z3 = x(3) - alfa2;

Upbar = R2*z2hat - dalfa2 + k3*z3 - x(6)*x(2) + 0.5*r2*z3;
Up = -(z3*x(7)^2*Upbar^2)/((z3^2*x(7)^2*Upbar^2 + epsilon31^2)^0.5);


% 输出状态
sys(1) = x(1);
sys(2) = x(2);
sys(3) = x(3);
% 输出误差
sys(4) = z1;
sys(5) = z2;
sys(6) = z3;
% 输出性能函数
sys(7) = -detaud1*ro1;
sys(8) = detaba1*ro1;
sys(9) = -detaud2*ro2;
sys(10) = detaba2*ro2;
% 控制信号
sys(11) = Uv;
sys(12) = Up;

% end mdlOutputs
