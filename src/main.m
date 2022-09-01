% Main function
%
% Example on the application of Lie group integrators for the solution of 
% the equations of motion of a multibody system, August 2022, Andrea Leone @ NTNU
%
% Ref.: E. Celledoni, E. Çokaj, A. Leone, D. Murari, B. Owren.
% "Lie group integrators for mechanical systems",
% International Journal of Computer Mathematics, 99:1, 58-88.

clear all;
close all;
clc;

% add paths
addpath('Lie_group_functions'); addpath('integrators');
addpath('pre_post_processing'); addpath('control_functions');

%% ------------------------------------------------------------------------
% Parameters of the mechanical system, initial conditions, desired trajectory

T = 5; % integration time
N = 2000; % number of time steps

J1 = diag([0.0820, 0.0845, 0.1377]);
J2 = J1;
my = 0.4;
m1 = 0.755;
m2 = 0.755;
L1 = 0.6;
L2 = 0.8;
g = 9.81;
e3 = [0;0;1];

% Initial conditions
R1 = eye(3);
R2 = R1;
y0 = [-2;0;0];
q10 = [sin(deg2rad(80)); 0; cos(deg2rad(80))];
q20 = [-sin(deg2rad(80)); 0; cos(deg2rad(80))];
v = zeros(3,1);
w1=v; 
w2 = w1;
Omega1 = v; 
Omega2 = v;

z0 = [y0, v, R1, Omega1, R2, Omega2, q10, w1, q20, w2];

% Desidered trajectory (e.g. parabolic trajectory) with derivatives 
% needed to define controls
yd = @(t) [t-2.5;0;(t-2.5).^2-7];
yd_dot = @(t) [1;0;2*(t-2.5)];
yd_ddot = @(t) [0;0;2];
yd_dddot = @(t) [0;0;0];
yd_ddddot = @(t) [0;0;0];
yd_dddddot = @(t) [0;0;0];
yd_ddddddot = @(t) [0;0;0];

d_trajectory = @(t) [yd(t),yd_dot(t),yd_ddot(t),yd_dddot(t),yd_ddddot(t),yd_dddddot(t),yd_ddddddot(t)];

%% ------------------------------------------------------------------------
% group action, exponential, dexpinv

% Representation of a generic element of the phase space TR3x(TSO(3)^2)x(TS2^2):
% P = [y, v, R1, Omega1, R2, Omega2, q1, w1, q2, w2] in R^{3x14} 
% v, y, v Omega1, Omega2, q1, w1, q2, w2 are 3x1 vecotrs
% R1 and R2 are 3x3 rotation matrices

% Representation of a generic element of the group G = R6x(TSO(3)^2)x(SE(3))^2: 
% B = [a1,a2,   B1,b1,B2,b2,   C1,c1,C2,c2] in R^{3x18} 
% a1, a2 3x1 vectors; B1, B2, C1, C2 3x3 matrices; b1, b2, c1, c2 3x1 vectors

% Representation of a generic element of the Lie algebra g: 
% v = [xi1;xi2;eta1;eta2;eta3;eta4;mu1;mu2;mu3;mu4] in R^30 
% all the entries are 3x1 vectors

action = @(B,P) [P(:,1) + B(:,1), ...
                 B(:,2) + P(:,2), ... 
                 B(:,3:5)*P(:,3:5), ...
                 B(:,6) + P(:,6), ...
                 B(:,7:9)*P(:,7:9), ...
                 B(:,10) + P(:,10), ... 
                 B(:,11:13)*P(:,11), ...
                 B(:,11:13)*P(:,12) + hat(B(:,14))*B(:,11:13)*P(:,11), ...
                 B(:,15:17)*P(:,13), ...
                 B(:,15:17)*P(:,14) + hat(B(:,18))*B(:,15:17)*P(:,13)]; 
 
exponential = @(v) [v(1:3), v(4:6), expSO3(v(7:9)), v(10:12),...
                    expSO3(v(13:15)), v(16:18), expSE3(v(19:24)), expSE3(v(25:30))];
      
dexpinv = @(v,u) [u(1:3); u(4:6); dexpinvSO3(v(7:9),u(7:9)) ; u(10:12); ...
                  dexpinvSO3(v(13:15),u(13:15)) ; u(16:18); ...
                  dexpinvSE3(v(19:24),u(19:24)); dexpinvSE3(v(25:30),u(25:30)) ];

%% ------------------------------------------------------------------------
% Vector field on the Lie algebra 

Mq = @(P) my*eye(3) + m1*P(:,11)*P(:,11)' + m2*P(:,13)*P(:,13)';

A = @(P) [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
          zeros(3), Mq(P), zeros(3), zeros(3), zeros(3), zeros(3) ;
          zeros(3), zeros(3), J1, zeros(3), zeros(3), zeros(3);
          zeros(3), zeros(3), zeros(3), J2,  zeros(3), zeros(3);
          zeros(3), -1/L1*hat(P(:,11)), zeros(3),  zeros(3), eye(3), zeros(3);
          zeros(3), -1/L2*hat(P(:,13)), zeros(3),  zeros(3), zeros(3), eye(3)];

% Support functions to define the vector field 

getvi = @(v,i) v(3*i-2:3*i); 
GetPerp = @(qi,u) -hat(qi)^2*u;
GetPar  = @(qi,u) (qi'*u)*qi;
component = @(v,i) v(3*i-2:3*i);
flat = @(M) M(:);

upar = @(P,trajectory,i) component(getUpar(m1,m2,my,L1,L2,P,trajectory),i);
uperp = @(P,trajectory,i) component(getUperp(m1,m2,my,L1,L2,P,trajectory),i);

M = @(P,trajectory,i) ones(3,1);

h1 = @(P,trajectory) P(:,2);
h2 = @(P,trajectory) -m1*L1*norm(P(:,12),2)^2*P(:,11)-m2*L2*norm(P(:,14),2)^2*P(:,13) + Mq(P)*g*e3 + upar(P,trajectory,1) + upar(P,trajectory,2);
h3 = @(P,trajectory) -hat(P(:,6))*(J1*P(:,6)) + M(P,trajectory,1);
h4 = @(P,trajectory) -hat(P(:,10))*(J2*P(:,10)) + M(P,trajectory,2);
h5 = @(P,trajectory) -1/L1 * g * hat(P(:,11))*e3 - 1/(m1*L1) * hat(P(:,11)) * uperp(P,trajectory,1);
h6 = @(P,trajectory) -1/L2 * g * hat(P(:,13))*e3 - 1/(m2*L2) * hat(P(:,13)) * uperp(P,trajectory,2);

h = @(P,trajectory) [h1(P,trajectory); h2(P,trajectory); h3(P,trajectory); h4(P,trajectory); h5(P,trajectory); h6(P,trajectory)];
hbar = @(P,trajectory) A(P)\h(P,trajectory);

f = @(B,trajectory) [getvi(hbar(B,trajectory),1); getvi(hbar(B,trajectory),2); B(:,3:5)*B(:,6); getvi(hbar(B,trajectory),3); B(:,7:9)*B(:,10); getvi(hbar(B,trajectory),4); ...
            B(:,12); hat(B(:,11))*getvi(hbar(B,trajectory),5); B(:,14); hat(B(:,13))*getvi(hbar(B,trajectory),6)];

% vector field
vector_field = @(sigma,B,trajectory) dexpinv(sigma,f(action(exponential(sigma),B),trajectory));

%% ------------------------------------------------------------------------
% Reference solution with ODE45 (solRef)

eqR1 = @(B) flat(B(:,3:5)*hat(B(:,6)));
eqR2 = @(B) flat(B(:,7:9)*hat(B(:,10)));

getBi = @(B,i) B(:,i);

odeFunc = @(t,B) [getvi(hbar(reshape(B,3,14),d_trajectory(t)),1); getvi(hbar(reshape(B,3,14),d_trajectory(t)),2); eqR1(reshape(B,3,14)); getvi(hbar(reshape(B,3,14),d_trajectory(t)),3); eqR2(reshape(B,3,14)); getvi(hbar(reshape(B,3,14),d_trajectory(t)),4); ...
            hat(getBi(reshape(B,3,14),12))*getBi(reshape(B,3,14),11); getvi(hbar(reshape(B,3,14),d_trajectory(t)),5); hat(getBi(reshape(B,3,14),14))*getBi(reshape(B,3,14),13); getvi(hbar(reshape(B,3,14),d_trajectory(t)),6)];
   
ic = z0(:);
options = odeset('AbsTol',1e-14,'RelTol',1e-10);
[timeSol,zC] = ode45(odeFunc, [0 T], ic, options);
zC = zC';

q1ref = zeros(3, length(timeSol));
q2ref = zeros(3, length(timeSol));
R1ref = zeros(3,3, length(timeSol));
R2ref = zeros(3,3, length(timeSol));
detR1ref = zeros(1, length(timeSol));
detR2ref = zeros(1, length(timeSol));

for i = 1:length(timeSol)
    zRef = zC(:,i);
    solRef = reshape(zRef,3,14);
    q1ref(:,i) = solRef(:,11);
    q2ref(:,i) = solRef(:,13);
    R1ref(:,:,i) = solRef(:,3:5);
    R2ref(:,:,i) = solRef(:,7:9);
    detR1ref(i) = det(R1ref(:,:,i));
    detR2ref(i) = det(R2ref(:,:,i));
end

%% ------------------------------------------------------------------------
% Sol with Lie group integrators and plots

time = linspace(0,T,N);
dt = time(2)-time(1);

% Initialization of the solutions with the 5 different integrators
z1 = z0;
z2 = z0;
z3 = z0;
z4 = z0;
z5 = z0;

sigma0 = zeros(30,1);

% Preprocessing for errors, controls and plots
preprocess;

for i = 1:N-1

    t=time(i);
    
%     z1 = LieEuler(vector_field,exponential,action,z1,dt,sigma0,d_trajectory,t);
%     z2 = RKMK2Heun(vector_field,exponential,action,z2,dt,sigma0,d_trajectory,t);
%     z3 = RKMK3(vector_field,exponential,action,z2,dt,sigma0,d_trajectory,t);
    z4 = RKMK4(vector_field,exponential,action,z4,dt,sigma0,d_trajectory,t);
%     z5 = CFree4(f,action,exponential,dt,z5,sigma0,d_trajectory,t);

    sol=z4; % solution at time t=time(i) with integrator RKMK4 in this case

% Postprocessing for errors, controls and plots
    postprocess;

end

% Plots of trajectories and errors 
plots;
