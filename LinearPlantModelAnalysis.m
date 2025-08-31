close all
clear all
clc

%% calculating linearized model of plant (plant = wind turbine)

rho = 1.25;
Ar = 4;
R = 1;

Vop = 6.4;
Wop = 47.1;
Bop = 90;

K = .5*rho*Ar*R;
lambdaop = (R*Wop)/Vop;

%linearization of wind turbine model
K11 = ((K*Vop^3)/(R*Wop))*(.44-.0167*Bop)*((pi*R)/(Vop*(15-.3*Bop)))*cos(pi*((lambdaop-3)/(15-.3*Bop)));
K12 = -((K*Vop^3)/(R*Wop))*(.44-.0167*Bop)*sin(pi*((lambdaop-3)/(15-.3*Bop)));
K13 = -.00184*K*(Bop*Vop^2+((3*Bop*Vop^3)/(R*Wop^2)));

K21 = (.44-.0167*Bop)*((3*K*Vop^2)/(R*Wop))*sin(pi*((lambdaop-3)/(15-.3*Bop)));
K22 = -(.44-.0167*Bop)*((K*Vop^3)/(R*Wop))*((pi*lambdaop)/(Vop^2*(15-.3*Bop)))*cos(pi*((lambdaop-3)/(15-.3*Bop)));
K23 = -.00184*K*(2*Vop*Bop-((9*Bop*Vop)/lambdaop));

K31 = ((-.0167*K*Vop^2)/lambdaop)*sin(pi*((lambdaop-3)/(15-.3*Bop)));
K32 = ((.0167*K*Vop^2)/lambdaop)*(.44-.0167*Bop)*(.3*pi*((lambdaop-3)/(15-.3*Bop)^2))*cos(pi*((lambdaop-3)/(15-.3*Bop)));
K33 = (-.00184*K*(lambdaop-3)*Vop^2)/lambdaop;

gamma = K11+K12+K13;
xi = K21+K22+K23;
delt = K31+K32+K33;

Jt=1;
D = gamma/Jt;
tau = .5;

Wref = 47.1;
Bref = 10;
Wt = 0;

delWt = Wref-Wt;

s = tf('s');

delB = (delWt+Bref)-Bop;

H = (delt/Jt)*delB;
Geq = (1/(s-D))/(1+(1/(s-D))*H);

dWdV = (xi/Jt)*Geq;

%% 4b. analyzing plant

% open loop state space model
[A,B,C,D] = tf2ss(cell2mat(dWdV.Numerator),cell2mat(dWdV.Denominator));

% eigenvalues
e = eig(A);
%e = -30.2996 & -27.1176 --> poles in LHS of plane RE(eig)<0 --> stable

% lyapunov stability
Q = eye(2); 
P = lyap(A',Q);
ep = eig(P);

%% 4c. 

% similarity transform
[T,d] = eig(A);
A_diag = inv(T)*A*T;
B_diag = inv(T)*B;
C_diag = C*T;

% TF of system
[num, denom] = ss2tf(A_diag,B_diag,C_diag,D);
G = tf(num,denom);

% discretized SS models
dt = .1;
Gs = ss(A_diag,B_diag,C_diag,D);
Gz = c2d(Gs,dt,'zoh');
[num, denom] = ss2tf(Gz.A,Gz.B,Gz.C,Gz.D);
TFz = tf(num,denom,dt);

%% System Response

%% setting inputs
t_in = 0:.1:10;

f1 = 10;
f2 = 50;
f3 = 100;
f4 = 500;
V = 4.7+(8.1-4.7)*rand(length(t_in),1); %realistic wind speed
Vsin1 = sin(f1*t_in);
Vsin2 = sin(f2*t_in);
Vsin3 = sin(f3*t_in);
Vsin4 = sin(f4*t_in);

delV= V-Vop;
delVsin1 = Vsin1-Vop;
delVsin2 = Vsin2-Vop;
delVsin3 = Vsin3-Vop;
delVsin4 = Vsin4-Vop;

%% plot system response 

%plot wind speed input
figure;
subplot(2,1,1);
plot(t_in,V);
xlabel('time');
ylabel('Wind Speed [m/s]')
title('Input: Wind Speed')

[dW,t] = lsim(dWdV,delV,t_in);

subplot(2,1,2);
plot(t,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Wind Speed Response of Plant')

%plot sin input
figure;
subplot(2,1,1);
plot(t_in,Vsin1);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f1) 'Hz'])

[dW,t] = lsim(dWdV,delVsin1,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response of Plant')

figure;
subplot(2,1,1);
plot(t_in,Vsin2);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f2) 'Hz'])

[dW,t] = lsim(dWdV,delVsin2,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response of Plant')

figure;
subplot(2,1,1);
plot(t_in,Vsin3);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f3) 'Hz'])

[dW,t] = lsim(dWdV,delVsin3,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response of Plant')

figure;
subplot(2,1,1);
plot(t_in,Vsin4);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f4) 'Hz'])

[dW,t] = lsim(dWdV,delVsin4,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response of Plant')

%plot step input
[dW,t] = step(dWdV);

figure;
subplot(2,1,1);
plot(heaviside(t))
xlabel('time');
ylabel('Amplitude')
title('Input: Step Input')

subplot(2,1,2);
plot(t,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Step Response of Plant')

%plot impulse input
[dW,t] = impulse(dWdV);

figure;
subplot(2,1,1);
plot(dirac(t))
xlabel('time');
ylabel('Amplitude')
title('Input: Impulse Input')

subplot(2,1,2);
plot(t,dW+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Impulse Response of Plant')

