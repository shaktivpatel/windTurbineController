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

Ka = .5*rho*Ar*R;
lambdaop = (R*Wop)/Vop;

%linearization of wind turbine model
K11 = ((Ka*Vop^3)/(R*Wop))*(.44-.0167*Bop)*((pi*R)/(Vop*(15-.3*Bop)))*cos(pi*((lambdaop-3)/(15-.3*Bop)));
K12 = -((Ka*Vop^3)/(R*Wop))*(.44-.0167*Bop)*sin(pi*((lambdaop-3)/(15-.3*Bop)));
K13 = -.00184*Ka*(Bop*Vop^2+((3*Bop*Vop^3)/(R*Wop^2)));

K21 = (.44-.0167*Bop)*((3*Ka*Vop^2)/(R*Wop))*sin(pi*((lambdaop-3)/(15-.3*Bop)));
K22 = -(.44-.0167*Bop)*((Ka*Vop^3)/(R*Wop))*((pi*lambdaop)/(Vop^2*(15-.3*Bop)))*cos(pi*((lambdaop-3)/(15-.3*Bop)));
K23 = -.00184*Ka*(2*Vop*Bop-((9*Bop*Vop)/lambdaop));

K31 = ((-.0167*Ka*Vop^2)/lambdaop)*sin(pi*((lambdaop-3)/(15-.3*Bop)));
K32 = ((.0167*Ka*Vop^2)/lambdaop)*(.44-.0167*Bop)*(.3*pi*((lambdaop-3)/(15-.3*Bop)^2))*cos(pi*((lambdaop-3)/(15-.3*Bop)));
K33 = (-.00184*Ka*(lambdaop-3)*Vop^2)/lambdaop;

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

%% state feedback controller

% open loop state-space model
[A,B,C,D] = tf2ss(cell2mat(dWdV.Numerator),cell2mat(dWdV.Denominator));

% pole placement
p1 = -2000;
p2 = -3000;

% observor pole placement
op1 = -30;
op2 = -7000;

%calculate gain and matrices
K = place(A,B,[p1 p2]);
L = place(A',C',[op1 op2])';

dWdV_sfc = ss(A-B*K,B,C,D);

At = [A-B*K B*K; zeros(size(A)) A-L*C];
Bt = [B; zeros(size(B))];
Ct = [C zeros(size(C))];

dWdV_osfc = ss(At,Bt,Ct,D);
%% System Response

%% setting inputs
t_in = 0:1.5E-9:.1;

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

[dW,t] = lsim(dWdV_sfc,delV,t_in);
[dW1,t1] = lsim(dWdV_osfc,delV,t_in);
[dW2,t2] = lsim(dWdV,delV,t_in);

subplot(2,1,2);
plot(t,dW+Wop)
hold on
plot(t1,dW1+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Wind Speed Response')
legend('state-feedback','observer state-feedback','plant')

%plot sin input
figure;
subplot(2,1,1);
plot(t_in,Vsin1);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f1) 'Hz'])

[dW,t] = lsim(dWdV_sfc,delVsin1,t_in);
[dW1,t1] = lsim(dWdV_osfc,delVsin1,t_in);
[dW2,t2] = lsim(dWdV,delVsin1,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
hold on
plot(t1,dW1+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response')
legend('state-feedback','observer state-feedback','plant')

figure;
subplot(2,1,1);
plot(t_in,Vsin2);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f2) 'Hz'])

[dW,t] = lsim(dWdV_sfc,delVsin2,t_in);
[dW1,t1] = lsim(dWdV_osfc,delVsin2,t_in);
[dW2,t2] = lsim(dWdV,delVsin2,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
hold on
plot(t1,dW1+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response')
legend('state-feedback','observer state-feedback','plant')

figure;
subplot(2,1,1);
plot(t_in,Vsin3);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f3) 'Hz'])

[dW,t] = lsim(dWdV_sfc,delVsin3,t_in);
[dW1,t1] = lsim(dWdV_osfc,delVsin3,t_in);
[dW2,t2] = lsim(dWdV,delVsin3,t_in);

subplot(2,1,2);
plot(t_in,dW+Wop)
hold on
plot(t1,dW1+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response')
legend('state-feedback','observer state-feedback','plant')

figure;
subplot(2,1,1);
plot(t_in,Vsin4);
xlabel('time');
ylabel('Amplitude')
title(['Input: Sinusoidal Signal with Frequency ' num2str(f4) 'Hz'])

[dW,t] = lsim(dWdV_sfc,delVsin4,t_in);
[dW1,t1] = lsim(dWdV_osfc,delVsin4,t_in);
[dW2,t2] = lsim(dWdV,delVsin4,t_in);

subplot(2,1,2);
plot(t,dW+Wop)
hold on
plot(t1,dW1+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Sin Response')
legend('state-feedback','observer state-feedback','plant')

%plot step input
[dW,t] = step(dWdV_sfc);
[dW1,t1] = step(dWdV_osfc);
[dW2,t2] = step(dWdV);

figure;
subplot(2,1,1);
plot(heaviside(t))
xlabel('time');
ylabel('Amplitude')
title('Input: Step Input')

subplot(2,1,2);
plot(t,dW+Wop)
hold on
plot(t,dW+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Step Response')
legend('state-feedback','observer state-feedback','plant')

%plot impulse input
[dW,t] = impulse(dWdV_sfc);
[dW1,t1] = impulse(dWdV_osfc);
[dW2,t2] = impulse(dWdV);

figure;
subplot(2,1,1);
plot(dirac(t))
xlabel('time');
ylabel('Amplitude')
title('Input: Impulse Input')

subplot(2,1,2);
plot(t,dW+Wop)
hold on
plot(t,dW+Wop)
hold on
plot(t2,dW2+Wop)
xlabel('time');
ylabel('Rotor speed (Wr) [rad/s]');
title('Impulse Response')
legend('state-feedback','observer state-feedback','plant')
