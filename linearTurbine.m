close all
clear all
clc

rho = 1.25;
A = 4;
R = 1;

Vop = 6.4;
Wop = 47.1;
Bop = 90;

K = .5*rho*A*R;
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
delta = K31+K32+K33;

Jt=1;
D = gamma/Jt;
tau = .5;

Wref = 47.1;
Bref = 10;
Wt = 0; %IC
V = 4.7+(8.1-4.7)*rand(100,1); %input, need to test different inputs for project

delV= V-Vop;
delWt = Wref-Wt;


Kp = 0;
Ki = 0;
Kd = .5;


s = tf('s');

C = Kp + Ki/s + Kd*s;

%delB = (C*delWt+Bref)/(tau*s+1)-Bop;
for i =1:length(V)

    Wt = (1/(s-D))*((xi/Jt)*delV+(delta/Jt)*(C*delWt+Bref)/(tau*s+1)-Bop)+Wop;
    delWt = Wt-Wref;

end


Wt.InputDelay = 5;
