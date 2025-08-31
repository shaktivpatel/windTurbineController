
close all
clear all
clc

tau = .5;
J = .75;
B = 2;

s =tf('s');

%w/o PID control
G_pc = tf([1],[tau 1],'InputDelay',1.5);
G_dt = tf([1],[J B],'InputDelay',1.5);

Geq = ((G_pc*G_dt)/(1+G_dt))/(1+(G_pc*G_dt)/(1+G_dt));

[A,B,C,D] = ssdata(Geq);

figure;
[y,t] =step(Geq,50);
plot(t,y)



%w/ PID control
Kp = 1;
Ki = .9;
Kd = .1;
C_PI = Kp + Ki/s + Kd*s;

Geq_PI = ((C_PI*G_pc*G_dt)/(1+G_dt))/(1+(C_PI*G_pc*G_dt)/(1+G_dt));

figure;
[y,t] =step(Geq_PI,50);
plot(t,y)

