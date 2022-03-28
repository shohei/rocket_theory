clear; close all; clc;

global C; global tau_m; global mu; 

rf = 6491.14e3; %m
vf = 7832; %m/s
vrf = 0; %m/s

vx0 = 465.1; %m/s 水平方向初速. 赤道上、真東への打ち上げ
C = 4413; %m/s
Isp = 450; %s 真空中比推力
tau_m = 320; %s tau_m = m0/beta. betaは推進薬質量流量. m0はロケット初期質量
% aT = T/m 真空中推力加速度
aT_0 = 13.8; %m2/s. 1.4G
mu = 3.986e14; %m3/s2

%l10,l20,l30,l40,tfを求める




return;
tf = ; %後で求める. シミュレーション時間.
tspan = linspace(0,tf,10000);
l10 = ; %todo
l20= ; %todo
l30 = ; %todo
l40 = ; %todo
[t,x] = ode23s(@adjoint_ode, tspan, [0,0,vx0,0,l10,l20,l30,l40]);


function dx = adjoint_ode(t,xprev)
global C; global tau_m; global mu; 
x = xprev(1);
y = xprev(2);
vx = xprev(3);
vy = xprev(4);
l1 = xprev(5);
l2 = xprev(6);
l3 = xprev(7);
l4 = xprev(8);
r = sqrt(x^2+y^2);
aT = C/(tau_m-t);
cos_th = l3/sqrt(l3^2+l4^2);
sin_th = l4/sqrt(l3^2+l4^2);
dx = [vx;
    vy;
    aT*cos_th-mu/r^3*x;
    aT*sin_th-mu/r^3*y;
    mu/r^3*(1-3*x^2/r^2)*l3-3*mu/r^5*x*y*l4;
    mu/r^3*(1-3*y^2/r^2)*l4-3*mu/r^5*x*y*l3;
    -l1;
    -l2];

end






