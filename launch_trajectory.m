clear; close all; clc;

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






function dlambda = adjoint_ode(t,x)
l1 = x(1);
l2 = x(2);
l3 = x(3);
l4 = x(4);
dlambda = [mu/r^3*(1-3*x^2/r^2)*l3-3*mu/r^5*x*y*l4;
           mu/r^3*(1-3*y^2/r^2)*l4-3*mu/r^5*x*y*l3;
           -l1;
           -l2];
end






