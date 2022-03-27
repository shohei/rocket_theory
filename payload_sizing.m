clear; close all; clc;

C = 4413; %m/s
Isp = 450; %s
deltaV = 9000; %m/s

epsilon = linspace(0.88,0.92,20);
m0_mu = epsilon./(exp(-deltaV/C)+epsilon-1);
plot(epsilon,m0_mu);
ylim([0,Inf]);
title(['単段ロケットのサイジング']);
xlabel('推進薬搭載率 ε');
ylabel('全備質量/ペイロード質量=M0/Mu');
% ylabel('ペイロード質量/全備質量=Mu/M0');

big;