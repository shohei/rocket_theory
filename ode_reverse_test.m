clear; close all; clc;

x0 = 0;
v0 = 10; %m/s
t_final = 2; %s
tspan = linspace(0,t_final,20);
% 順方向に解く
[t,xnext] = ode45(@myode,tspan,[x0,v0]);
x = xnext(:,1);
v = xnext(:,2);
figure();
plot(t,x);
title('ODEを時間順方向に解いた場合');
big;

x_final = x(end);
v_final = v(end);
t_final = t(end);

% 逆方向に解く。時間を逆にセットし、末端での境界条件を与える。
[t,xnext] = ode45(@myode,flip(tspan),[x_final,v_final]);
x = xnext(:,1);
v = xnext(:,2);
figure();
plot(t,x);
title('ODEを時間逆方向に解いた場合');
big;

function dx = myode(t,xprev)
 g = 9.8;
 x = xprev(1);
 v = xprev(2);
 dx = [v;
      -g];
end
