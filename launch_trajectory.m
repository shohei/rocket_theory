clear; close all; clc;
global count;
count=0;

global C; global tau_m; global mu;

rf = 6491.14e3; %m
vf = 7832; %m/s
vrf = 0; %m/s

vy0 = 0; %m/s 垂直方向初速
vx0 = 465.1; %m/s 水平方向初速. 赤道上、真東への打ち上げ
re = 6378.142e3; %m
x0 = 0;
y0 = re;

C = 4413; %m/s
Isp = 450; %s 真空中比推力
tau_m = 320; %s tau_m = m0/beta. betaは推進薬質量流量. m0はロケット初期質量
% aT = T/m 真空中推力加速度
aT_0 = 13.8; %m2/s. 1.4G
mu = 3.986e14; %m3/s2

tf = 274.1; %[s] エンジン燃焼時間
tspan = linspace(0,tf,10000);
l10 = 1; %todo
l20= 1; %todo
l30 = 1; %todo
l40 = 1; %todo
[t_,xnext] = ode23s(@adjoint_ode, tspan, [x0,y0,vx0,vy0,l10,l20,l30,l40]);
x_ = xnext(:,1);
y_ = xnext(:,2);
vx_ = xnext(:,3);
vy_ = xnext(:,4);
l1_ = xnext(:,5);
l2_ = xnext(:,6);
l3_ = xnext(:,7);
l4_ = xnext(:,8);
aT_ = C./(tau_m-t_);
x_final = x_(end);
y_final = y_(end);
vx_final = vx_(end);
vy_final = vy_(end);
l1_final = l1_(end);
l2_final = l2_(end);
l3_final = l3_(end);
l4_final = l4_(end);
aT_final = C/(tau_m-t_(end));
r_final = sqrt(x_final^2+y_final^2);
E1 = r_final-rf;
E2 = sqrt(vx_final^2+vy_final^2) - vf;
E3 = x_final/r_final*vx_final+y_final/r_final*vy_final - vrf;
E4 = l1_final*y_final - l2_final*x_final;
cos_th_final = l3_final/sqrt(l3_final^2+l4_final^2);
sin_th_final = l4_final/sqrt(l3_final^2+l4_final^2);
E5 = 1+l1_final*vx_final+...
    l2_final*vy_final+...n
    l3_final*(aT_final*cos_th_final-mu/r_final^3*x_final)+...
    l4_final*(aT_final*sin_th_final-mu/r_final^3*y_final);
% if max(E1,E2,E3,E4,E5) < 1e-2
%   break;
% end

syms x y vx vy l1 l2 l3 l4 aT;
f1 = vx;
f2 = vy;
f3 = aT*sqrt(vx^2/(vx^2+vy^2))-mu/sqrt(x^2+y^2)^3*x;
f4 = aT*sqrt(vy^2/(vx^2+vy^2))-mu/sqrt(x^2+y^2)^3*y;
H = 1 + l1*f1 + l2*f2 + l3*f3 + l4*f4;
df1_dx = diff(f1,x);
df1_dy = diff(f1,y);
df1_dvx = diff(f1,vx);
df1_dvy = diff(f1,vy);
df2_dx = diff(f2,x);
df2_dy = diff(f2,y);
df2_dvx = diff(f2,vx);
df2_dvy = diff(f2,vy);
df3_dx = diff(f3,x);
df3_dy = diff(f3,y);
df3_dvx = diff(f3,vx);
df3_dvy = diff(f3,vy);
df4_dx = diff(f4,x);
df4_dy = diff(f4,y);
df4_dvx = diff(f4,vx);
df4_dvy = diff(f4,vy);
df_dx = [df1_dx df1_dy df1_dvx df1_dvy;
    df2_dx df2_dy df1_dvx df2_dvy;
    df3_dx df3_dy df3_dvx df3_dvy;
    df4_dx df4_dy df4_dvx df4_dvy];
d2H_dxdx = diff(diff(H,x),x);
d2H_dxdy = diff(diff(H,x),y);
d2H_dxdvx = diff(diff(H,x),vx);
d2H_dxdvy = diff(diff(H,x),vy);
d2H_dydx = diff(diff(H,y),x);
d2H_dydy = diff(diff(H,y),y);
d2H_dydvx = diff(diff(H,y),vx);
d2H_dydvy = diff(diff(H,y),vy);
d2H_dvxdx = diff(diff(H,vx),x);
d2H_dvxdy = diff(diff(H,vx),y);
d2H_dvxdvx = diff(diff(H,vx),vx);
d2H_dvxdvy = diff(diff(H,vx),vy);
d2H_dvydx = diff(diff(H,vy),x);
d2H_dvydy = diff(diff(H,vy),y);
d2H_dvydvx = diff(diff(H,vy),vx);
d2H_dvydvy = diff(diff(H,vy),vy);

d2H_dx2 = [d2H_dxdx d2H_dxdy d2H_dxdvx d2H_dxdvy;
    d2H_dydx d2H_dydy d2H_dydvx d2H_dydvy;
    d2H_dvxdx d2H_dvxdy d2H_dvxdvx d2H_dvxdvy;
    d2H_dvydx d2H_dvydy d2H_dvydvx d2H_dvydvy];
At = df_dx;
Bt = 0;
Ct = d2H_dx2;
    r = sqrt(x^2+y^2);
cos_theta = vx/sqrt(vx^2+vy^2);
sin_theta = vy/sqrt(vx^2+vy^2);
dtheta_dt = cos_theta^2*(-vy/vx^2*(aT*cos_theta-mu/r^3*x)+...
    1/vx*(aT*sin_theta-mu/r^3*y));
dr_dt = x/sqrt(x^2+y^2)*vx + y/sqrt(x^2+y^2)*vy;
phi1 = r - rf;
phi2 = r*dtheta_dt - vf;
phi3 = dr_dt - vrf;
phi = phi1+phi2+phi3; % Terminal cost: scalar function
d2phi_dxdx = diff(diff(phi,x),x);
d2phi_dxdy = diff(diff(phi,x),y);
d2phi_dxdvx = diff(diff(phi,x),vx);
d2phi_dxdvy = diff(diff(phi,x),vy);
d2phi_dydx = diff(diff(phi,y),x);
d2phi_dydy = diff(diff(phi,y),y);
d2phi_dydvx = diff(diff(phi,y),vx);
d2phi_dydvy = diff(diff(phi,y),vy);
d2phi_dvxdx = diff(diff(phi,vx),x);
d2phi_dvxdy = diff(diff(phi,vx),y);
d2phi_dvxdvx = diff(diff(phi,vx),vx);
d2phi_dvxdvy = diff(diff(phi,vx),vy);
d2phi_dvydx = diff(diff(phi,vy),x);
d2phi_dvydy = diff(diff(phi,vy),y);
d2phi_dvydvx = diff(diff(phi,vy),vx);
d2phi_dvydvy = diff(diff(phi,vy),vy);
S = [d2phi_dxdx d2phi_dxdy d2phi_dxdvx d2phi_dxdvy;
    d2phi_dydx d2phi_dydy d2phi_dydvx d2phi_dydvy;
    d2phi_dvxdx d2phi_dvxdy d2phi_dvxdvx d2phi_dvxdvy;
    d2phi_dvydx d2phi_dvydy d2phi_dvydvx d2phi_dvydvy];
S_final = eval(subs(S,[x,y,vx,vy,aT],[x_final,y_final,vx_final,vy_final,aT_final]));
S_final = reshape(S_final, [16,1]);
[t,Snext] = ode23s(@(t,Sprev) riccati_ode1(...
    t,Sprev,At,Bt,Ct,t_,x_,y_,vx_,vy_,l1_,l2_,l3_,l4_,aT_),...
    flip(tspan),S_final);


function dS = riccati_ode1(t,S,At,Bt,Ct,t_,x_,y_,vx_,vy_,l1_,l2_,l3_,l4_,aT_)

t_idx = find(t_==t);
x_t = x_(t_idx); y_t = y_(t_idx); vx_t = vx_(t_idx); vy_t = vy_(t_idx);
l1_t = l1_(t_idx); l2_t = l2_(t_idx); l3_t = l3_(t_idx); l4_t = l4_(t_idx);
aT_t = aT_(t_idx);
syms x y vx vy aT l1 l2 l3 l4;
A_t = eval(subs(At,[x y vx vy aT],[x_t,y_t,vx_t,vy_t,aT_t]));
B_t = Bt;
C_t = eval(subs(Ct,[x y vx vy l1 l2 l3 l4 aT],...
    [x_t,y_t,vx_t,vy_t,l1_t,l2_t,l3_t,l4_t,aT_t]));
S = reshape(S,[4,4]);
dS = -A_t'*S-S*A_t-S*B_t*S+C_t;
dS = reshape(dS,[16,1]);

end


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

