clear; close all; clc;

% 計算条件の定義
global C; global tau_m; global mu;

rf = 6491.14e3; %m
vf = 7832; %m/s
vrf = 0; %m/s

vy0 = 1e-3; %m/s 垂直方向初速. ゼロにするとエラーになるので微小の値を入れた.
vx0 = 465.1; %m/s 水平方向初速. 赤道上、真東への打ち上げ
re = 6378.142e3; %m 地球の半径
x0 = 0;
y0 = re;

C = 4413; %m/s
Isp = 450; %s 真空中比推力
tau_m = 320; %s tau_m = m0/beta. betaは推進薬質量流量. m0はロケット初期質量
% aT = T/m 真空中推力加速度
aT_0 = 13.8; %m2/s. 1.4G
mu = 3.986e14; %m3/s2

tf = 274.1; %[s] エンジン燃焼時間
div_time = 1000; %時間刻み
tspan = linspace(0,tf,div_time);

% シューティング法（Backward sweep）による繰り返し計算
% λの初期値変数の初期値を仮定
% l1_0 = 1; l2_0 = 1; l3_0 = 1; l4_0 = 1;
l1_0 = 6; l2_0 = 6; l3_0 = 1; l4_0 = -1;

% λの初期値変数の定義と初期化
global l1i; global l2i; global l3i; global l4i; 
l1i = l1_0; l2i = l2_0; l3i = l3_0; l4i = l4_0;
% λの修正量δλの定義と初期化
global dl1; global dl2; global dl3; global dl4;
dl1 = 0; dl2 = 0; dl3 = 0; dl4 = 0;

figure(1);
hold on;
title('Error curve');
xlabel('Iteration');
ylabel('Error');

iter = 0;
Es = [];
iters = [];

%% メイン計算部
tic;
while true %無限ループ
% λの初期値を仮定
l1i = l1i + dl1; 
l2i = l2i + dl2; 
l3i = l3i + dl3; 
l4i = l4i + dl4;
[t_,xnext] = ode23s(@adjoint_ode, tspan, [x0,y0,vx0,vy0,l1i,l2i,l3i,l4i]);
x_ = xnext(:,1);
y_ = xnext(:,2);
vx_ = xnext(:,3);
vy_ = xnext(:,4);
l1_ = xnext(:,5);
l2_ = xnext(:,6);
l3_ = xnext(:,7);
l4_ = xnext(:,8);
aT_ = C./(tau_m-t_);
r_ = sqrt(x_.^2+y_.^2);

x_final = x_(end);
y_final = y_(end);
vx_final = vx_(end);
vy_final = vy_(end);
l1_final = l1_(end);
l2_final = l2_(end);
l3_final = l3_(end);
l4_final = l4_(end);
aT_final = aT_(end);
r_final = r_(end);

syms x y vx vy l1 l2 l3 l4 aT;
f1 = vx;
f2 = vy;
f3 = aT*sqrt(vx^2/(vx^2+vy^2))-mu/sqrt(x^2+y^2)^3*x;
f4 = aT*sqrt(vy^2/(vx^2+vy^2))-mu/sqrt(x^2+y^2)^3*y;
H = 1 + l1*f1 + l2*f2 + l3*f3 + l4*f4;
% 状態方程式fのヤコビ行列の計算
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
% Hのヘッセ行列の計算
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
A_sym = df_dx;
B_sym = 0;
C_sym = d2H_dx2;

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

% ターミナルコストΦ(終端条件)のグラディエントの計算
dphi_dx = [diff(phi,x); diff(phi,y); diff(phi,vx); diff(phi,vy)];
dphi_dx_transpose_final = eval(subs((dphi_dx)',[x,y,vx,vy,aT],...
                         [x_final,y_final,vx_final,vy_final,aT_final]));

% エラーの計算
E = [l1_final;l2_final;l3_final;l4_final] -...
     dphi_dx_transpose_final * [x_final;y_final;vx_final;vy_final];
error = max(abs(E));
rf_error = abs(r_final - rf);
vrf_error = abs(...
              eval(subs((r*dtheta_dt),[x,y,vx,vy,aT],...
                 [x_final,y_final,vx_final,vy_final,aT_final])) - vf);
vr_error = abs(eval(subs(dr_dt,[x,y,vx,vy],...
               [x_final,y_final,vx_final,vy_final])) - vrf);

fprintf('Loop: %d, Error: %.3f, Rf_error: %.2f[km], Vrf_error: %.2f[m/s], Vr_error: %.2f[m/s] \n',...
                  iter,error,rf_error/1000,vrf_error/1000,vr_error/1000);
figure(1);
iters(end+1) = iter;
Es(end+1) = error;
plot(iters,Es,'k');
set(gca, 'YScale', 'log');
drawnow();
if rf_error < 10e3  
    disp('Simulation end');
    break;
end
if error < 0.001
  disp('Simulation end');
    break;
end

%初回に一度make_ABCを実行してAt,Bt,Ctを作成しておく.
%  [At,Bt,Ct] = make_ABC(A_sym,B_sym,C_sym,t_,x_,y_,vx_,vy_,aT_,l1_,l2_,l3_,l4_);
%  save('ABC.mat','At','Bt','Ct');

%2回目以降はAt,Bt,Ctを呼び出して使う
load('ABC','At','Bt','Ct');
[t,Sprev] = ode23s(@(t,Snext) riccati_ode1(t,Snext,At,Bt,Ct,t_),flip(tspan),S_final);

c_final = -E;
[t,cprev] = ode23s(@(t,cnext) riccati_ode2(t,cnext,At,Bt,Ct,Sprev,t_),flip(tspan),c_final);
% tは時系列的に逆順に入っているので、cprevの最後の要素がc(0)になる。
c_zero = cprev(end,:);
dl1 = c_zero(1);
dl2 = c_zero(2);
dl3 = c_zero(3);
dl4 = c_zero(4);

iter = iter + 1;
end
toc;

subplot(221);
plot(x_/1000,y_/1000);
title('赤道面上での軌跡');
xlabel('水平位置 x [km]')
ylabel('垂直距離 y [km]');

subplot(222);
% 地表からの高度
h_ = (r_ - re); %m
plot(t_,h_/1000);
title('高度');
xlabel('Time [s]');
ylabel('高度 h [km]');

subplot(223);
% x軸から測った機軸姿勢角θ
theta_ = atan2(vy_,vx_);
plot(t_,theta_/pi*180);
title('推力軸姿勢角θ');
xlabel('Time [s]');
ylabel('θ [deg]');

% 飛行レンジ角Φ
phi_ = atan2(x_,y_);
% 局所水平面から測った機軸姿勢角θL
theta_L_ = theta_ + phi_;
subplot(224);
plot(t_,theta_L_/pi*180);
title('局所水平面に対する機軸姿勢角θL');
xlabel('Time [s]');
ylabel('θL [deg]');

big;
grid;


function dc = riccati_ode2(t,c,At,Bt,Ct,St,t_all)
N = length(t_all);
St = reshape(St, [4,4,N]);
[val,t_idx]=min(abs(t_all-t));
A_t = At(:,:,t_idx);
B_t = Bt(:,:,t_idx);
C_t = Ct(:,:,t_idx);
S_t = St(:,:,t_idx);
dc = (S_t*B_t-A_t')*c;
end


function [At, Bt, Ct] = make_ABC(A_sym,B_sym,C_sym,t_all,x_,y_,vx_,vy_,aT_,l1_,l2_,l3_,l4_)
tic;
syms x y vx vy aT l1 l2 l3 l4;
N = length(t_all);
At = zeros(4,4,N);
Bt = zeros(4,4,N);
Ct = zeros(4,4,N);
for t_idx=1:length(t_all)
    t_idx
    x_t = x_(t_idx); y_t = y_(t_idx); vx_t = vx_(t_idx); 
    vy_t = vy_(t_idx); aT_t = aT_(t_idx); 
    l1_t = l1_(t_idx); l2_t = l2_(t_idx); l3_t = l3_(t_idx); 
    l4_t = l4_(t_idx);
    A = eval(subs(A_sym,[x y vx vy aT],[x_t,y_t,vx_t,vy_t,aT_t]));
    B = B_sym;
    C = eval(subs(C_sym,[x y vx vy l1 l2 l3 l4 aT],...
        [x_t,y_t,vx_t,vy_t,l1_t,l2_t,l3_t,l4_t,aT_t]));
    At(:,:,t_idx) = A;
    Bt(:,:,t_idx) = B;
    Ct(:,:,t_idx) = C;
end
toc;
end

function dS = riccati_ode1(t,S,At,Bt,Ct,t_all)
S = reshape(S,[4,4]);
[val,t_idx]=min(abs(t_all-t));
A_t = At(:,:,t_idx);
B_t = Bt(:,:,t_idx);
C_t = Ct(:,:,t_idx);
dS = -A_t'*S-S*A_t-S*B_t*S+C_t;
dS = dS(:);
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

