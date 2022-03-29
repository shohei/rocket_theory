clear; close all; clc;

C = 4413; %m/s
lambdas = [0.025, 0.05]; %ペイロード比 2.5%, 5%
epsilon = linspace(0.85,0.95,10); %推進薬搭載率 85%~95%
figure();
hold on;
for idx=1:length(lambdas)
    lambda = lambdas(idx);
    deltaV = C*log(1./(1-(1-lambda).*epsilon));
    plot(epsilon,deltaV);
end
legend(arrayfun(@(x) sprintf('lambda=%.3f',x), lambdas,UniformOutput=false));
title('単段式ロケットの推進薬搭載率と理想速度増分');
xlabel('推進薬搭載率ε');
ylabel('速度増分⊿V[m/s]');
big;

figure();
hold on;
a0_g = linspace(1,3.5,50);
Lambdas = [2 2.5 3 3.5 4]; %質量比の逆数

y_pks = [];
a0_g_peaks = [];
for idx=1:length(Lambdas)
    Lambda = Lambdas(idx);
    y_g_Isp2 = 1./a0_g.*(1+...
                -1/Lambda*(1+log(Lambda))...
                -0.5*(1./a0_g)*(1-1./Lambda).^2);   
    [y_pk,loc] = findpeaks(y_g_Isp2);
    a0_g_peak = a0_g(loc);
    y_pks(end+1) = y_pk;
    a0_g_peaks(end+1) = a0_g_peak;
    plot(a0_g, y_g_Isp2);
end
plot(a0_g_peaks, y_pks, 'xk');
legend(arrayfun(@(x) sprintf('Lambda=%.1f',x), Lambdas,UniformOutput=false));
title('垂直に上昇飛行するロケットの無次元到達高度');
xlabel('発射時の加速度倍数G=a0/g0');
ylabel('$y/g I_{sp}^2$','Interpreter','latex');
big;
