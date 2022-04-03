clear; close all; clc;

% Input parameters
D = 18e-2; %m
L = 60e-2; %m
R = D/2;
x = linspace(0,L,20);

% Tangent ogive
rho = (R^2+L^2)/(2*R);
y = sqrt(rho^2-(L-x).^2) + R - rho;
subplot(331);
plot(x,y,'k',x,-y,'k');
title('Tangent ogive');

% Conic
y = x*R/L;
subplot(332);
plot(x,y,'k',x,-y,'k');
title('Conic');

% Bi-conic
% mixing_rate1 = 0.35;
% mixing_rate2 = 0.8;
% L1 = mixing_rate1*L;
% L2 = L-L1;
% R2 = R;
% R1 = mixing_rate2*R2;
% y = zeros(size(x));
% for idx=1:length(x)
%   if x(idx) <=L1
%       y(idx) = x(idx)*R1/L1;
%   else
%       y(idx) = R1+(x(idx)-L1)*(R2-R1)/L2;
%   end
% end
% subplot(332);
% plot(x,y,'k',x,-y,'k');
% title('Bi-conic');


% Eliptical
% y = R*sqrt(1-x.^2/L^2);
% subplot(334);
% plot(-x,y,'k',-x,-y,'k');
% title('Eliptical');

% Haack series
theta = acos(1-2*x/L);
C = 1/3;
y = R/sqrt(pi)*sqrt(theta-sin(2*theta)/2+C*sin(theta).^3);
subplot(333);
plot(x,y,'k',x,-y,'k');
title('LV Haack');

C = 0;
y = R/sqrt(pi)*sqrt(theta-sin(2*theta)/2+C*sin(theta).^3);
subplot(334);
plot(x,y,'k',x,-y,'k');
title('Von Kármán');

% Parabolic
Kprime = 0;
y = R*(2*(x/L)-Kprime*(x/L).^2)/(2-Kprime);
subplot(335);
plot(x,y,'k',x,-y,'k');
title('Parabola');

Kprime = 3/4;
y = R*(2*(x/L)-Kprime*(x/L).^2)/(2-Kprime);
subplot(336);
plot(x,y,'k',x,-y,'k');
title('3/4 Parabolic');

Kprime = 1/2;
y = R*(2*(x/L)-Kprime*(x/L).^2)/(2-Kprime);
subplot(337);
plot(x,y,'k',x,-y,'k');
title('1/2 Parabolic');

% Power series
n = 3/4;
y = R*(x/L).^n;
subplot(338);
plot(x,y,'k',x,-y,'k');
title('x^{3/4} Power');

n = 1/2;
y = R*(x/L).^n;
subplot(339);
plot(x,y,'k',x,-y,'k');
title('x^{1/2} Power');



big;
