clear; close all; clc;

deltaV_total = 9000; %m/s
C1 = 4413; %m/s
C2 = 4413; %m/s
e1s = [0.75 0.8 0.85 0.9];
e2s = [0.75 0.8 0.85 0.9];

figure();
hold on;
for idx=1:length(e1s)
    e1 = e1s(idx);
    xs = [];
    ys = [];
    for jdx=1:length(e2s)
        e2 = e2s(jdx);
        L1 = sqrt(exp(deltaV_total/C1)*(1-e2)/(1-e1));
        L2 = sqrt(exp(deltaV_total/C1)*(1-e1)/(1-e2));
        m2_mu = (L2-1)/(1-(1-e2)*L2);
        m1_mu = (L1-1)*L2*e2/((1-(1-e1)*L1)*(1-(1-e2)*L2));
        m2_m1 = m2_mu/m1_mu;
        m0_mu = L1*L2*e1*e2/((1-(1-e1)*L1)*(1-(1-e2)*L2));
        xs(end+1) = m0_mu;
        ys(end+1) = m2_m1;
    end
    plot(xs,ys,'r','DisplayName','1段推進薬搭載率一定');
end

for idx=1:length(e1s)
    e2 = e2s(idx);
    xs = [];
    ys = [];
    for jdx=1:length(e2s)
        e1 = e1s(jdx);
        L1 = sqrt(exp(deltaV_total/C1)*(1-e2)/(1-e1));
        L2 = sqrt(exp(deltaV_total/C1)*(1-e1)/(1-e2));
        m2_mu = (L2-1)/(1-(1-e2)*L2);
        m1_mu = (L1-1)*L2*e2/((1-(1-e1)*L1)*(1-(1-e2)*L2));
        m2_m1 = m2_mu/m1_mu;
        m0_mu = L1*L2*e1*e2/((1-(1-e1)*L1)*(1-(1-e2)*L2));
        xs(end+1) = m0_mu;
        ys(end+1) = m2_m1;
    end
    plot(xs,ys,'b','DisplayName','2段推進薬搭載率一定');
end
legend();
% legend('1段推進薬搭載率一定','2段推進薬搭載率一定');

big;



