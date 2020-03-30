close all
clear
clc

t_end = 300;  % simulation end time
tRange = [0 t_end];
Y0 = [0.01; 0.01; 0.01; 0.01; 0.01; 0.1];  % initial densities for RS, RL, J1, A1, J2, A2

[tSol, YSol] = ode45(@model_III, tRange, Y0);

RS = YSol(:,1);
RL = YSol(:,2);
J1 = YSol(:,3);
A1 = YSol(:,4);
J2 = YSol(:,5);
A2 = YSol(:,6);

figure
subplot(3,1,1)
plot(tSol, RS, '-g', 'LineWidth', 1)
xlim([0 tRange(end)])
hold on
plot(tSol, RL, '-g', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Algal';'density (mg/L)'})
legend('RS','RL')
hold off

subplot(3,1,2)
plot(tSol, J1, '-k', 'LineWidth', 1)
xlim([0 tRange(end)])
ylim([0 inf])
hold on
plot(tSol, A1, '-k', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('J1','A1')
hold off

subplot(3,1,3)
plot(tSol, J2, '-k', 'LineWidth', 1)
xlim([0 tRange(end)])
ylim([0 inf])
hold on
plot(tSol, A2, '-k', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('J2','A2')
hold off