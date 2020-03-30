close all
clear
clc

t_end = 300;  % simulation end time
tRange = [0 t_end];
Y0 = [0.01; 0.01; 0.01; 0.01];  % initial densities for RS, RL, J, A

[tSol, YSol] = ode45(@model_II, tRange, Y0);

RS = YSol(:,1);
RL = YSol(:,2);
J  = YSol(:,3);
A  = YSol(:,4);

figure
subplot(2,1,1)
plot(tSol, RS, '-g', 'LineWidth', 1)
xlim([0 tRange(end)])
hold on
plot(tSol, RL, '-g', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('RS','RL')
hold off

subplot(2,1,2)
plot(tSol, J, '-k', 'LineWidth', 1)
xlim([0 tRange(end)])
ylim([0 inf])
hold on
plot(tSol, A, '-k', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('J','A')
hold off