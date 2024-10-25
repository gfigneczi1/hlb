clear; close all;

a = load("C:\git\KDP\publications\JAT_2024\Data\testI.mat");
b = load("C:\git\KDP\publications\JAT_2024\Data\testII.mat");
c = load("C:\git\KDP\publications\JAT_2024\Data\testIII.mat");


f = figure();
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
set(f, "Position", [100,100,1250,850]);
subplot(3,3,1);
% error plot
plot(a.Ts:a.Ts:a.N*a.Ts, a.errorInstinctivePathToPath, 'DisplayName','Planned lane offset', 'LineWidth',1.5, 'Color','k', 'LineStyle','--');
grid on; hold on;
plot(a.Ts:a.Ts:a.N*a.Ts, a.error, 'DisplayName','Real lane offset', 'LineWidth',1.5,'Color','k');
ylim([-1,1]); xlabel('time(s)'); ylabel('$\delta(m)$');
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
xlim([a.Ts,a.N*a.Ts]);
set(gca,'FontSize',14);
legend('Location','best', 'FontSize',11);
title("\textbf{Test I.}");

subplot(3,3,4);
% input plot
plot(a.Ts:a.Ts:a.N*a.Ts, a.u_saved, 'DisplayName', 'Steering angle', 'LineWidth',1.5, 'color', 'k');
grid on;
xlabel('time(s)'); ylabel('$\alpha_f(rad)$');
ylim([-0.05,0.05]); xlim([a.Ts,a.N*a.Ts]);
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
set(gca,'FontSize',14);
%legend('Location','best', 'FontSize',11);

subplot(3,3,7);
% cost plot
plot(a.Ts:a.Ts:a.N*a.Ts, a.fval(:,1), 'DisplayName', 'Total cost', 'LineWidth',2,'color','k');
hold on; grid on;
plot(a.Ts:a.Ts:a.N*a.Ts, a.fval(:,2), 'DisplayName', 'Output error cost', 'LineWidth',1.5,'color','k', 'LineStyle', '--');
plot(a.Ts:a.Ts:a.N*a.Ts, a.fval(:,3), 'DisplayName', 'Input cost', 'LineWidth',1.5,'color','b', 'LineStyle','--', 'LineStyle', ':');
plot(a.Ts:a.Ts:a.N*a.Ts, a.fval(:,5), 'DisplayName', 'Traffic cost', 'LineWidth',1.5,'color','b','LineStyle', '-');
plot(a.Ts:a.Ts:a.N*a.Ts, a.fval(:,6), 'DisplayName', 'Run-off cost', 'LineWidth',1.5,'color','r', 'LineStyle', '-');
ylim([0,0.025]);
xlabel('time(s)'); ylabel('cost(-)');
xlim([a.Ts,a.N*a.Ts]);
set(gca,'FontSize',14);
%legend('Location','southoutside', 'FontSize',11, 'Orientation','horizontal','NumColumns',2);

%% B plot - test II
subplot(3,3,2);
% error plot
plot(b.Ts:b.Ts:b.N*b.Ts, b.errorInstinctivePathToPath, 'DisplayName','Planned lane offset', 'LineWidth',1.5, 'Color','k', 'LineStyle','--');
grid on; hold on;
plot(b.Ts:b.Ts:b.N*b.Ts, b.error, 'DisplayName','Real lane offset', 'LineWidth',1.5,'Color','k');
ylim([-1,1]); xlabel('time(s)'); ylabel('$\delta(m)$');
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
xlim([b.Ts,b.N*b.Ts]);
set(gca,'FontSize',14);
legend('Location','best', 'FontSize',11);

title("\textbf{Test II.}");


subplot(3,3,5);
% input plot
plot(b.Ts:b.Ts:b.N*b.Ts, b.u_saved, 'DisplayName', 'Steering angle', 'LineWidth',1.5, 'color', 'k');
grid on;
xlabel('time(s)'); ylabel('$\alpha_f(rad)$');
ylim([-0.05,0.05]); xlim([b.Ts,b.N*b.Ts]);
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
set(gca,'FontSize',14);
%legend('Location','best', 'FontSize',11);

subplot(3,3,8);
% cost plot
plot(b.Ts:b.Ts:b.N*b.Ts, b.fval(:,1), 'DisplayName', 'Total', 'LineWidth',2,'color','k');
hold on; grid on;
plot(b.Ts:b.Ts:b.N*b.Ts, b.fval(:,2), 'DisplayName', 'Output error', 'LineWidth',1.5,'color','k', 'LineStyle', '--');
plot(b.Ts:b.Ts:b.N*b.Ts, b.fval(:,3), 'DisplayName', 'Input', 'LineWidth',1.5,'color','b', 'LineStyle','--', 'LineStyle', ':');
plot(b.Ts:b.Ts:b.N*b.Ts, b.fval(:,5), 'DisplayName', 'Traffic', 'LineWidth',1.5,'color','b','LineStyle', '-');
plot(b.Ts:b.Ts:b.N*b.Ts, b.fval(:,6), 'DisplayName', 'Run-off', 'LineWidth',1.5,'color','r', 'LineStyle', '-');
ylim([0,0.025]);
xlabel('time(s)'); ylabel('cost(-)');
xlim([b.Ts,b.N*b.Ts]);
set(gca,'FontSize',14);
legend('Location','north', 'FontSize',11, 'Orientation','horizontal','NumColumns',2);

%% C plot - test III.
subplot(3,3,3);
% error plot
plot(c.Ts:c.Ts:c.N*c.Ts, c.errorInstinctivePathToPath, 'DisplayName','Planned lane offset', 'LineWidth',1.5, 'Color','k', 'LineStyle','--');
grid on; hold on;
plot(c.Ts:c.Ts:c.N*c.Ts, c.error, 'DisplayName','Real lane offset', 'LineWidth',1.5,'Color','k');
ylim([-1,1]); xlabel('time(s)'); ylabel('$\delta(m)$');
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
xlim([c.Ts,c.N*c.Ts]);
set(gca,'FontSize',14);
legend('Location','best', 'FontSize',11);
title("\textbf{Test III.}");

subplot(3,3,6);
% input plot
plot(c.Ts:c.Ts:c.N*c.Ts, c.u_saved, 'DisplayName', 'Steering angle', 'LineWidth',1.5, 'color', 'k');
grid on;
xlabel('time(s)'); ylabel('$\alpha_f(rad)$');
ylim([-0.05,0.05]); xlim([c.Ts,c.N*c.Ts]);
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
set(gca,'FontSize',14);
%legend('Location','best', 'FontSize',11);

subplot(3,3,9);
% cost plot
plot(c.Ts:c.Ts:c.N*c.Ts, c.fval(:,1), 'DisplayName', 'Total cost', 'LineWidth',2,'color','k');
hold on; grid on;
plot(c.Ts:c.Ts:c.N*c.Ts, c.fval(:,2), 'DisplayName', 'Output error cost', 'LineWidth',1.5,'color','k', 'LineStyle', '--');
plot(c.Ts:c.Ts:c.N*c.Ts, c.fval(:,3), 'DisplayName', 'Input cost', 'LineWidth',1.5,'color','b', 'LineStyle','--', 'LineStyle', ':');
plot(c.Ts:c.Ts:c.N*c.Ts, c.fval(:,5), 'DisplayName', 'Traffic cost', 'LineWidth',1.5,'color','b','LineStyle', '-');
plot(c.Ts:c.Ts:c.N*c.Ts, c.fval(:,6), 'DisplayName', 'Run-off cost', 'LineWidth',1.5,'color','r', 'LineStyle', '-');
ylim([0,0.5]);
xlabel('time(s)'); ylabel('cost(-)');
xlim([c.Ts,c.N*c.Ts]);
set(gca,'FontSize',14);
annotation('textbox', [0.755, 0.22, 0.1, 0.1], 'String', "Different y-axis scale!", 'BackgroundColor', 'y');
%legend('Location','eastoutside', 'FontSize',11, 'Orientation','vertical','NumColumns',1);
