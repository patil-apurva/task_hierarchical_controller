clear 
close all
load('for_plotting_PID_single_agent')
load('for_plotting_PI_single_agent')

figure(12)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';

curve1 = dist_from_goal_PI_mean + dist_from_goal_PI_std;
curve2 = dist_from_goal_PI_mean - dist_from_goal_PI_std;

t_steps = 0:h:T-h;
t_steps2 = [t_steps, fliplr(t_steps)];

inBetween = [curve1, fliplr(curve2)];
lightBlue = [91, 207, 244] / 255;
fill(t_steps2, inBetween, lightBlue, 'HandleVisibility', 'off');

% plot(0:h:T-h, dist_from_goal, 'r', 'Linewidth', 2)

plot(t_steps, dist_from_goal_PI_mean, 'b', 'LineWidth', 2);
xlabel("$t$", 'Interpreter','latex', 'FontSize', 30);
ylabel('Distance from goal', 'Interpreter','latex', 'FontSize', 30)
% legend ('Without PI', 'With PI')
grid on