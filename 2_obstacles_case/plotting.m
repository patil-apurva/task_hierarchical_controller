clear 
close all
load('for_plotting_PID.mat')
load('for_plotting_PI.mat')

T = 20;
h = 0.1;
t0 = 0;

figure(4)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
plot(t0:h:T-h, agents_dist_PID, 'Linewidth', 3)
plot(t0:h:T-h, agents_dist_PI, 'Linewidth', 3)
xlabel('t')
ylabel('Distance between two agents')
plot(t0:h:T-h, desired_dist_bet_agents, '--', 'Linewidth', 3);
legend ('Without PI', 'With PI', 'Desired dist.')

figure(5)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
plot(t0:h:T-h, dist_from_goal_PID, 'Linewidth', 3)
plot(t0:h:T-h, dist_from_goal_PI, 'Linewidth', 3)
xlabel('t')
ylabel('Dist. between centroid and goal')
legend ('Without PI', 'With PI')