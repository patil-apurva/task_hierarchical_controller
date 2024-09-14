%time parameters
T = 10;
h = 0.1;

% gains for u3 - maintain distance ell between the agents
% kp3 = 0.95;
% ki3 = 0;
% kd3 = 2.0015;
kp3 = 1;
ki3 = 0;
kd3 = 1;
 
%errors for u3
error3_integral = 0;
last_error3 = 0;

% %gains for u2 - centroid should reach to the goal position
% kp2 = 0.025;
% ki2 = 0;
% kd2 = 0;
% 
% %errors for u2
% error2_integral = zeros(2,1);
% last_error2 = zeros(2,1);

%gains for u1 - avoid collision
kp1 = 3;
ki1 = 0;
kd1 = 0;

%errors for u1 - agent 1
error1_1_integral = 0;
last_error1_1 = 0;

%errors for u1 - agent 2
error1_2_integral = 0;
last_error1_2 = 0;

%obstacle radius and center
r = 8;
r_ext = 12;
cx = 20;
cy = 0;
c = [cx; cy];

%when to activate u1
threshold = 10;

%goal position
pg = [40; -5];

pgx = pg(1);
pgy = pg(2);

%initial states of two agents
x1_0 = 0;
y1_0 = -7;
v1_0 = 1;
theta_1_0 = 0;

x2_0 = 0;
y2_0 = -3;
v2_0 = 1;
theta_2_0 = 0;

%desired dist between the two agents
ell = 4;

p1_0 = [x1_0; y1_0];
p2_0 = [x2_0; y2_0];

s2 = 0.1; %Sigma=s*g_tilde, s2=s^2
alpha = 10; %R=alpha.I % try reducing this
b = 0.04; %V=b.||x(t)-goal||^2
e = 0; %cost associated with maintaining a specified distance between agents
eta1 = 0; %agent 1's cost associated with colliding with obstacle
eta2 = 0; %agent 2's cost associated with colliding with obstacle

runs = 10000; %MC runs
traj_num = 1; %number of trajectroies to plot

lambda = alpha*s2; %PDE linearization constant
s = sqrt(s2); %Sigma=s.G

figure(2)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
xlabel('$x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$y$', 'Interpreter','latex','FontSize', 30); 
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
axis equal;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal