%time parameters
T = 10;
h = 0.1;

% gains for u3 - maintain distance ell between the agents
kp3 = 10;
ki3 = 0;
kd3 = 5;
% kp3 = 0.5;
% ki3 = 0;
% kd3 = 0.2;
 
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
kp1 = 1.7;
ki1 = 0;
kd1 = 0;

% kp1 = 2;
% ki1 = 0;
% kd1 = 0;

%errors for u1 - agent 1
error1_1_integral = 0;
last_error1_1 = 0;

%errors for u1 - agent 2
error1_2_integral = 0;
last_error1_2 = 0;

%gain for the drift term
e_tilde = 0;
k1 = -e_tilde/T;

k2 = 0;
k3 = 0;
%obstacle radius and center
r = 1;
r_ext = 1.5;
cx = -2;
cy = 0.2;
c = [cx; cy];

%when to activate u1
threshold = 1.5;

%goal position
pg = [0; 0];

pgx = pg(1);
pgy = pg(2);

%initial states of two agents
x1_0 = -4.5;
y1_0 = 0;
v1_0 = 0.1;
theta_1_0 = -pi/4;

x2_0 = -4;
y2_0 = 0;
v2_0 = 0.1;
theta_2_0 = -pi/4;

%desired dist between the two agents
ell = 0.5;

p1_0 = [x1_0; y1_0];
p2_0 = [x2_0; y2_0];

s2 = 0.1; %Sigma=s*g_tilde, s2=s^2
alpha = 10; %R=alpha.I % try reducing this
b = 0.21; %V=b.||x(t)-goal||^2
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
plot(pgx, pgy, '.', 'MarkerSize',100,'MarkerEdgeColor','g','MarkerFaceColor','g')
axis equal