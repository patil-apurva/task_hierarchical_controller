%time parameters
T = 20;
h = 0.1;
t0 = 0;

% gains for acc_2
kp2 = 0.03031;
ki2 = 0;
kd2 = 0.344;

%errors for acc_2
error2_integral = zeros(2,1);
last_error2 = zeros(2,1);

% % gains for acc_3
% kp3 = 0.2031;
% ki3 = 0;
% kd3 = 0.744;

% gains for acc_3
kp3 = 0.9031;
ki3 = 0;
kd3 = 0.844;

%errors for acc_3
error3_integral = 0;
last_error3 = 0;

% %gains for acc_1
% kp1 = 0.70;
% ki1 = 0;
% kd1 = 1;

%gains for acc_1
kp1 = 0.1;
ki1 = 0;
kd1 = 1;


%errors for acc_1 - agent 1
error1_1_integral = 0;
last_error1_1 = 0;

%errors for acc_1 - agent 2
error1_2_integral = 0;
last_error1_2 = 0;

%obstacle radius and center
r = 1;
r_ext = 2;
% r_ext = r;
cx = 3;
cy = 0;
c = [cx; cy];

%goal position
pg = [6; 0];

pgx = pg(1);
pgy = pg(2);

%when to activate a1
threshold = 2;

%initial states of two agents
qx1_0 = 0;
qy1_0 = 0.5;
vx1_0 = 0;
vy1_0 = 0;

qx2_0 = 0;
qy2_0 = -0.5;
vx2_0 = 0;
vy2_0 = 0;

%desired dist between the two agents
ell = 1;

s2 = 0.01; %Sigma=s*g_tilde, s2=s^2
alpha = 100; %R=alpha.I % try reducing this
b = 1; %V=b.||x(t)-goal||^2
e = 0; %cost associated with maintaining a specified distance between agents
eta1 = 0; %agent 1's cost associated with colliding with obstacle
eta2 = 0; %agent 2's cost associated with colliding with obstacle

runs = 10000; %MC runs
traj_num = 1; %number of trajectroies to plot

lambda = alpha*s2; %PDE linearization constant
s = sqrt(s2); %Sigma=s.G

%initial position
q1_0 = [qx1_0; qy1_0];
q2_0 = [qx2_0; qy2_0];
q_0 = [q1_0; q2_0];

v1_0 = [vx1_0; vy1_0];
v2_0 = [vx2_0; vy2_0];
v_0 = [v1_0; v2_0];

figure(1)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
xlabel('$q_x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$q_y$', 'Interpreter','latex','FontSize', 30); 
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
axis equal;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal