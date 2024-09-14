%time parameters
T = 20;
h = 0.1;
t0 = 0;

% gains for acc_3
kp3 = 0.03031;
ki3 = 0;
kd3 = 0.344;

%errors for acc_3
error3_integral = zeros(2,1);
last_error3 = zeros(2,1);

% % gains for acc_4
% kp4 = 0.2031;
% ki4 = 0;
% kd4 = 0.744;

% gains for acc_4
kp4 = 0.7031;
ki4 = 0;
kd4 = 0;

%errors for acc_4
error4_integral = 0;
last_error4 = 0;

%gains for acc_1
kp1 = 0.7;
ki1 = 0;
kd1 = 1;

%gains for acc_2
kp2 = kp1;
ki2 = ki1;
kd2 = kd1;

%errors for acc_1 (obs A) - agent 1
error1_1_integral = 0;
last_error1_1 = 0;

%errors for acc_1 (obs A) - agent 2
error1_2_integral = 0;
last_error1_2 = 0;

%errors for acc_2 (obs B) - agent 1
error2_1_integral = 0;
last_error2_1 = 0;

%errors for acc_2 (obs B) - agent 2
error2_2_integral = 0;
last_error2_2 = 0;

%obstacles radius and center
rA = 1;
r_extA = 1.5;
% r_ext = r;
cxA = 3;
cyA = 1.4;
cA = [cxA; cyA];

rB = 1;
r_extB = 1.5;
% r_ext = r;
cxB = 3;
cyB = -1.4;
cB = [cxB; cyB];

%goal position
pg = [7; 0];

pgx = pg(1);
pgy = pg(2);

%when to activate v1 and v2
threshold = 1.5;

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

figure(4)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
xlabel('$q_x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$q_y$', 'Interpreter','latex','FontSize', 30); 
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
axis equal;
p = nsidedpoly(1000, 'Center', [cxA cyA], 'Radius', rA);
plot(p, 'FaceColor', 'r','HandleVisibility','off')
p = nsidedpoly(1000, 'Center', [cxB cyB], 'Radius', rB);
plot(p, 'FaceColor', 'r')
axis equal