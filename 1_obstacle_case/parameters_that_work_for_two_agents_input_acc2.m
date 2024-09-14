%time parameters
T = 10;
h = 0.1;
t0 = 0;

% % gains for acc_2
% kp2 = 0.000731;
% ki2 = 0;
% kd2 = 0.144;
% 
% %errors for acc_2
% error2_integral = 0;
% last_error2 = 0;

% gains for acc_3
kp3 = 0.00731;
ki3 = 0;
kd3 = 0.144;

%errors for acc_3
error3_integral = 0;
last_error3 = 0;

%gains for acc_1
kp1 = 0.55;
ki1 = 0;
kd1 = 0.1;

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

%when to activate v1
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
ell = 2;

s2 = 1; %Sigma=s*g_tilde, s2=s^2
alpha = 1; %R=alpha.I
b = 0.01; %V=b.||x(t)-goal||^2
e = 0.01; %cost associated with maintaining a specified distance between agents
eta1 = 0.1; %agent 1's cost associated with colliding with obstacle
eta2 = 0.1; %agent 2's cost associated with colliding with obstacle

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
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal