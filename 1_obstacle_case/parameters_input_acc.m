%time parameters
T = 5;
h = 0.01;
t0 = 0;

%gains for task 1
kp1 = 0.035;
ki1 = 0.000306;
kd1 = 1;

% gains for acc_2
kp2 = 0.000731;
ki2 = 0;
kd2 = 0.144;

error1_integral = 0;
last_error1 = 0;

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

s2 = 1; %Sigma=s*g_tilde, s2=s^2
alpha = 1; %R=alpha.I
b = 0.1; %V=b.||x(t)-goal||^2
d = 1; %psi=d.||x(T)-goal||^2

eta = 0.1; %lagrange multiplier
runs = 10000; %MC runs
traj_num = 10;

q0 = [0; 0];
v0 = [0; 0];
a1_0 = [0; 0];
N2 = eye(2);

lambda = alpha*s2; %PDE linearization constant
s = sqrt(s2); %Sigma=s.I

figure(4)
hold on;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal