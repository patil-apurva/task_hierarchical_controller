%time parameters
T = 5;
h = 0.01;
t0 = 0;

%gains
k1 = 10; %gain for task 1
k2 = eye(2);

%obstacle radius and center
r = 1;
cx = 3;
cy = 0;
c = [cx; cy];

% r = 8;
% cx = 20;
% cy = -5;
% c = [cx; cy];

%goal position
pg = [6; 0];
% pg = [40; -3];

pgx = pg(1);
pgy = pg(2);

%when to activate v1
threshold = 1.5;
% threshold = 10;

s2 = 1; %Sigma=s*g_tilde, s2=s^2
alpha = 1; %R=alpha.I
b = 1; %V=b.||x(t)-goal||^2
d = 1; %psi=d.||x(T)-goal||^2

eta = 1; %lagrange multiplier
runs = 10000; %MC runs
traj_num = 1;

x0 = [0; 0];
v1_0 = [0; 0];
N_J1 = eye(2);

lambda = alpha*s2; %PDE linearization constant
s = sqrt(s2); %Sigma=s.I

figure(4)
hold on;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal