%time parameters
T = 10;
h = 0.01;
t0 = 0;

%gains
kp1 =  3; %proportional gain for task 1
ki1 = 0; %integral gain for task 1
kd1 = 0; %derivative gain for task 1

kp2 = 0.04; %proportional gain for task 2
ki2 = 0; %integral gain for task 2
kd2 = 0.05; %derivative gain for task 2

%errors for control input 1    
error1_integral = 0;
last_error1 = 0;

%errors for control input 2
last_error2 = [0; 0]; %for D controller of task 2
error2_integral = [0; 0]; %for I controller of task 2

%obstacle radius and center
r = 8;
r_ext = 12;
cx = 20;
cy = 0;
c = [cx; cy];

%when to activate u1
threshold = 10;

%goal position
pg = [40; 0];

pgx = pg(1);
pgy = pg(2);

s2 = 0.1; %Sigma=s*g_tilde, s2=s^2
alpha = 10; %R=alpha.I
b = 0.07; %V=b.||x(t)-goal||^2

eta = 0; %collision cost
runs = 10000; %MC runs
traj_num = 1;

x_0 = 0;
y_0 = 0;
v_0 = 1;
theta_0 = 0;

q_0 = [x_0; y_0; v_0; theta_0];
p_0 = [x_0; y_0];

lambda = alpha*s2; %PDE linearization constant
s = sqrt(s2); %Sigma=s.I

figure(4)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
xlabel('$x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$y$', 'Interpreter','latex','FontSize', 30); 
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
axis equal;