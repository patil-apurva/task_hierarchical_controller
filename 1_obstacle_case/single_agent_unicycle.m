close all;
clear;

%time parameters
T = 10;
h = 0.01;

%gains
% ku = 1; %ultimate gain
% Tu = 6.3; %oscillation time

kp1 = 3; %proportional gain for task 1
ki1 = 0; %integral gain for task 1
kd1 = 0; %derivative gain for task 1

% kp2 = 0.6*ku*[1, 0; 0, 1]; %proportional gain for task 2
% ki2 = (1.2*ku/Tu)*[1, 0; 0, 1]; %integral gain for task 2
% kd2 = (0.075*ku*Tu)*[1, 0; 0, 1]; %derivative gain for task 2

kp2 = 0.04; %proportional gain for task 2
ki2 = 0; %integral gain for task 2
kd2 = 0.05; %derivative gain for task 2

%obstacle radius and center
r = 8;
r_ext = 12;
cx = 20;
cy = 0;
c = [cx; cy];

%goal position
pg = [40; 0];

%when to activate u1
threshold = 10;

%initial state
% x_0 = 11.38;
% y_0 = -2.32;
% v_0 = 0.001;
% theta_0 = -77*pi/180;

x_0 = 0;
y_0 = 0;
v_0 = 1;
theta_0 = 0;

q_0 = [x_0; y_0; v_0; theta_0];
p_0 = [x_0; y_0];

is_task1_active = 0; 

q_t = q_0;
q_t_w_flag = [q_t; is_task1_active];
p_t = p_0;

Q_t = [q_t];
Q_t_w_flag = [q_t_w_flag];
P_t = [p_t];

last_error1 = 0; %for D controller of task 1
error1_integral = 0; %for I controller of task 1

last_error2 = [0; 0]; %for D controller of task 2
error2_integral = [0; 0]; %for I controller of task 2

ERROR2 = [];
OMEGA = [];
dist_from_goal = []; %distance between the agent and goal

for t=0:h:T-h
    
    %current state
    x_t = q_t(1);
    y_t = q_t(2);
    v_t = q_t(3);
    theta_t = q_t(4);
    
%     v_t = min(v_t,0.1);
    %compute u2
    K2 = [cos(theta_t), -v_t*sin(theta_t); sin(theta_t), v_t*cos(theta_t)];
    K2_dagger = K2.'/(K2*K2.');
    
    error2 = pg - p_t; %for P controller
    error2_derivative = (error2 - last_error2) / h; %for D controller
    error2_integral = error2_integral + error2*h; %for I controller
    
    u2 = K2_dagger*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
%     u2 = K2.'*(kp2*error2 + kd2*error2_derivative);
    
    last_error2 = error2; 
    ERROR2 = [ERROR2, error2];

    uc_obs_angle = atan((cy - y_t)/(cx - x_t)); %angle of the line connecting the unicycle and the obstacle
    is_uc_in_ObsDir = cos(uc_obs_angle - theta_t) > 0; %is the unicycle in the obstacle direction
     
%     %compute u1
    if(norm(p_t - c)<threshold && is_uc_in_ObsDir)
        K1 = [((x_t - cx)/norm(p_t-c))*cos(theta_t) + ((y_t - cy)/norm(p_t-c))*sin(theta_t) , ((y_t - cy)*v_t/norm(p_t-c))*cos(theta_t) - ((x_t - cx)*v_t/norm(p_t-c))*sin(theta_t)];
        K1_dagger = K1.' / (K1 * K1.');
        
        error1 = r_ext - norm(p_t - c); %for P controller
        error1_derivative = (error1 - last_error1) / h; %for D controller
        error1_integral = error1_integral + error1*h; %for I controller
        
        error1_tilde = ((y_t - cy)^2*v_t^2*(cos(theta_t))^2 + (x_t - cx)^2*v_t^2*(sin(theta_t))^2)/((norm(p_t-c))^3);
        
        u1 = K1_dagger*((kp1*error1 + ki1*error1_integral + kd1*error1_derivative) - error1_tilde);
%         u1 = K1.'*((kp1*error1 + ki1*error1_integral + kd1*error1_derivative) - error1_tilde);
        
        last_error1 = error1; 
         
        N_K1 = eye(2) - K1_dagger*K1;
%         N_K1 = eye(2) - K1.'*K1;
        
        J1 = [(x_t - cx)/norm(p_t - c), (y_t - cy)/norm(p_t - c), 0];
        J1_dagger = J1.' / (J1 * J1.');
        N_J1 = eye(3) - J1_dagger*J1;
        is_task1_active = 1;
    else
        u1 = zeros(2,1);
        N_K1 = eye(2);
        N_J1 = eye(3);
        is_task1_active = 0;
    end
    
    ud = u1 + N_K1*u2; %with task hierarchy
%     ud = u1 + N_J1(1:2, 1:2)*u2; %with task hierarchy
%     ud = u1 + u2; %without task hierarchy
%     ud = u2;
    
    ad = ud(1);
    omega_d = ud(2);
    
    OMEGA = [OMEGA, omega_d*180/pi];
    
    % integrate the dynamics forward
    x_tp1 = x_t + v_t*cos(theta_t)*h;
    y_tp1 = y_t + v_t*sin(theta_t)*h;
    v_tp1 = v_t + ad*h;
    theta_tp1 = theta_t + omega_d*h;
    
    x_t = x_tp1;
    y_t = y_tp1;
    v_t = v_tp1;
    theta_t = theta_tp1;
    
    q_t = [x_t; y_t; v_t; theta_t];
    q_t_w_flag = [q_t; is_task1_active];
    
    p_t = [x_t; y_t];
    Q_t = [Q_t, q_t];
    
    Q_t_w_flag = [Q_t_w_flag, q_t_w_flag];
    P_t = [P_t, p_t];
    
    dist_from_goal = [dist_from_goal, norm(p_t - pg)];
end

Q_t(4,:) = Q_t(4,:)*180/pi;
Q_t_w_flag(4,:) = Q_t_w_flag(4,:)*180/pi;

figure(5)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
xlabel('$x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$y$', 'Interpreter','latex','FontSize', 30); 
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 2;
ax.Color = 'w';
axis equal;
plot(pg(1), pg(2), '.', 'MarkerSize',100,'MarkerEdgeColor','g','MarkerFaceColor','g')
plot(P_t(1,:), P_t(2,:), 'b', 'LineWidth', 2)
legend('obstacle', 'goal position', 'agent')
grid on

figure(6)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
plot(0:h:T-h, ERROR2(1,:), 'b', 'LineWidth', 2)
grid on
% title("Distance from goal");
xlabel("$t$", 'Interpreter','latex', 'FontSize', 30);
ylabel("Distance from goal", 'Interpreter','latex', 'FontSize', 30)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
% axis equal;

% figure(7)
% hold on;
% hold on;
% set(gca, 'FontName', 'Arial', 'FontSize', 18)
% plot(0:h:T-h, ERROR2(2,:))
% grid on
% title("error in y");
% xlabel("t");
% ylabel("error in y")
% set(gca,'LineWidth',1)
% ax = gca;
% ax.LineWidth = 1;
% ax.Color = 'w';
% % axis equal;

% figure(8)
% hold on;
% plot(OMEGA)
% grid on
% title("omega");
% xlabel("t");
% ylabel("omega")
% 
% figure(9)
% hold on;
% plot(Q_t(4,:))
% grid on
% title("theta");
% xlabel("t");
% ylabel("theta")


save("for_plotting_PID_single_agent.mat", "dist_from_goal", "T", "h")