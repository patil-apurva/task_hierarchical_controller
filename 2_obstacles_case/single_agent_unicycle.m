close all;
clear;

%time parameters
T = 20;
h = 0.1;
t0 = 0;

% gains for acc_3
kp3 = 0.006;
ki3 = 0;
kd3 = 0;

last_error3 = [0; 0]; %for D controller of task 2
error3_integral = [0; 0]; %for I controller of task 2

%gains for acc_1
kp1 = 1;
ki1 = 0;
kd1 = 0;

%gains for acc_2
kp2 = kp1;
ki2 = ki1;
kd2 = kd1;

%errors for acc_1 (obs A) - agent 1
error1_integral = 0;
last_error1 = 0;

%errors for acc_2 (obs B) - agent 1
error2_integral = 0;
last_error2 = 0;

%obstacles radius and center
rA = 1;
r_extA = 1.5;
cxA = 3;
cyA = 0.5;
cA = [cxA; cyA];

rB = 1;
r_extB = 1.5;
cxB = 3;
cyB = -0.5;
cB = [cxB; cyB];

%goal position
pg = [7; 0];

pgx = pg(1);
pgy = pg(2);

%when to activate v1 and v2
threshold = 1.5;

x_0 = 0;
y_0 = 0;
v_0 = 0;
theta_0 = 0;

q_0 = [x_0; y_0; v_0; theta_0];
p_0 = [x_0; y_0];

q_t = q_0;
p_t = p_0;

Q_t = [q_t];
P_t = [p_t];

ERROR3 = [];
dist_from_goal = []; %distance between the agent and goal

for t=0:h:T-h
    
    %current state
    x_t = q_t(1);
    y_t = q_t(2);
    v_t = q_t(3);
    theta_t = q_t(4);
    
    %compute u3
    K3 = [cos(theta_t), -v_t*sin(theta_t); sin(theta_t), v_t*cos(theta_t)];
%     K3_dagger = K3.'/(K3*K3.');
    
    error3 = pg - p_t; %for P controller
    error3_derivative = (error3 - last_error3) / h; %for D controller
    error3_integral = error3_integral + error3*h; %for I controller
    
%     u3 = K3_dagger*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
    u3 = K3.'*(kp3*error3 + kd3*error3_derivative);
    
    last_error3 = error3; 
    ERROR3 = [ERROR3, error3];
    
    %compute u1 (obs A)
    uc_obs_angle = atan((cyA - y_t)/(cxA - x_t)); %angle of the line connecting the unicycle and the obstacle
    is_uc_in_ObsDir = cos(uc_obs_angle - theta_t) > 0; %is the unicycle in the obstacle direction
    
    if(norm(p_t - cA)<threshold && is_uc_in_ObsDir)
        K1 = [((x_t - cxA)/norm(p_t - cA))*cos(theta_t) + ((y_t - cyA)/norm(p_t - cA))*sin(theta_t) , ((y_t - cyA)*v_t/norm(p_t - cA))*cos(theta_t) - ((x_t - cxA)*v_t/norm(p_t - cA))*sin(theta_t)];
        K1_dagger = K1.' / (K1 * K1.');
        
        error1 = r_extA - norm(p_t - cA); %for P controller
        error1_derivative = (error1 - last_error1) / h; %for D controller
        error1_integral = error1_integral + error1*h; %for I controller
        
        error1_tilde = ((y_t - cyA)^2*v_t^2*(cos(theta_t))^2 + (x_t - cxA)^2*v_t^2*(sin(theta_t))^2)/((norm(p_t - cA))^3);
        
%         u1 = K1_dagger*((kp1*error1 + ki1*error1_integral + kd1*error1_derivative) - error1_tilde);
        u1 = K1.'*((kp1*error1 + ki1*error1_integral + kd1*error1_derivative) - error1_tilde);
        
        last_error1 = error1; 
         
        N_K1 = eye(2) - K1_dagger*K1;
%         N_K1 = eye(2) - K1.'*K1;
    else
        u1 = zeros(2,1);
        N_K1 = eye(2);
    end
    
    %compute u2 (obs B)
    uc_obs_angle = atan((cyB - y_t)/(cxB - x_t)); %angle of the line connecting the unicycle and the obstacle
    is_uc_in_ObsDir = cos(uc_obs_angle - theta_t) > 0; %is the unicycle in the obstacle direction
    
    if(norm(p_t - cB)<threshold && is_uc_in_ObsDir)
        K2 = [((x_t - cxB)/norm(p_t - cB))*cos(theta_t) + ((y_t - cyB)/norm(p_t - cB))*sin(theta_t) , ((y_t - cyB)*v_t/norm(p_t - cB))*cos(theta_t) - ((x_t - cxB)*v_t/norm(p_t - cB))*sin(theta_t)];
        K2_dagger = K2.' / (K2 * K2.');
        
        error2 = r_extB - norm(p_t - cB); %for P controller
        error2_derivative = (error2 - last_error2) / h; %for D controller
        error2_integral = error2_integral + error2*h; %for I controller
        
        error2_tilde = ((y_t - cyB)^2*v_t^2*(cos(theta_t))^2 + (x_t - cxB)^2*v_t^2*(sin(theta_t))^2)/((norm(p_t - cB))^3);
        
%         u2 = K2_dagger*((kp2*error2 + ki2*error2_integral + kd2*error2_derivative) - error2_tilde);
        u2 = K2.'*((kp2*error2 + ki2*error2_integral + kd2*error2_derivative) - error2_tilde);
        
        last_error2 = error2; 
         
        N_K2 = eye(2) - K2_dagger*K2;
%         N_K2 = eye(2) - K2.'*K2;
    else
        u2 = zeros(2,1);
        N_K2 = eye(2);
    end
    
    ud = u1 + N_K1*u2 + N_K1*N_K2*u3; %with task hierarchy
%     ud = u3; 
    
    ad = ud(1);
    omega_d = ud(2);
    
     % integrate the dynamics forward
    x_t = x_t + v_t*cos(theta_t)*h;
    y_t = y_t + v_t*sin(theta_t)*h;
    v_t = v_t + ad*h;
    theta_t = theta_t + omega_d*h;
    
    q_t = [x_t; y_t; v_t; theta_t];
    
    p_t = [x_t; y_t];
    Q_t = [Q_t, q_t];
    P_t = [P_t, p_t];
    
    dist_from_goal = [dist_from_goal, norm(p_t - pg)];
end

figure(5)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
p = nsidedpoly(1000, 'Center', [cxA cyA], 'Radius', rA);
plot(p, 'FaceColor', 'r')
p = nsidedpoly(1000, 'Center', [cxB cyB], 'Radius', rB);
plot(p, 'FaceColor', 'r')
xlabel('$x$', 'Interpreter','latex', 'FontSize', 30); ylabel('$y$', 'Interpreter','latex','FontSize', 30); 
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
axis equal;
plot(pg(1), pg(2), '.', 'MarkerSize',100,'MarkerEdgeColor','g','MarkerFaceColor','g')
plot(P_t(1,:), P_t(2,:), 'k', 'LineWidth', 1)
legend('obstacle', 'goal position', 'agent')
grid on

figure(6)
hold on;
plot(ERROR3(1,:))
grid on
title("error in x");
xlabel("t");
ylabel("error in x")

figure(7)
hold on;
plot(ERROR3(2,:))
grid on
title("error in y");
xlabel("t");
ylabel("error in y")
