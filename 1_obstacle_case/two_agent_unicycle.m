close all;
clear;

%time parameters
T = 10;
h = 0.1;

% gains for u3 - maintain distance ell between the agents
% kp3 = 0.95;
% ki3 = 0;
% kd3 = 2.0015;
kp3 = 10;
ki3 = 0;
kd3 = 5;

%errors for u3
error3_integral = 0;
last_error3 = 0;

%gains for u2 - centroid should reach to the goal position
% kp2 = 0.02;
% ki2 = 0;
% kd2 = 0;
kp2 = 0.05;
ki2 = 0;
kd2 = 0;

%errors for u2
error2_integral = zeros(2,1);
last_error2 = zeros(2,1);

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

%gain for the drift term
e_tilde = 0;
k1 = -e_tilde/T;

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

p1 = p1_0;
p2 = p2_0;

x1 = x1_0;
y1 = y1_0;
v1 = v1_0;
theta_1 = theta_1_0;

x2 = x2_0;
y2 = y2_0;
v2 = v2_0;
theta_2 = theta_2_0;

P1 = [p1];
P2 = [p2];

agents_dist = []; %distance between 2 agents
dist_from_goal = []; %distance between centroid of two agents and goal
desired_dist_bet_agents = [];
centroid = []; %centorid of the two agents

for t=0:h:T-h
   %compute u3
    K3 = [(x1 - x2)*cos(theta_1) + (y1 - y2)*sin(theta_1), ...
         (x2 - x1)*v1*sin(theta_1) + (y1 - y2)*v1*cos(theta_1), ...
         (x2 - x1)*cos(theta_2) + (y2 - y1)*sin(theta_2), ...
         (x1 - x2)*v2*sin(theta_2) + (y2 - y1)*v2*cos(theta_2)];
     
    K3_dagger = K3.'/(K3*K3.');
  
    sigma3 = 0.5*((x1 - x2)^2 + (y1 - y2)^2);
    sigma3_desire = ell^2/2;
  
    error3 = sigma3_desire - sigma3;
    error3_integral = error3_integral + error3*h; %for I controller
    error3_derivative = (error3 - last_error3) / h; %for D controller
  
    error3_tilde = v1^2 + v2^2 - 2*v1*v2*(cos(theta_1)*cos(theta_2) + sin(theta_1)*sin(theta_2));
  
    u3 = K3_dagger*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde);
%     u3 = K3.'*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde);
  
    last_error3 = error3;
    N4 = eye(4) - K3_dagger*K3;
  
    %compute u2
    K2 = 0.5*[cos(theta_1) -v1*sin(theta_1) cos(theta_2) -v2*sin(theta_2); ...
            sin(theta_1) v1*cos(theta_1) sin(theta_2) v2*cos(theta_2)];
        
    K2_dagger = K2.'/(K2*K2.');
  
    sigma2 = [(x1 + x2)/2; (y1 + y2)/2];
    sigma2_desire = pg;
    
    error2 = sigma2_desire - sigma2;
    error2_integral = error2_integral + error2*h; %for I controller
    error2_derivative = (error2 - last_error2) / h; %for D controller
  
    u2 = K2_dagger*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
%     u2 = K2.'*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
    
    last_error2 = error2;  
    N3 = eye(4) - K2_dagger*K2;
  
    %compute u1 for agent 1
    uc_obs_angle_1 = atan((cy - y1)/(cx - x1)); %angle of the line connecting the unicycle 1 and the obstacle
    is_uc_in_ObsDir_1 = cos(uc_obs_angle_1 - theta_1) > 0; %is the unicycle 1 in the obstacle direction

    if(norm(p1 - c)<threshold && is_uc_in_ObsDir_1)
        K1_1 = [((x1 - cx)/norm(p1 - c))*cos(theta_1) + ((y1 - cy)/norm(p1 - c))*sin(theta_1), ...  
                ((y1 - cy)*v1/norm(p1 - c))*cos(theta_1) - ((x1 - cx)*v1/norm(p1 - c))*sin(theta_1)];
            
        K1_1_dagger = K1_1.' / (K1_1 * K1_1.');
        
        error1_1 = r_ext - norm(p1 - c); %for P controller
        error1_1_derivative = (error1_1 - last_error1_1) / h; %for D controller
        error1_1_integral = error1_1_integral + error1_1*h; %for I controller
        
        error1_1_tilde = ((y1 - cy)^2*v1^2*(cos(theta_1))^2 + (x1 - cx)^2*v1^2*(sin(theta_1))^2)/((norm(p1-c))^3);
        
        u1_1 = K1_1_dagger*((kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative) - error1_1_tilde);
%         u1_1 = K1_1.'*((kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative) - error1_1_tilde);
        
        last_error1_1 = error1_1; 
    else
        u1_1 = zeros(2,1);
        K1_1 = zeros(1,2);
        K1_1_dagger = zeros(2,1);
    end
 %================================================================    
    %compute u1 for agent 2
    uc_obs_angle_2 = atan((cy - y2)/(cx - x2)); %angle of the line connecting the unicycle 2 and the obstacle
    is_uc_in_ObsDir_2 = cos(uc_obs_angle_2 - theta_2) > 0; %is the unicycle 2 in the obstacle direction

    if(norm(p2 - c)<threshold && is_uc_in_ObsDir_2)
        K1_2 = [((x2 - cx)/norm(p2 - c))*cos(theta_2) + ((y2 - cy)/norm(p2 - c))*sin(theta_2), ...  
                ((y2 - cy)*v2/norm(p2 - c))*cos(theta_2) - ((x2 - cx)*v2/norm(p2 - c))*sin(theta_2)];
            
        K1_2_dagger = K1_2.' / (K1_2 * K1_2.');
        
        error1_2 = r_ext - norm(p2 - c); %for P controller
        error1_2_derivative = (error1_2 - last_error1_2) / h; %for D controller
        error1_2_integral = error1_2_integral + error1_2*h; %for I controller
        
        error1_2_tilde = ((y2 - cy)^2*v2^2*(cos(theta_2))^2 + (x2 - cx)^2*v2^2*(sin(theta_2))^2)/((norm(p2-c))^3);
        
        u1_2 = K1_2_dagger*((kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative) - error1_2_tilde);
%         u1_2 = K1_2.'*((kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative) - error1_2_tilde);
        
        last_error1_2 = error1_2; 
    else
        u1_2 = zeros(2,1);
        K1_2 = zeros(1,2);
        K1_2_dagger = zeros(2,1);
    end
  
    K1 = blkdiag(K1_1, K1_2);
    K1_dagger = blkdiag(K1_1_dagger, K1_2_dagger);
    
    N2 = eye(4) - K1_dagger*K1;
  
    u1 = [u1_1; u1_2];
    
    u = u1 + N2*(u2 + N3*u3); %with task hierarchy
%     u = u1 + u2 + u3; %without task hierarchy
%     u = u2 + N3*u3;

%     u = u1 + N2*(u3 + N4*u2); %with task hierarchy
%     u = u3 + N4*u2;

    a1 = u(1);
    omega_1 = u(2);
    a2 = u(3);
    omega_2 = u(4);
    
    % integrate the dynamics forward
    x1 = x1 + k1*x1*h + v1*cos(theta_1)*h;
    y1 = y1 + k1*(y1)*h + v1*sin(theta_1)*h;
    v1 = v1 + a1*h;
    theta_1 = theta_1 + omega_1*h;
    
    x2 = x2 + k1*x2*h + v2*cos(theta_2)*h;
    y2 = y2 + k1*(y2)*h + v2*sin(theta_2)*h;
    v2 = v2 + a2*h;
    theta_2 = theta_2 + omega_2*h;
    
    p1 = [x1; y1];
    p2 = [x2; y2];
    
    P1 = [P1, p1];
    P2 = [P2, p2];
    
    centroid = [centroid, (p1 + p2)/2];
    agents_dist = [agents_dist, norm(p1 - p2)];
    dist_from_goal = [dist_from_goal, norm((p1 + p2)/2 - pg)];
    desired_dist_bet_agents = [desired_dist_bet_agents, ell];
        
end

figure(5)
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
plot(pg(1), pg(2), '.', 'MarkerSize',100,'MarkerEdgeColor','g','MarkerFaceColor','g')
plot(P1(1,:), P1(2,:), 'k', 'LineWidth', 2)
plot(P2(1,:), P2(2,:), 'b', 'LineWidth', 2)
plot (centroid(1, :), centroid(2, :), '--m', 'LineWidth',2);
legend('obstacle', 'goal position', 'agent1', 'agent2', 'centroid')
grid on

figure(6)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
plot(0:h:T-h, agents_dist, 'b', 'LineWidth', 2)
grid on
% title("Distance from goal");
xlabel("$t$", 'Interpreter','latex', 'FontSize', 30);
ylabel("Dist. between agents", 'Interpreter','latex', 'FontSize', 30)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';
plot(0:h:T-h, desired_dist_bet_agents, '--k', 'LineWidth', 2);
legend('Dist. between agents', 'Desired dist. between agents');
 
figure(7)
hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 18)
plot(0:h:T-h, dist_from_goal, 'b', 'LineWidth', 2)
grid on
% title("Distance from goal");
xlabel("$t$", 'Interpreter','latex', 'FontSize', 30);
ylabel("Distance from goal", 'Interpreter','latex', 'FontSize', 30)
set(gca,'LineWidth',1)
ax = gca;
ax.LineWidth = 1;
ax.Color = 'w';

% save("for_plotting_PID.mat","agents_dist", "dist_from_goal", "T", "h")