clearvars
close all
run('parameters_unicycle_2agents.m')
seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator

for traj_itr=1:1
    traj_itr
    
    P1 = []; %to store all positions of agent 1's trajectory
    P2 = []; %to store all positions of agent 2's trajectory
    
    P1 = [P1, p1_0]; %stack the initial position of agent 1
    P2 = [P2, p2_0]; %stack the initial position of agent 2
    
    %start the agent 1's state from the given initial state
    x1 = x1_0;
    y1 = y1_0;
    v1 = v1_0;
    theta_1 = theta_1_0;
    
    p1 = p1_0;

    %start the agent 2's state from the given initial state
    x2 = x2_0;
    y2 = y2_0;
    v2 = v2_0;
    theta_2 = theta_2_0;
    
    p2 = p2_0;
    
    agents_dist_PI = []; %distance between 2 agents
    dist_from_goal_PI = []; %distance between centroid of two agents and goal
    desired_dist_bet_agents = [];
    centroid = []; %centorid of the two agents
    
    for t = 0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
%          t
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
%         u3 = K3.'*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde); 
        
        last_error3 = error3;
        N4 = eye(4) - K3_dagger*K3;
        
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
        
        %compute u2
        K2 = 0.5*[cos(theta_1) -v1*sin(theta_1) cos(theta_2) -v2*sin(theta_2); ...
                sin(theta_1) v1*cos(theta_1) sin(theta_2) v2*cos(theta_2)];

        K2_dagger = K2.'/(K2*K2.');
        
        N3 = eye(4) - K2_dagger*K2;

%         sigma2 = [(x1 + x2)/2; (y1 + y2)/2];
%         sigma2_desire = pg;
% 
%         error2 = sigma2_desire - sigma2;
%         error2_integral = error2_integral + error2*h; %for I controller
%         error2_derivative = (error2 - last_error2) / h; %for D controller
%         
%         last_error2 = error2;  
%         
%         u2 = K2_dagger*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative); 
        
        eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
        eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt
        eps_t_all_3 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_3(t) at the start of each sample path starting at time t and state xt
        eps_t_all_4 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_4(t) at the start of each sample path starting at time t and state xt
        
        S_tau_all = arrayfun(@simulateMC_unicycle_2agents, eps_t_all_1, eps_t_all_2, eps_t_all_3, eps_t_all_4, ...
            x1, y1, v1, theta_1, x2, y2, v2, theta_2, ...
            u1_1(1), u1_1(2), u1_2(1), u1_2(2), u3(1), u3(2), u3(3), u3(4), ... 
            t, h, T, b, e, s, cx, cy, r, r_ext, threshold, pgx, pgy, ell, eta1, eta2, kp1, ki1, kd1, kp3, ki3, kd3,...
            N2(1,1), N2(1,2), N2(1,3), N2(1,4), N2(2,1), N2(2,2), N2(2,3), N2(2,4), N2(3,1), N2(3,2), N2(3,3), N2(3,4), N2(4,1), N2(4,2), N2(4,3), N2(4,4),...
            N3(1,1), N3(1,2), N3(1,3), N3(1,4), N3(2,1), N3(2,2), N3(2,3), N3(2,4), N3(3,1), N3(3,2), N3(3,3), N3(3,4), N3(4,1), N3(4,2), N3(4,3), N3(4,4));
        
        eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2; eps_t_all_3; eps_t_all_4]); %convert GPU array into a regular array
        
        denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size: (4 X 1))
        denom = sum(denom_i); %scalar
        
        cali_G = (N2.')/(N2*N2.');

        u2 = (s/sqrt(h))*cali_G*N2*numer/denom; %the control input
       
        if(any(isnan(u2(:))))
            fprintf("error!")
            return
        end
        
        u = u1 + N2*(u2 + N3*u3); %with task hierarchy
    %     u = u1 + u2 + u3; %without task hierarchy
%         u = u1 + N2*u2;

%          u = u1 + N2*(u3 + N4*u2); %with task hierarchy
%     u = u3 + N4*u2;

        a1 = u(1);
        omega_1 = u(2);
        a2 = u(3);
        omega_2 = u(4);

        % integrate the dynamics forward
        x1 = x1 + v1*cos(theta_1)*h;
        y1 = y1 + v1*sin(theta_1)*h;
        v1 = v1 + a1*h;
        theta_1 = theta_1 + omega_1*h;

        x2 = x2 + v2*cos(theta_2)*h;
        y2 = y2 + v2*sin(theta_2)*h;
        v2 = v2 + a2*h;
        theta_2 = theta_2 + omega_2*h;

        p1 = [x1; y1];
        p2 = [x2; y2];

        P1 = [P1, p1];
        P2 = [P2, p2];
        
        centroid = [centroid, (p1 + p2)/2];
        agents_dist_PI = [agents_dist_PI, norm(p1 - p2)];
        dist_from_goal_PI = [dist_from_goal_PI, norm((p1 + p2)/2 - pg)];
        desired_dist_bet_agents = [desired_dist_bet_agents, ell];
    end  
    plot(P1(1,:), P1(2,:), 'k', 'LineWidth', 1)
    plot(P2(1,:), P2(2,:), 'b', 'LineWidth', 1)
    plot (centroid(1, :), centroid(2, :), '--', 'LineWidth',1)
end

plot(pgx, pgy, '.', 'MarkerSize',50,'MarkerEdgeColor','g','MarkerFaceColor','g')
% plot(P1(1,:), P1(2,:), 'k', 'LineWidth', 1)
% plot(P2(1,:), P2(2,:), 'b', 'LineWidth', 1)
% plot (centroid(1, :), centroid(2, :), '--', 'LineWidth',1);
grid on;
legend('obstacle', 'goal position', 'agent1', 'agent2', 'centroid')

figure(3)
hold on;
plot(0:h:T-h, agents_dist_PI)
xlabel('t')
ylabel('distance between two agents')
plot(0:h:T-h, desired_dist_bet_agents);
title("Path Integral")
grid on

figure(4)
hold on;
plot(0:h:T-h, dist_from_goal_PI)
xlabel('t')
grid on
ylabel('distance between centroid of two agents and goal')
title("Path Integral")

save("for_plotting_PI.mat","agents_dist_PI", "dist_from_goal_PI", "desired_dist_bet_agents")