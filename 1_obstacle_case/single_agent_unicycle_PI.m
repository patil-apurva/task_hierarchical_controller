close all;
clearvars;
run('parameters_unicycle.m')

seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator

for traj_itr=1:10
    traj_itr
    Q_t = []; %to store all states of this trajectory
    Q_t = [Q_t, q_0]; %stack the initial state
    
    P_t = []; %to store all positions of this trajectory
    P_t = [P_t, p_0]; %stack the initial position
    
    x_t = x_0;
    y_t = y_0;
    v_t = v_0;
    theta_t = theta_0;
    
    q_t = [x_t; y_t; v_t; theta_t];
    p_t = [x_t; y_t];
    
    dist_from_goal_PI = []; %distance between the agent and goal

    for t=0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
%         t
        uc_obs_angle = atan((cy - y_t)/(cx - x_t)); %angle of the line connecting the unicycle and the obstacle
        is_uc_in_ObsDir = cos(uc_obs_angle - theta_t) > 0; %is the unicycle in the obstacle direction

        %%compute u1
        if(norm(p_t - c)<threshold && is_uc_in_ObsDir)
            K1 = [((x_t - cx)/norm(p_t-c))*cos(theta_t) + ((y_t - cy)/norm(p_t-c))*sin(theta_t) , ((y_t - cy)*v_t/norm(p_t-c))*cos(theta_t) - ((x_t - cx)*v_t/norm(p_t-c))*sin(theta_t)];
            K1_dagger = K1.' / (K1 * K1.');

            error1 = r_ext - norm(p_t - c); %for P controller
            error1_derivative = (error1 - last_error1) / h; %for D controller
            error1_integral = error1_integral + error1*h; %for I controller

            error1_tilde = ((y_t - cy)^2*v_t^2*(cos(theta_t))^2 + (x_t - cx)^2*v_t^2*(sin(theta_t))^2)/((norm(p_t-c))^3);

            u1_t = K1_dagger*((kp1*error1 + ki1*error1_integral + kd1*error1_derivative) - error1_tilde);

            last_error1 = error1; 

            N_K2 = eye(2) - K1_dagger*K1;

            J1 = [(x_t - cx)/norm(p_t - c), (y_t - cy)/norm(p_t - c), 0];
            J1_dagger = J1.' / (J1 * J1.');
            N_J2 = eye(3) - J1_dagger*J1;
        else
            u1_t = zeros(2,1);
            N_K2 = eye(2);
            N_J2 = eye(3);
        end
        
%         N2 = N_J2(1:2, 1:2);
        N2 = N_K2;
        
         %%compute u2
%         K2 = [cos(theta_t), -v_t*sin(theta_t); sin(theta_t), v_t*cos(theta_t)];
%         K2_dagger = K2.'/(K2*K2.');
% 
%         error2 = pg - p_t; %for P controller
%         error2_derivative = (error2 - last_error2) / h; %for D controller
%         error2_integral = error2_integral + error2*h; %for I controller
% 
%         u2_t = K2_dagger*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
%     %     u2_t = K2.'*(kp2*error2 + kd2*error2_derivative);
% 
%         last_error2 = error2; 

        eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
        eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt
        
        S_tau_all = arrayfun(@simulateMC_unicycle, eps_t_all_1, eps_t_all_2, x_t, y_t, v_t, theta_t, ...
    u1_t(1), u1_t(2), t, h, T, b, s, cx, cy, r, r_ext, threshold, pgx, pgy, eta, kp1, ki1, kd1, ...
    N2(1,1), N2(1,2), N2(2,1), N2(2,2));
        
        eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2]); %concatenate eps_t_all_arr_1 and eps_t_all_arr_2 in an array

        denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size:.rr. (2 X 1))
        denom = sum(denom_i); %scalar
        
        cali_G = (N2.')/(N2*N2.');
        
        u2_t = (s/sqrt(h))*cali_G*N2*numer/denom; %the control input
        
        if(any(isnan(u2_t(:))))
            fprintf("error!")
%             continue;
            return
        end
        
        ud = u1_t + N2*u2_t; %with task hierarchy
%         ud = u1_t + u2_t; %without task hierarchy
%         ud = u2_t;

        ad = ud(1);
        omega_d = ud(2);

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
        p_t = [x_t; y_t];

        Q_t = [Q_t, q_t];
        P_t = [P_t, p_t];
        
        dist_from_goal_PI = [dist_from_goal_PI, norm(p_t - pg)];
    end
    plot(P_t(1,:), P_t(2,:), 'b', 'LineWidth', 1)
end

plot(pgx, pgy, '.', 'MarkerSize',100,'MarkerEdgeColor','g','MarkerFaceColor','g')

grid on

% figure(6)
% hold on;
% plot(ERROR2(1,:))
% grid on
% title("error in x");
% xlabel("t");
% ylabel("error in x")
% 
% figure(7)
% hold on;
% plot(ERROR2(2,:))
% grid on
% title("error in y");
% xlabel("t");
% ylabel("error in y")
save("for_plotting_PI_single_agent.mat", "dist_from_goal_PI")