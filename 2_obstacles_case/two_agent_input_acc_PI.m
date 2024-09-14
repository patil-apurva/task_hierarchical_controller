clearvars
close all
run('parameters_input_acc_2agents.m')
seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator

tic
for traj_itr=1:traj_num
    traj_itr
    Q1 = []; %to store all positions of agent 1's trajectory
    Q2 = []; %to store all positions of agent 2's trajectory
    agents_dist_PI = []; %distance between 2 agents
    dist_from_goal_PI = []; %distance between centroid of two agents and goal
    desired_dist_bet_agents = [];
    centroid = []; %centorid of the two agents
    
    Q1 = [Q1, q1_0]; %stack the initial position of agent 1
    Q2 = [Q2, q2_0]; %stack the initial position of agent 2
    
    q1_t = q1_0; %start the agent 1's position from the given initial position
    q2_t = q2_0; %start the agent 2's position from the given initial position
    q_t = [q1_t; q2_t];
    
    v1_t = v1_0; %start the agent 1's velocity from the given initial velocity
    v2_t = v2_0; %start the agent 2's velocity from the given initial velocity
    v_t = [v1_t; v2_t];
    
    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
         t
         
        %compute a4_t
        J4 = [(q1_t - q2_t).', (q2_t - q1_t).'];
        J4_dagger = J4.'/(J4*J4.');

        sigma4 = 1/2*(q1_t - q2_t).'*(q1_t - q2_t);
        sigma4_desire = ell^2/2;

        error4 = sigma4_desire - sigma4;
        error4_integral = error4_integral + error4*h; %for I controller
        error4_derivative = (error4 - last_error4) / h; %for D controller

        a4_t = J4_dagger*(kp4*error4 + ki4*error4_integral + kd4*error4_derivative);
        
        a4_1_t = a4_t(1:2);
        a4_2_t = a4_t(3:4);

        last_error4 = error4;
        %==================================================================
        %compute a1_t (obs A) for agent 1
        v_low_1 = v_t(1:2);
        is_v_low_1_obsdir = dot(v_low_1, (cA-q1_t)) > 0;
        
        error1_1 = (r_extA - norm(q1_t - cA));
        error1_1_integral = error1_1_integral + error1_1*h;
        error1_1_derivative = (error1_1 - last_error1_1) / h;

        if(norm(q1_t - cA)<threshold && is_v_low_1_obsdir)
            J1_1 = (q1_t - cA).' ./ norm(q1_t - cA);
            J1_1_dagger = J1_1.' / (J1_1 * J1_1.');
            a1_1_t = J1_1_dagger*(kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative);
        else
            a1_1_t = zeros(2,1);
            J1_1 = zeros(1,2);
            J1_1_dagger = zeros(2,1);
        end

        last_error1_1 = error1_1;
        
        %==================================================================
        %compute a1_t (obs A) for agent 2
        v_low_2 = v_t(3:4);
        is_v_low_2_obsdir = dot(v_low_2, (cA-q2_t)) > 0;
        
        error1_2 = (r_extA - norm(q2_t - cA));
        error1_2_integral = error1_2_integral + error1_2*h;
        error1_2_derivative = (error1_2 - last_error1_2) / h;

        if(norm(q2_t - cA)<threshold && is_v_low_2_obsdir)
            J1_2 = (q2_t - cA).' ./ norm(q2_t - cA);
            J1_2_dagger = J1_2.' / (J1_2 * J1_2.');
            a1_2_t = J1_2_dagger*(kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative);
        else
            a1_2_t = zeros(2,1);
            J1_2 = zeros(1,2);
            J1_2_dagger = zeros(2,1);
        end

        last_error1_2 = error1_2;

        J1 = blkdiag(J1_1, J1_2);
        J1_dagger = blkdiag(J1_1_dagger, J1_2_dagger);

        N2 = eye(4) - J1_dagger*J1;

        a1_t = [a1_1_t; a1_2_t];
        %==================================================================
        %==================================================================
         %compute a2_t (obs B) for agent 1
        v_low_1 = v_t(1:2);
        is_v_low_1_obsdir = dot(v_low_1, (cB-q1_t)) > 0;
        
        error2_1 = (r_extB - norm(q1_t - cB));
        error2_1_integral = error2_1_integral + error2_1*h;
        error2_1_derivative = (error2_1 - last_error2_1) / h;

        if(norm(q1_t - cB)<threshold && is_v_low_1_obsdir)
            J2_1 = (q1_t - cB).' ./ norm(q1_t - cB);
            J2_1_dagger = J2_1.' / (J2_1 * J2_1.');
            a2_1_t = J2_1_dagger*(kp1*error2_1 + ki1*error2_1_integral + kd1*error2_1_derivative);
        else
            a2_1_t = zeros(2,1);
            J2_1 = zeros(1,2);
            J2_1_dagger = zeros(2,1);
        end

        last_error2_1 = error2_1;
        
        %==================================================================
        %compute a2_t (obs B) for agent 2
        v_low_2 = v_t(3:4);
        is_v_low_2_obsdir = dot(v_low_2, (cB-q2_t)) > 0;
        
        error2_2 = (r_extB - norm(q2_t - cB));
        error2_2_integral = error2_2_integral + error2_2*h;
        error2_2_derivative = (error2_2 - last_error2_2) / h;

        if(norm(q2_t - cB)<threshold && is_v_low_2_obsdir)
            J2_2 = (q2_t - cB).' ./ norm(q2_t - cB);
            J2_2_dagger = J2_2.' / (J2_2 * J2_2.');
            a2_2_t = J2_2_dagger*(kp1*error2_2 + ki1*error2_2_integral + kd1*error2_2_derivative);
        else
            a2_2_t = zeros(2,1);
            J2_2 = zeros(1,2);
            J2_2_dagger = zeros(2,1);
        end

        last_error2_2 = error2_2;

        J2 = blkdiag(J2_1, J2_2);
        J2_dagger = blkdiag(J2_1_dagger, J2_2_dagger);

        N3 = eye(4) - J2_dagger*J2;

        a2_t = [a2_1_t; a2_2_t];
        
         %compute a3
        J3 = 1/2*[1, 0, 1, 0; 0, 1, 0, 1];
        J3_dagger = J3.'/(J3*J3.');
        
        N4 = eye(4) - J3_dagger*J3;
         
        sigma3 = 1/2*(q1_t + q2_t);
        sigma3_desire = pg;

        error3 = sigma3_desire - sigma3;
        error3_integral = error3_integral + error3*h; %for I controller
        error3_derivative = (error3 - last_error3) / h; %for D controller

        a3_t = J3_dagger*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
% 
        last_error3 = error3;

%         eps_t_all_1 = randn(1, runs); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_2 = randn(1, runs); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_3 = randn(1, runs); %GPU array that stores eps_3(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_4 = randn(1, runs); %GPU array that stores eps_4(t) at the start of each sample path starting at time t and state xt
        
%         eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_3 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_3(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_4 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_4(t) at the start of each sample path starting at time t and state xt
%         
%         S_tau_all = arrayfun(@simulateMC_input_acc_2agents, eps_t_all_1, eps_t_all_2, eps_t_all_3, eps_t_all_4, ...
%             q1_t(1)*ones(1,runs), q1_t(2)*ones(1,runs), q2_t(1)*ones(1,runs), q2_t(2)*ones(1,runs), v1_t(1)*ones(1,runs), v1_t(2)*ones(1,runs), v2_t(1)*ones(1,runs), v2_t(2)*ones(1,runs), ...
%             a1_1_t(1)*ones(1,runs), a1_1_t(2)*ones(1,runs), a1_2_t(1)*ones(1,runs), a1_2_t(2)*ones(1,runs),...
%             a2_1_t(1)*ones(1,runs), a2_1_t(2)*ones(1,runs), a2_2_t(1)*ones(1,runs), a2_2_t(2)*ones(1,runs),...
%             a4_1_t(1)*ones(1,runs), a4_1_t(2)*ones(1,runs), a4_2_t(1)*ones(1,runs), a4_2_t(2)*ones(1,runs),...
%             t*ones(1,runs), h*ones(1,runs), T*ones(1,runs), b*ones(1,runs), e*ones(1,runs), s*ones(1,runs), ...
%             cxA*ones(1,runs), cyA*ones(1,runs), rA*ones(1,runs), r_extA*ones(1,runs), cxB*ones(1,runs), cyB*ones(1,runs), rB*ones(1,runs), r_extB*ones(1,runs),...
%             threshold*ones(1,runs), pgx*ones(1,runs), pgy*ones(1,runs), ell*ones(1,runs), eta1*ones(1,runs), eta2*ones(1,runs), ...
%             kp1*ones(1,runs), ki1*ones(1,runs), kd1*ones(1,runs), kp2*ones(1,runs), ki2*ones(1,runs), kd2*ones(1,runs), kp4*ones(1,runs), ki4*ones(1,runs), kd4*ones(1,runs),...
%             N2(1, 1)*ones(1,runs), N2(1, 2)*ones(1,runs), N2(1, 3)*ones(1,runs), N2(1, 4)*ones(1,runs), N2(2, 1)*ones(1,runs), N2(2, 2)*ones(1,runs), N2(2, 3)*ones(1,runs), N2(2, 4)*ones(1,runs), ...
%             N2(3, 1)*ones(1,runs), N2(3, 2)*ones(1,runs), N2(3, 3)*ones(1,runs), N2(3, 4)*ones(1,runs), N2(4, 1)*ones(1,runs), N2(4, 2)*ones(1,runs), N2(4, 3)*ones(1,runs), N2(4, 4)*ones(1,runs), ...
%             N3(1, 1)*ones(1,runs), N3(1, 2)*ones(1,runs), N3(1, 3)*ones(1,runs), N3(1, 4)*ones(1,runs), N3(2, 1)*ones(1,runs), N3(2, 2)*ones(1,runs), N3(2, 3)*ones(1,runs), N3(2, 4)*ones(1,runs), ...
%             N3(3, 1)*ones(1,runs), N3(3, 2)*ones(1,runs), N3(3, 3)*ones(1,runs), N3(3, 4)*ones(1,runs), N3(4, 1)*ones(1,runs), N3(4, 2)*ones(1,runs), N3(4, 3)*ones(1,runs), N3(4, 4)*ones(1,runs),...
%             N4(1, 1)*ones(1,runs), N4(1, 2)*ones(1,runs), N4(1, 3)*ones(1,runs), N4(1, 4)*ones(1,runs), N4(2, 1)*ones(1,runs), N4(2, 2)*ones(1,runs), N4(2, 3)*ones(1,runs), N4(2, 4)*ones(1,runs), ...
%             N4(3, 1)*ones(1,runs), N4(3, 2)*ones(1,runs), N4(3, 3)*ones(1,runs), N4(3, 4)*ones(1,runs), N4(4, 1)*ones(1,runs), N4(4, 2)*ones(1,runs), N4(4, 3)*ones(1,runs), N4(4, 4)*ones(1,runs));
%         
%         eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2; eps_t_all_3; eps_t_all_4]); %convert GPU array into a regular array
% 
%         denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
%         numer = eps_t_all_arr*(denom_i.'); %(size: (4 X 1))
%         denom = sum(denom_i); %scalar
%         
%         if(denom==0)
%             denom = 10^-3;
%         end
%         
%         N2N3 = N2*N3;
%         
%         cali_G = (N2N3.')/(N2N3*N2N3.');
% 
%         a3_t = (s/sqrt(h))*cali_G*N2N3*numer/denom; %the control input
%        
%         if(any(isnan(a3_t(:))))
%             fprintf("error!")
% %             continue;
%             return
%         end
        
        %==================================================================
        ad = a1_t + N2*a2_t + N2*N3*a3_t + N2*N3*N4*a4_t;
%         ad = a3_t;
        %move the trajectory forward
        q_t = q_t + v_t*h; 
        v_t = v_t + ad*h;
        
        q1_t = q_t(1:2);
        q2_t = q_t(3:4);
        
        v1_t = v_t(1:2);
        v2_t = v_t(3:4);
        
        Q1 = [Q1, q1_t]; %stack the new position
        Q2 = [Q2, q2_t]; %stack the new position
        centroid = [centroid, (q1_t + q2_t)/2];
        agents_dist_PI = [agents_dist_PI, norm(q1_t - q2_t)];
        dist_from_goal_PI = [dist_from_goal_PI, norm((q1_t + q2_t)/2 - pg)];
        desired_dist_bet_agents = [desired_dist_bet_agents, ell];
    end
    
    plot (Q1(1, :), Q1(2, :), 'b', 'LineWidth',1);
    plot (Q2(1, :), Q2(2, :), 'k', 'LineWidth',1);
    plot (centroid(1, :), centroid(2, :), '--', 'LineWidth',1);
    
    plot(pgx, pgy, '.', 'MarkerSize',100,'MarkerEdgeColor','g','MarkerFaceColor','g')
    legend('obstacle', 'agent 1', 'agent 2', 'centroid', 'goal position');
    
    figure(5)
    hold on;
    plot(t0:h:T-h, agents_dist_PI)
    xlabel('t')
    ylabel('distance between two agents')
    plot(t0:h:T-h, desired_dist_bet_agents);

    figure(6)
    hold on;
    plot(t0:h:T-h, dist_from_goal_PI)
    xlabel('t')
    ylabel('distance between centroid of two agents and goal')
end
toc
% save("for_plotting_PI.mat", "agents_dist_PI", "dist_from_goal_PI", "desired_dist_bet_agents")
save("for_plotting_PID.mat","agents_dist_PID", "dist_from_goal_PID")
