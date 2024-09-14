clearvars
close all
run('parameters_input_acc_2agents.m')
seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator

for traj_itr=1:traj_num
    traj_itr
    Q1 = []; %to store all positions of agent 1's trajectory
    Q1 = [Q1, q1_0]; %stack the initial position of agent 1
    
    q1_t = q1_0; %start the agent 1's position from the given initial position
    q_t = [q1_t];
    
    v1_t = v1_0; %start the agent 1's velocity from the given initial velocity
    v_t = [v1_t];
    
    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
         t

%         eps_t_all_1 = randn(1, runs); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_2 = randn(1, runs); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_3 = randn(1, runs); %GPU array that stores eps_3(t) at the start of each sample path starting at time t and state xt
%         eps_t_all_4 = randn(1, runs); %GPU array that stores eps_4(t) at the start of each sample path starting at time t and state xt
        
        eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
        eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt
        
        S_tau_all = arrayfun(@simulateMC_input_acc_single_agent_only_PI, eps_t_all_1, eps_t_all_2, ...
            q1_t(1)*ones(1,runs), q1_t(2)*ones(1,runs), v1_t(1)*ones(1,runs), v1_t(2)*ones(1,runs),  ...
            t*ones(1,runs), h*ones(1,runs), T*ones(1,runs), b*ones(1,runs), s*ones(1,runs), cx*ones(1,runs), cy*ones(1,runs), r*ones(1,runs), pgx*ones(1,runs), pgy*ones(1,runs), eta1*ones(1,runs));
        
        eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2]); %convert GPU array into a regular array

        denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size: (4 X 1))
        denom = sum(denom_i); %scalar
        
%         if(denom==0)
%             denom = 10^-3;
%         end

        cali_G = eye(2);

        a2_t = (s/sqrt(h))*cali_G*numer/denom; %the control input
        
        if(any(isnan(a2_t(:))))
            fprintf("error!")
%             continue;
            return
        end
        
        ad = a2_t;
        %move the trajectory forward
        q_t = q_t + v_t*h; 
        v_t = v_t + ad*h;
        
        q1_t = q_t(1:2);
        
        v1_t = v_t(1:2);
        
        Q1 = [Q1, q1_t]; %stack the new position
        
    end
    
    plot (Q1(1, :), Q1(2, :), 'b', 'LineWidth',1);
end