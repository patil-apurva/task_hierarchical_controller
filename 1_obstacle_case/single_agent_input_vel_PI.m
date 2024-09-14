clearvars
% close all
run('parameters_input_vel.m')

seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator.

for traj_itr=1:traj_num
    traj_itr
    X = []; %to store all positions of this trajectory
    X = [X, x0]; %stack the initial position

    xt = x0; %start the state from the given initial position
    f_xt = v1_0; %initial f_xt
    
    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
%          t
        eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
        eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt

        S_tau_all = arrayfun(@simulateMC_input_vel, eps_t_all_1, eps_t_all_2, xt(1), xt(2), f_xt(1), f_xt(2), t, h, T, b, s, cx, cy, r, threshold, pgx, pgy, eta, d, k1, N_J1(1,1), N_J1(1,2), N_J1(2,1), N_J1(2,2)); %an array that stores S(tau) of each sample path starting at time t and state xt

        eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2]); %concatenate eps_t_all_arr_1 and eps_t_all_arr_2 in an array

        denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size: (2 X 1))
        denom = sum(denom_i); %scalar

        cali_G = (N_J1.')/(N_J1*N_J1.');
        
        v2_t = (s/sqrt(h))*cali_G*N_J1*numer/denom; %the control input2
        
%         %compute v2_t
%         J2 = [1, 0; 0, 1];
%         J2_dagger = J2.'/(J2*J2.');
% 
%         v2_t = J2_dagger*k2*(pg - xt);
    
        
        if(any(isnan(v2_t(:))))
            fprintf("error!")
            return
        end
        
        is_v2_t_ObsDir = dot(v2_t, (c - xt) ) > 0; %check whether v2 is in the obstacle direction or not

%         if(norm(xt - c)<=r)
%             break;
%         end

%         eps = randn(2,1);
%         xt = xt + f_xt*h + N_J1*(v2_t*h + s*eps*sqrt(h)); %update the position with the control input ut=> x(t+h) = x(t) + f.h + g.u(t).h + sigma*dw
        xt = xt + f_xt*h + N_J1*v2_t*h; %move the trajectory forward
        X = [X,xt]; %stack the new position

%         compute v1_t and N_J1. Will be used in the next iteration.
        if(norm(xt - c)<threshold && is_v2_t_ObsDir)
            J1 = (xt - c).' ./ norm(xt - c);
            J1_dagger = J1.' / (J1 * J1.');

            v1_t = J1_dagger*k1*(r - norm(xt - c));
            N_J1 = eye(2) - J1_dagger*J1;
        else
            v1_t = zeros(2,1);
            N_J1 = eye(2);
        end
        
        f_xt = v1_t; %update f(x(t))

    end
    plot (X(1, :), X(2, :), 'b', 'LineWidth',1);
end
 
plot(x0(1), x0(2), 'p', 'MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','y')
plot(pgx, pgy, 'p', 'MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','g')

 