clearvars
% close all
run('parameters_input_acc.m')

seed = 1234;
rng( seed ); % Reset the CPU random number generator.
gpurng( seed ); % Reset the GPU random number generator

all_path_len = 0;

for traj_itr=1:traj_num
    traj_itr
    Q = []; %to store all positions of this trajectory
    Q = [Q, q0]; %stack the initial position

    qt = q0; %start the state from the given initial position
    vt = v0; %start the state from the given initial position
    a1_t = a1_0;
    
    path_len = 0; 
    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
%          t
        eps_t_all_1 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_1(t) at the start of each sample path starting at time t and state xt
        eps_t_all_2 = randn(1, runs, 'gpuArray'); %GPU array that stores eps_2(t) at the start of each sample path starting at time t and state xt

        S_tau_all = arrayfun(@simulateMC_input_acc, eps_t_all_1, eps_t_all_2, qt(1)*ones(1, runs), qt(2)*ones(1, runs), vt(1)*ones(1, runs), vt(2)*ones(1, runs), a1_t(1)*ones(1, runs), a1_t(2)*ones(1, runs), t*ones(1, runs), h*ones(1, runs), T*ones(1, runs), b*ones(1, runs), s*ones(1, runs), cx*ones(1, runs), cy*ones(1, runs), r*ones(1, runs), r_ext*ones(1, runs), threshold*ones(1, runs), pgx*ones(1, runs), pgy*ones(1, runs), eta*ones(1, runs), d*ones(1, runs), kp1*ones(1, runs), ki1*ones(1, runs), kd1*ones(1, runs), N2(1,1)*ones(1, runs), N2(1,2)*ones(1, runs), N2(2,1)*ones(1, runs), N2(2,2)*ones(1, runs)); %an array that stores S(tau) of each sample path starting at time t and state xt

        eps_t_all_arr = gather([eps_t_all_1; eps_t_all_2]); %concatenate eps_t_all_arr_1 and eps_t_all_arr_2 in an array

        denom_i = exp(-S_tau_all/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size: (2 X 1))
        denom = sum(denom_i); %scalar
        
%         if(denom==0)
%             denom = 10^-3;
%         end

        cali_G = (N2.')/(N2*N2.');

        a2_t = (s/sqrt(h))*cali_G*N2*numer/denom; %the control input

%         %compute acc_2
%         J2 = [1, 0; 0, 1];
% 
%         error2 = (pg - qt);
%         error2_integral = error2_integral + error2*h; %for I controller
%         error2_derivative = (error2 - last_error2) / h; %for D controller
%         a2_t = (J2.')*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
% 
%         last_error2 = error2;
        
        if(any(isnan(a2_t(:))))
            fprintf("error!")
%             continue;
            return
        end
        
        is_vt_ObsDir = dot(vt, (c - qt) ) > 0;
        
%         eps = randn(2,1);
%         xt = xt + f_xt*h + N_J1*(v2_t*h + s*eps*sqrt(h)); %update the position with the control input ut=> x(t+h) = x(t) + f.h + g.u(t).h + sigma*dw
        prev_qt = qt;
        
        qt = qt + vt*h; 
        vt = vt + a1_t*h + N2*a2_t*h;
        Q = [Q,qt]; %stack the new position
        
        path_len = path_len + norm(qt - prev_qt);
        
%         if(norm(xt - c)<=r)
%             break;
%         end
        
        % compute a1_t and N2. Will be used in the next iteration 
        error1 = (r_ext - norm(qt - c));
        error1_integral = error1_integral + error1*h;
        error1_derivative = (error1 - last_error1) / h;
    
        last_error1 = error1;
%     
        if(norm(qt - c)<threshold && is_vt_ObsDir)
            J1 = (qt - c).' ./ norm(qt - c);
            a1_t = (J1.')*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative);
            N2 = eye(2) - (J1.')*((pinv(J1)).');
        else 
            a1_t = zeros(2,1);
            N2 = eye(2);
        end
    end
    plot (Q(1, :), Q(2, :), 'b', 'LineWidth',1);
    all_path_len = all_path_len + path_len;
end

path_len_avg = all_path_len/traj_num

plot(q0(1), q0(2), 'p', 'MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','y')
plot(pgx, pgy, 'p', 'MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','g')