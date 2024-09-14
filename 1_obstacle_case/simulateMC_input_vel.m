function S_tau_all = simulateMC_input_vel(eps_t_all_1, eps_t_all_2, xt_prime_1, xt_prime_2, f_xt_prime_1, f_xt_prime_2, t, h, T, b, s, cx, cy, r, threshold, pgx, pgy, eta, d, k1, N_J1_11, N_J1_12, N_J1_21, N_J1_22)

    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t
    
    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
%     safe_flag_tau = 1;
    
    for t_prime = t:h:T-h % this loop is to compute S(tau_i)
        
        if ((xt_prime_1 - cx)*(xt_prime_1 - cx) + (xt_prime_2 - cy)*(xt_prime_2 - cy) <= r*r)%if yes means the robot state is in the obstacle
            S_tau = S_tau + eta; 
%             safe_flag_tau = 0;
%             break; %end this tau
        end
        
        S_tau = S_tau + h*b*((xt_prime_1 - pgx) *(xt_prime_1 - pgx) + (xt_prime_2 - pgy)*(xt_prime_2 - pgy)); %add the state dependent running cost
        
        xt_prime_1 = xt_prime_1 + f_xt_prime_1*h + s*sqrt(h)*(N_J1_11*eps_t_prime_1 + N_J1_12*eps_t_prime_2); %move tau ahead
        xt_prime_2 = xt_prime_2 + f_xt_prime_2*h + s*sqrt(h)*(N_J1_21*eps_t_prime_1 + N_J1_22*eps_t_prime_2); %move tau ahead
        
%         is_eps_t_prime_ObsDir = eps_t_prime_1*(cx - xt_prime_1) + eps_t_prime_2*(cy - xt_prime_2) > 0;
        is_eps_t_prime_ObsDir = 1;
%         dot(v2_t, (c - xt) ) > 0
        
        if((xt_prime_1 - cx)*(xt_prime_1 - cx) + (xt_prime_2 - cy)*(xt_prime_2 - cy) <= threshold*threshold && is_eps_t_prime_ObsDir)
            J1_1 = (xt_prime_1 - cx)/sqrt((xt_prime_1 - cx)*(xt_prime_1 - cx) + (xt_prime_2 - cy)*(xt_prime_2 - cy));
            J1_2 = (xt_prime_2 - cy)/sqrt((xt_prime_1 - cx)*(xt_prime_1 - cx) + (xt_prime_2 - cy)*(xt_prime_2 - cy));
            
            J1_dagger_1 = J1_1/((J1_1)^2 + (J1_2)^2);
            J1_dagger_2 = J1_2/((J1_1)^2 + (J1_2)^2);  

            v1_t_prime_1 = J1_dagger_1*k1*(r - sqrt((xt_prime_1 - cx)*(xt_prime_1 - cx) + (xt_prime_2 - cy)*(xt_prime_2 - cy)));
            v1_t_prime_2 = J1_dagger_2*k1*(r - sqrt((xt_prime_1 - cx)*(xt_prime_1 - cx) + (xt_prime_2 - cy)*(xt_prime_2 - cy)));
            
            N_J1_11 = 1 - J1_dagger_1*J1_1; %N_J1_11 at new t_prime. 
            N_J1_12 = -J1_dagger_1*J1_2; %N_J1_12 at new t_prime. 
            N_J1_21 = -J1_dagger_2*J1_1; %N_J1_21 at new t_prime. 
            N_J1_22 = 1-J1_dagger_2*J1_2; %N_J1_22 at new t_prime. 
        else
            v1_t_prime_1 = 0;
            v1_t_prime_2 = 0;
            
            N_J1_11 = 1; %N_J1_11 at new t_prime. 
            N_J1_12 = 0; %N_J1_12 at new t_prime. 
            N_J1_21 = 0; %N_J1_21 at new t_prime.
            N_J1_22 = 1; %N_J1_22 at new t_prime. 
        end
        
        f_xt_prime_1 = v1_t_prime_1; %f_xt_prime_1 at new t_prime. Will be used in the next iteration 
        f_xt_prime_2 = v1_t_prime_2; %f_xt_prime_2 at new t_prime. Will be used in the next iteration 
        
        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
    
%     if(safe_flag_tau==1) %if tau has not collided 
        S_tau = S_tau + d*((xt_prime_1 - pgx) *(xt_prime_1 - pgx) + (xt_prime_2 - pgy)*(xt_prime_2 - pgy)); %add the terminal cost to S_tau
%     end
   
    S_tau_all = S_tau;