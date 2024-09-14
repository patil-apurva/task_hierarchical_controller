function S_tau_all = simulateMC_unicycle(eps_t_all_1, eps_t_all_2, xt_prime, yt_prime, vt_prime, theta_t_prime, ...
    a1_t_prime, omega1_t_prime, t, h, T, b, s, cx, cy, r, r_ext, threshold, pgx, pgy, eta, kp1, ki1, kd1, ...
    N2_11, N2_12, N2_21, N2_22)

    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t

    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
    
    %parameters of PID controller
    error1_integral = 0;
    last_error1 = 0;
    
    for t_prime = t:h:T % this loop is to compute S(tau_i)
        
        if ((xt_prime - cx)^2 + (yt_prime - cy)^2 <= r^2)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta; 
        end
        
        S_tau = S_tau + h*b*((xt_prime - pgx)^2 + (yt_prime - pgy)^2); %add the state dependent running cost
        
        xt_prime = xt_prime + vt_prime*cos(theta_t_prime)*h;
        yt_prime = yt_prime + vt_prime*sin(theta_t_prime)*h;
        vt_prime = vt_prime + a1_t_prime*h + s*sqrt(h)*(N2_11*eps_t_prime_1 + N2_12*eps_t_prime_2);
        theta_t_prime = theta_t_prime + omega1_t_prime*h + s*sqrt(h)*(N2_21*eps_t_prime_1 + N2_22*eps_t_prime_2);
        
        %compute u1_t_prime
        uc_obs_angle = atan((cy - yt_prime)/(cx - xt_prime)); %angle of the line connecting the unicycle and the obstacle
        is_uc_in_ObsDir = cos(uc_obs_angle - theta_t_prime) > 0; %is the unicycle in the obstacle direction
        
        if((xt_prime - cx)^2 + (yt_prime - cy)^2 <= threshold*threshold && is_uc_in_ObsDir)
            K1_1 = ((xt_prime - cx)/sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2))*cos(theta_t_prime) + ((yt_prime - cy)/sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2))*sin(theta_t_prime);
            K1_2 = ((yt_prime - cy)*vt_prime/sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2))*cos(theta_t_prime) - ((xt_prime - cx)*vt_prime/sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2))*sin(theta_t_prime);
            
            K1_dagger_1 = K1_1/(K1_1^2 + K1_2^2);
            K1_dagger_2 = K1_2/(K1_1^2 + K1_2^2);
            
            error1 = r_ext - sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2);
            error1_derivative = (error1 - last_error1) / h; %for D controller
            error1_integral = error1_integral + error1*h; %for I controller
            
            error1_tilde = ((yt_prime - cy)^2*vt_prime^2*(cos(theta_t_prime))^2 + (xt_prime - cx)^2*vt_prime^2*(sin(theta_t_prime))^2)/((sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2))^3);
            
            a1_t_prime = K1_dagger_1*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative - error1_tilde);
            omega1_t_prime = K1_dagger_2*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative - error1_tilde);
            
            last_error1 = error1; 
            
            NK2_11 = 1 - K1_1 * K1_dagger_1;
            NK2_12 = -K1_2 * K1_dagger_1;
            
            NK2_21 = -K1_1 * K1_dagger_2;
            NK2_22 = 1 - K1_2 * K1_dagger_2;
        
            J1_1 = (xt_prime - cx)/sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2);
            J1_2 = (yt_prime - cy)/sqrt((xt_prime - cx)^2 + (yt_prime - cy)^2);
            
            J1_dagger_1 = J1_1/(J1_1^2 + J1_2^2);
            J1_dagger_2 = J1_2/(J1_1^2 + J1_2^2);
            
            NJ2_11 = 1 - J1_1 * J1_dagger_1;
            NJ2_12 = -J1_2 * J1_dagger_1;
            
            NJ2_21 = -J1_1 * J1_dagger_2;
            NJ2_22 = 1 - J1_2 * J1_dagger_2;
            
        else
            a1_t_prime = 0;
            omega1_t_prime = 0;
           
            NK2_11 = 1;
            NK2_12 = 0;
            NK2_21 = 0;
            NK2_22 = 1;
            
            NJ2_11 = 1;
            NJ2_12 = 0;
            NJ2_21 = 0;
            NJ2_22 = 1;
        end
        
%         N2_11 = NJ2_11;
%         N2_12 = NJ2_12;
%         N2_21 = NJ2_21;
%         N2_22 = NJ2_22;
        
        N2_11 = NK2_11;
        N2_12 = NK2_12;
        N2_21 = NK2_21;
        N2_22 = NK2_22;

       eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
       eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        
    end
    
    S_tau_all = S_tau;