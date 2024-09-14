function S_tau_all = simulateMC_input_acc(eps_t_all_1, eps_t_all_2, qt_prime_1, qt_prime_2, vt_prime_1, vt_prime_2, a1_t_prime_1, a1_t_prime_2, t, h, T, b, s, cx, cy, r, r_ext, threshold, pgx, pgy, eta, d, kp1, ki1, kd1, N2_11, N2_12, N2_21, N2_22)
    
    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t
    
    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
%     safe_flag_tau = 1;
    
    error1_integral = 0;
    last_error1 = 0;

    for t_prime = t:h:T-h % this loop is to compute S(tau_i)
        
        if ((qt_prime_1 - cx)*(qt_prime_1 - cx) + (qt_prime_2 - cy)*(qt_prime_2 - cy) <= r*r)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta; 
%             safe_flag_tau = 0;
%             break; %end this tau
        end
        
        S_tau = S_tau + h*b*((qt_prime_1 - pgx) *(qt_prime_1 - pgx) + (qt_prime_2 - pgy)*(qt_prime_2 - pgy)); %add the state dependent running cost
        
        qt_prime_1 = qt_prime_1 + vt_prime_1*h; 
        qt_prime_2 = qt_prime_2 + vt_prime_2*h;
        
        vt_prime_1 = vt_prime_1 + a1_t_prime_1*h + s*sqrt(h)*(N2_11*eps_t_prime_1 + N2_12*eps_t_prime_2); %move tau ahead
        vt_prime_2 = vt_prime_2 + a1_t_prime_2*h + s*sqrt(h)*(N2_21*eps_t_prime_1 + N2_22*eps_t_prime_2); %move tau ahead
       
%         is_eps_t_prime_ObsDir = eps_t_prime_1*(cx - qt_prime_1) + eps_t_prime_2*(cy - qt_prime_2) > 0;
        is_eps_t_prime_ObsDir = 1;
        
        error1 = (r_ext - sqrt((qt_prime_1 - cx)*(qt_prime_1 - cx) + (qt_prime_2 - cy)*(qt_prime_2 - cy)));
        error1_integral = error1_integral + error1*h;
        error1_derivative = (error1 - last_error1) / h;
        
        last_error1 = error1;
        
        if((qt_prime_1 - cx)*(qt_prime_1 - cx) + (qt_prime_2 - cy)*(qt_prime_2 - cy) <= threshold*threshold && is_eps_t_prime_ObsDir)
            J1_1 = (qt_prime_1 - cx)/sqrt((qt_prime_1 - cx)*(qt_prime_1 - cx) + (qt_prime_2 - cy)*(qt_prime_2 - cy));
            J1_2 = (qt_prime_2 - cy)/sqrt((qt_prime_1 - cx)*(qt_prime_1 - cx) + (qt_prime_2 - cy)*(qt_prime_2 - cy));
            
            J1_sharp_1 = J1_1/((J1_1)^2 + (J1_2)^2);
            J1_sharp_2 = J1_2/((J1_1)^2 + (J1_2)^2);

            a1_t_prime_1 = J1_1*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative);
            a1_t_prime_2 = J1_2*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative);
            
            N2_11 = 1 - J1_1*J1_sharp_1;
            N2_12 = - J1_1*J1_sharp_2;
            N2_21 = - J1_2*J1_sharp_1;
            N2_22 = 1 - J1_2*J1_sharp_2;
        else
            a1_t_prime_1 = 0;
            a1_t_prime_2 = 0;
            
            N2_11 = 1;
            N2_12 = 0;
            N2_21 = 0;
            N2_22 = 1;
        end
        
        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
    
%     if(safe_flag_tau==1) %if tau has not collided 
        S_tau = S_tau + d*((qt_prime_1 - pgx) *(qt_prime_1 - pgx) + (qt_prime_2 - pgy)*(qt_prime_2 - pgy)); %add the terminal cost to S_tau
%     end
   
    S_tau_all = S_tau;