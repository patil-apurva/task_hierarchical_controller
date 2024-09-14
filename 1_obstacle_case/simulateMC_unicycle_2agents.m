function S_tau_all = simulateMC_unicycle_2agents(eps_t_all_1, eps_t_all_2, eps_t_all_3, eps_t_all_4, ...
    x1_t_prime, y1_t_prime, v1_t_prime, theta_1_t_prime, x2_t_prime, y2_t_prime, v2_t_prime, theta_2_t_prime, ...
    a1_1_t_prime, omega1_1_t_prime, a1_2_t_prime, omega1_2_t_prime, a3_1_t_prime, omega3_1_t_prime, a3_2_t_prime, omega3_2_t_prime, ... 
    t, h, T, b, e, s, cx, cy, r, r_ext, threshold, pgx, pgy, ell, eta1, eta2, kp1, ki1, kd1, kp3, ki3, kd3,...
    N2_11, N2_12, N2_13, N2_14, N2_21, N2_22, N2_23, N2_24, N2_31, N2_32, N2_33, N2_34, N2_41, N2_42, N2_43, N2_44,...
    N3_11, N3_12, N3_13, N3_14, N3_21, N3_22, N3_23, N3_24, N3_31, N3_32, N3_33, N3_34, N3_41, N3_42, N3_43, N3_44)   

    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t
    eps_t_prime_3 = eps_t_all_3; %standard normal noise at t
    eps_t_prime_4 = eps_t_all_4; %standard normal noise at t
    
    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
    
    %parameters of PID controller
    error1_1_integral = 0;
    last_error1_1 = 0;
    
    error1_2_integral = 0;
    last_error1_2 = 0;
    
    error3_integral = 0;
    last_error3 = 0;
    
     for t_prime = t:h:T % this loop is to compute S(tau_i)
         
         if ((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2 <= r^2)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta1; 
         end
         
         if ((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2 <= r^2)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta2; 
         end
        
         S_tau = S_tau + h*b*(((x1_t_prime + x2_t_prime)/2 - pgx)^2 + ((y1_t_prime + y2_t_prime)/2 - pgy)^2)... %add the state dependent running cost
                 + h*e*0.5*(sqrt((x1_t_prime - x2_t_prime)^2 + (y1_t_prime - y2_t_prime)^2) - ell)^2;
             
%          S_tau = S_tau + h*b*((x1_t_prime - pgx)^2 + (y1_t_prime - pgy)^2)... %add the state dependent running cost
%                  + h*e*0.5*(sqrt((x1_t_prime - x2_t_prime)^2 + (y1_t_prime - y2_t_prime)^2) - ell)^2;
        
        %move the trajectory forward
        x1_t_prime = x1_t_prime + v1_t_prime*cos(theta_1_t_prime)*h;
        y1_t_prime = y1_t_prime + v1_t_prime*sin(theta_1_t_prime)*h;
        x2_t_prime = x2_t_prime + v2_t_prime*cos(theta_2_t_prime)*h;
        y2_t_prime = y2_t_prime + v2_t_prime*sin(theta_2_t_prime)*h;
        
        %N2*N3
        N2N3_11 = N2_11*N3_11 + N2_12*N3_21 + N2_13*N3_31 + N2_14*N3_41;
        N2N3_12 = N2_11*N3_12 + N2_12*N3_22 + N2_13*N3_32 + N2_14*N3_42;
        N2N3_13 = N2_11*N3_13 + N2_12*N3_23 + N2_13*N3_33 + N2_14*N3_43;
        N2N3_14 = N2_11*N3_14 + N2_12*N3_24 + N2_13*N3_34 + N2_14*N3_44;
        
        N2N3_21 = N2_21*N3_11 + N2_22*N3_21 + N2_23*N3_31 + N2_24*N3_41;
        N2N3_22 = N2_21*N3_12 + N2_22*N3_22 + N2_23*N3_32 + N2_24*N3_42;
        N2N3_23 = N2_21*N3_13 + N2_22*N3_23 + N2_23*N3_33 + N2_24*N3_43;
        N2N3_24 = N2_21*N3_14 + N2_22*N3_24 + N2_23*N3_34 + N2_24*N3_44;
        
        N2N3_31 = N2_31*N3_11 + N2_32*N3_21 + N2_33*N3_31 + N2_34*N3_41;
        N2N3_32 = N2_31*N3_12 + N2_32*N3_22 + N2_33*N3_32 + N2_34*N3_42;
        N2N3_33 = N2_31*N3_13 + N2_32*N3_23 + N2_33*N3_33 + N2_34*N3_43;
        N2N3_34 = N2_31*N3_14 + N2_32*N3_24 + N2_33*N3_34 + N2_34*N3_44;
        
        N2N3_41 = N2_41*N3_11 + N2_42*N3_21 + N2_43*N3_31 + N2_44*N3_41;
        N2N3_42 = N2_41*N3_12 + N2_42*N3_22 + N2_43*N3_32 + N2_44*N3_42;
        N2N3_43 = N2_41*N3_13 + N2_42*N3_23 + N2_43*N3_33 + N2_44*N3_43;
        N2N3_44 = N2_41*N3_14 + N2_42*N3_24 + N2_43*N3_34 + N2_44*N3_44;
        
        %move the trajectory forward
        v1_t_prime = v1_t_prime + a1_1_t_prime*h + N2N3_11*a3_1_t_prime*h + N2N3_12*omega3_1_t_prime*h + N2N3_13*a3_2_t_prime*h + N2N3_14*omega3_2_t_prime*h ...
                    + s*sqrt(h)*(N2_11*eps_t_prime_1 + N2_12*eps_t_prime_2 + N2_13*eps_t_prime_3 + N2_14*eps_t_prime_4); %move tau ahead
                
        theta_1_t_prime = theta_1_t_prime + omega1_1_t_prime*h + N2N3_21*a3_1_t_prime*h + N2N3_22*omega3_1_t_prime*h + N2N3_23*a3_2_t_prime*h + N2N3_24*omega3_2_t_prime*h ...
                         + s*sqrt(h)*(N2_21*eps_t_prime_1 + N2_22*eps_t_prime_2 + N2_23*eps_t_prime_3 + N2_24*eps_t_prime_4); %move tau ahead
                    
        v2_t_prime = v2_t_prime + a1_2_t_prime*h + N2N3_31*a3_1_t_prime*h + N2N3_32*omega3_1_t_prime*h + N2N3_33*a3_2_t_prime*h + N2N3_34*omega3_2_t_prime*h ...
                    + s*sqrt(h)*(N2_31*eps_t_prime_1 + N2_32*eps_t_prime_2 + N2_33*eps_t_prime_3 + N2_34*eps_t_prime_4); %move tau ahead
                
        theta_2_t_prime = theta_2_t_prime + omega1_2_t_prime*h + N2N3_41*a3_1_t_prime*h + N2N3_42*omega3_1_t_prime*h + N2N3_43*a3_2_t_prime*h + N2N3_44*omega3_2_t_prime*h ...
                        + s*sqrt(h)*(N2_41*eps_t_prime_1 + N2_42*eps_t_prime_2 + N2_43*eps_t_prime_3 + N2_44*eps_t_prime_4); %move tau ahead
                    
        %compute u3_t_prime
        K3_1 = (x1_t_prime - x2_t_prime)*cos(theta_1_t_prime) + (y1_t_prime - y2_t_prime)*sin(theta_1_t_prime);
        K3_2 = (x2_t_prime - x1_t_prime)*v1_t_prime*sin(theta_1_t_prime) + (y1_t_prime - y2_t_prime)*v1_t_prime*cos(theta_1_t_prime);
        K3_3 = (x2_t_prime - x1_t_prime)*cos(theta_2_t_prime) + (y2_t_prime - y1_t_prime)*sin(theta_2_t_prime);
        K3_4 = (x1_t_prime - x2_t_prime)*v2_t_prime*sin(theta_2_t_prime) + (y2_t_prime - y1_t_prime)*v2_t_prime*cos(theta_2_t_prime);
        
        K3_dagger_1 = K3_1/(K3_1^2 + K3_2^2 + K3_3^2 + K3_4^2);
        K3_dagger_2 = K3_2/(K3_1^2 + K3_2^2 + K3_3^2 + K3_4^2);
        K3_dagger_3 = K3_3/(K3_1^2 + K3_2^2 + K3_3^2 + K3_4^2);
        K3_dagger_4 = K3_4/(K3_1^2 + K3_2^2 + K3_3^2 + K3_4^2);
        
        sigma3 = 0.5*((x1_t_prime - x2_t_prime)^2 + (y1_t_prime - y2_t_prime)^2);
        sigma3_desire = ell^2/2;
        
        error3 = sigma3_desire - sigma3;
        error3_integral = error3_integral + error3*h; %for I controller
        error3_derivative = (error3 - last_error3) / h; %for D controller

        error3_tilde = v1_t_prime^2 + v2_t_prime^2 - 2*v1_t_prime*v2_t_prime*(cos(theta_1_t_prime)*cos(theta_2_t_prime) + sin(theta_1_t_prime)*sin(theta_2_t_prime));
        
        a3_1_t_prime = K3_dagger_1*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde);
        omega3_1_t_prime = K3_dagger_2*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde);
        a3_2_t_prime = K3_dagger_3*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde);
        omega3_2_t_prime = K3_dagger_4*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative - error3_tilde);
    
        last_error3 = error3;
        
        %compute u1_t_prime for agent 1
        uc_obs_angle_1 = atan((cy - y1_t_prime)/(cx - x1_t_prime)); %angle of the line connecting the unicycle 1 and the obstacle
        is_uc_in_ObsDir_1 = cos(uc_obs_angle_1 - theta_1_t_prime) > 0; %is the unicycle 1 in the obstacle direction

        if((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2 <= threshold*threshold && is_uc_in_ObsDir_1)
            K1_1_1 = ((x1_t_prime - cx)/sqrt((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2))*cos(theta_1_t_prime) + ((y1_t_prime - cy)/sqrt((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2))*sin(theta_1_t_prime);
            K1_1_2 = ((y1_t_prime - cy)*v1_t_prime/sqrt((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2))*cos(theta_1_t_prime) - ((x1_t_prime - cx)*v1_t_prime/sqrt((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2))*sin(theta_1_t_prime);
            
            K1_1_dagger_1 = K1_1_1/(K1_1_1^2 + K1_1_2^2);
            K1_1_dagger_2 = K1_1_2/(K1_1_1^2 + K1_1_2^2);
            
            error1_1 = r_ext - sqrt((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2);
            error1_1_derivative = (error1_1 - last_error1_1) / h; %for D controller
            error1_1_integral = error1_1_integral + error1_1*h; %for I controller
            
            error1_1_tilde = ((y1_t_prime - cy)^2*v1_t_prime^2*(cos(theta_1_t_prime))^2 + (x1_t_prime - cx)^2*v1_t_prime^2*(sin(theta_1_t_prime))^2)/((sqrt((x1_t_prime - cx)^2 + (y1_t_prime - cy)^2))^3);
            
            a1_1_t_prime = K1_1_dagger_1*(kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative - error1_1_tilde);
            omega1_1_t_prime = K1_1_dagger_2*(kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative - error1_1_tilde);
            
            last_error1_1 = error1_1;    
        else
            a1_1_t_prime = 0;
            omega1_1_t_prime = 0;
            
            K1_1_1 = 0;
            K1_1_2 = 0;
            
            K1_1_dagger_1 = 0;
            K1_1_dagger_2 = 0;
        end
        
        %compute u1_t_prime for agent 2
        uc_obs_angle_2 = atan((cy - y2_t_prime)/(cx - x2_t_prime)); %angle of the line connecting the unicycle 2 and the obstacle
        is_uc_in_ObsDir_2 = cos(uc_obs_angle_2 - theta_2_t_prime) > 0; %is the unicycle 2 in the obstacle direction

        if((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2 <= threshold*threshold && is_uc_in_ObsDir_2)
            K1_2_1 = ((x2_t_prime - cx)/sqrt((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2))*cos(theta_2_t_prime) + ((y2_t_prime - cy)/sqrt((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2))*sin(theta_2_t_prime);
            K1_2_2 = ((y2_t_prime - cy)*v2_t_prime/sqrt((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2))*cos(theta_2_t_prime) - ((x2_t_prime - cx)*v2_t_prime/sqrt((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2))*sin(theta_2_t_prime);
            
            K1_2_dagger_1 = K1_2_1/(K1_2_1^2 + K1_2_2^2);
            K1_2_dagger_2 = K1_2_2/(K1_2_1^2 + K1_2_2^2);
            
            error1_2 = r_ext - sqrt((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2);
            error1_2_derivative = (error1_2 - last_error1_2) / h; %for D controller
            error1_2_integral = error1_2_integral + error1_2*h; %for I controller
            
            error1_2_tilde = ((y2_t_prime - cy)^2*v2_t_prime^2*(cos(theta_2_t_prime))^2 + (x2_t_prime - cx)^2*v2_t_prime^2*(sin(theta_2_t_prime))^2)/((sqrt((x2_t_prime - cx)^2 + (y2_t_prime - cy)^2))^3);
            
            a1_2_t_prime = K1_2_dagger_1*(kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative - error1_2_tilde);
            omega1_2_t_prime = K1_2_dagger_2*(kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative - error1_2_tilde);
            
            last_error1_2 = error1_2;    
        else
            a1_2_t_prime = 0;
            omega1_2_t_prime = 0;
            
            K1_2_1 = 0;
            K1_2_2 = 0;
            
            K1_2_dagger_1 = 0;
            K1_2_dagger_2 = 0;
        end
        
        N2_11 = 1 - K1_1_1 * K1_1_dagger_1;
        N2_12 = -K1_1_2 * K1_1_dagger_1;
        N2_13 = 0;
        N2_14 = 0;
        
        N2_21 = -K1_1_1 * K1_1_dagger_2;
        N2_22 = 1 - K1_1_2 * K1_1_dagger_2;
        N2_23 = 0;
        N2_24 = 0;
        
        N2_31 = 0;
        N2_32 = 0;
        N2_33 = 1 - K1_2_1 * K1_2_dagger_1;
        N2_34 = -K1_2_2 * K1_2_dagger_1;
        
        N2_41 = 0;
        N2_42 = 0;
        N2_43 = -K1_2_1 * K1_2_dagger_2;
        N2_44 = 1 - K1_2_2 * K1_2_dagger_2;
        
        %==================================================================
        N3_11 = (v1_t_prime^2*v2_t_prime^2 + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2 + 2*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_12 = (v1_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v2_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);                                                          
        N3_13 = -(cos(theta_1_t_prime - theta_2_t_prime)*(2*v1_t_prime^2 + 2*v2_t_prime^2))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);                                   
        N3_14 = -(2*v2_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v1_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);


        N3_21 = (v1_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v2_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1); 
        N3_22 = (v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v2_t_prime^2 + 1)/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_23 = (2*v1_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v2_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_24 = -(4*v1_t_prime*v2_t_prime*cos(theta_1_t_prime - theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);


        N3_31 = -(cos(theta_1_t_prime - theta_2_t_prime)*(2*v1_t_prime^2 + 2*v2_t_prime^2))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_32 = (2*v1_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v2_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_33 = (v1_t_prime^2*v2_t_prime^2 + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 2*v1_t_prime^2 + v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_34 = -(v2_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v1_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);


        N3_41 = -(2*v2_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v1_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_42 = -(4*v1_t_prime*v2_t_prime*cos(theta_1_t_prime - theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_43 = -(v2_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v1_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
        N3_44 = (v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 1)/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);

        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_3 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_4 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
      
    S_tau_all = S_tau;
