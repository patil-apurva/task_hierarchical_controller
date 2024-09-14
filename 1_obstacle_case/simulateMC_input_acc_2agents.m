function S_tau_all = simulateMC_input_acc_2agents(eps_t_all_1, eps_t_all_2, eps_t_all_3, eps_t_all_4, ...
    q1_t_prime_1, q1_t_prime_2, q2_t_prime_1, q2_t_prime_2, v1_t_prime_1, v1_t_prime_2, v2_t_prime_1, v2_t_prime_2, ...
    a1_1_t_prime_1, a1_1_t_prime_2, a1_2_t_prime_1, a1_2_t_prime_2, a3_1_t_prime_1, a3_1_t_prime_2, a3_2_t_prime_1, a3_2_t_prime_2,...
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
        
        if ((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2 <= r^2)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta1; 
        end
        
        if ((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2 <= r^2)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta2; 
        end
        
        S_tau = S_tau + h*b*(((q1_t_prime_1 + q2_t_prime_1)/2 - pgx)^2 + ((q1_t_prime_2 + q2_t_prime_2)/2 - pgy)^2)...; %add the state dependent running cost
                + h*e*0.5*(sqrt((q1_t_prime_1 - q2_t_prime_1)^2 + (q1_t_prime_2 - q2_t_prime_2)^2) - ell)^2;
        
        %move the trajectory forward
        q1_t_prime_1 = q1_t_prime_1 + v1_t_prime_1*h; 
        q1_t_prime_2 = q1_t_prime_2 + v1_t_prime_2*h;
        q2_t_prime_1 = q2_t_prime_1 + v2_t_prime_1*h; 
        q2_t_prime_2 = q2_t_prime_2 + v2_t_prime_2*h;
        
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
        v1_t_prime_1 = v1_t_prime_1 + a1_1_t_prime_1*h +  N2N3_11*a3_1_t_prime_1*h + N2N3_12*a3_1_t_prime_2*h + N2N3_13*a3_2_t_prime_1*h + N2N3_14*a3_2_t_prime_2*h...
        + s*sqrt(h)*(N2_11*eps_t_prime_1 + N2_12*eps_t_prime_2 + N2_13*eps_t_prime_3 + N2_14*eps_t_prime_4); %move tau ahead
        
        v1_t_prime_2 = v1_t_prime_2 + a1_1_t_prime_2*h + N2N3_21*a3_1_t_prime_1*h + N2N3_22*a3_1_t_prime_2*h + N2N3_23*a3_2_t_prime_1*h + N2N3_24*a3_2_t_prime_2*h...
        + s*sqrt(h)*(N2_21*eps_t_prime_1 + N2_22*eps_t_prime_2 + N2_23*eps_t_prime_3 + N2_24*eps_t_prime_4); %move tau ahead
    
        v2_t_prime_1 = v2_t_prime_1 + a1_2_t_prime_1*h + N2N3_31*a3_1_t_prime_1*h + N2N3_32*a3_1_t_prime_2*h + N2N3_33*a3_2_t_prime_1*h + N2N3_34*a3_2_t_prime_2*h... 
        + s*sqrt(h)*(N2_31*eps_t_prime_1 + N2_32*eps_t_prime_2 + N2_33*eps_t_prime_3 + N2_34*eps_t_prime_4); %move tau ahead
        
        v2_t_prime_2 = v2_t_prime_2 + a1_2_t_prime_2*h + N2N3_41*a3_1_t_prime_1*h + N2N3_42*a3_1_t_prime_2*h + N2N3_43*a3_2_t_prime_1*h + N2N3_44*a3_2_t_prime_2*h...
        + s*sqrt(h)*(N2_41*eps_t_prime_1 + N2_42*eps_t_prime_2 + N2_43*eps_t_prime_3 + N2_44*eps_t_prime_4); %move tau ahead
        
        %compute a3_t_prime
        J3_1 = q1_t_prime_1 - q2_t_prime_1;
        J3_2 = q1_t_prime_2 - q2_t_prime_2;
        J3_3 = -q1_t_prime_1 + q2_t_prime_1;
        J3_4 = -q1_t_prime_2 + q2_t_prime_2;
    
        J3_dagger_1 = J3_1/(J3_1^2 + J3_2^2 + J3_3^2 + J3_4^2);
        J3_dagger_2 = J3_2/(J3_1^2 + J3_2^2 + J3_3^2 + J3_4^2);
        J3_dagger_3 = J3_3/(J3_1^2 + J3_2^2 + J3_3^2 + J3_4^2);
        J3_dagger_4 = J3_4/(J3_1^2 + J3_2^2 + J3_3^2 + J3_4^2);
    
        sigma3 = 0.5*((q1_t_prime_1 - q2_t_prime_1)^2 + (q1_t_prime_2 - q2_t_prime_2)^2);
        sigma3_desire = ell^2/2;       
    
        error3 = sigma3_desire - sigma3;
        error3_integral = error3_integral + error3*h; %for I controller
        error3_derivative = (error3 - last_error3) / h; %for D controller
    
        a3_1_t_prime_1 = J3_dagger_1*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
        a3_1_t_prime_2 = J3_dagger_2*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
        a3_2_t_prime_1 = J3_dagger_3*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
        a3_2_t_prime_2 = J3_dagger_4*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
    
        last_error3 = error3;
        
        %compute a1_1_t_prime
        is_v_low_1_obsdir = v1_t_prime_1*(cx - q1_t_prime_1) + v1_t_prime_2*(cy - q1_t_prime_2) > 0;
        
        error1_1 = r_ext - sqrt((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2);
        error1_1_integral = error1_1_integral + error1_1*h;
        error1_1_derivative = (error1_1 - last_error1_1) / h;
        
        if((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2 <= threshold*threshold && is_v_low_1_obsdir)
%         if((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2 <= threshold*threshold)
            J1_1_1 = (q1_t_prime_1 - cx)/sqrt((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2);
            J1_1_2 = (q1_t_prime_2 - cy)/sqrt((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2);
            
            J1_1_dagger_1 = J1_1_1/(J1_1_1^2 + J1_1_2^2);
            J1_1_dagger_2 = J1_1_2/(J1_1_1^2 + J1_1_2^2);

            a1_1_t_prime_1 = J1_1_dagger_1*(kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative);
            a1_1_t_prime_2 = J1_1_dagger_2*(kp1*error1_1 + ki1*error1_1_integral + kd1*error1_1_derivative);

        else
            a1_1_t_prime_1 = 0;
            a1_1_t_prime_2 = 0;
            
            J1_1_1 = 0;
            J1_1_2 = 0;
            
            J1_1_dagger_1 = 0;
            J1_1_dagger_2 = 0;
        end
        
        last_error1_1 = error1_1;
        
        %compute a1_2_t_prime
        is_v_low_2_obsdir = v2_t_prime_1*(cx - q2_t_prime_1) + v2_t_prime_2*(cy - q2_t_prime_2) > 0;
        
        error1_2 = r_ext - sqrt((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2);
        error1_2_integral = error1_2_integral + error1_2*h;
        error1_2_derivative = (error1_2 - last_error1_2) / h;
        
        if((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2 <= threshold*threshold && is_v_low_2_obsdir)
%         if((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2 <= threshold*threshold)
            J1_2_1 = (q2_t_prime_1 - cx)/sqrt((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2);
            J1_2_2 = (q2_t_prime_2 - cy)/sqrt((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2);
            
            J1_2_dagger_1 = J1_2_1/(J1_2_1^2 + J1_2_2^2);
            J1_2_dagger_2 = J1_2_2/(J1_2_1^2 + J1_2_2^2);
            
            a1_2_t_prime_1 = J1_2_dagger_1*(kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative);
            a1_2_t_prime_2 = J1_2_dagger_2*(kp1*error1_2 + ki1*error1_2_integral + kd1*error1_2_derivative);
        else 
            a1_2_t_prime_1 = 0;
            a1_2_t_prime_2 = 0;
            
            J1_2_1 = 0;
            J1_2_2 = 0;
            
            J1_2_dagger_1 = 0;
            J1_2_dagger_2 = 0;
        end
        last_error1_2 = error1_2;
        
        N2_11 = 1 - J1_1_1 * J1_1_dagger_1;
        N2_12 = -J1_1_2 * J1_1_dagger_1;
        N2_13 = 0;
        N2_14 = 0;
        
        N2_21 = -J1_1_1 * J1_1_dagger_2;
        N2_22 = 1 - J1_1_2 * J1_1_dagger_2;
        N2_23 = 0;
        N2_24 = 0;
        
        N2_31 = 0;
        N2_32 = 0;
        N2_33 = 1 - J1_2_1 * J1_2_dagger_1;
        N2_34 = -J1_2_2 * J1_2_dagger_1;
        
        N2_41 = 0;
        N2_42 = 0;
        N2_43 = -J1_2_1 * J1_2_dagger_2;
        N2_44 = 1 - J1_2_2 * J1_2_dagger_2;
        
        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_3 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_4 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
      
    S_tau_all = S_tau;