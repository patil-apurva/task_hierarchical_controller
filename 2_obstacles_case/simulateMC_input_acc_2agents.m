function S_tau_all = simulateMC_input_acc_2agents(eps_t_all_1, eps_t_all_2, eps_t_all_3, eps_t_all_4, ...
    q1_t_prime_1, q1_t_prime_2, q2_t_prime_1, q2_t_prime_2, v1_t_prime_1, v1_t_prime_2, v2_t_prime_1, v2_t_prime_2, ...
    a1_1_t_prime_1, a1_1_t_prime_2, a1_2_t_prime_1, a1_2_t_prime_2, ...
    a2_1_t_prime_1, a2_1_t_prime_2, a2_2_t_prime_1, a2_2_t_prime_2, ...
    a4_1_t_prime_1, a4_1_t_prime_2, a4_2_t_prime_1, a4_2_t_prime_2,...
    t, h, T, b, e, s, ...
    cxA, cyA, rA, r_extA, cxB, cyB, rB, r_extB, ...
    threshold, pgx, pgy, ell, eta1, eta2, ...
    kp1, ki1, kd1, kp2, ki2, kd2, kp4, ki4, kd4,...
    N2_11, N2_12, N2_13, N2_14, N2_21, N2_22, N2_23, N2_24, N2_31, N2_32, N2_33, N2_34, N2_41, N2_42, N2_43, N2_44,...
    N3_11, N3_12, N3_13, N3_14, N3_21, N3_22, N3_23, N3_24, N3_31, N3_32, N3_33, N3_34, N3_41, N3_42, N3_43, N3_44, ...
    N4_11, N4_12, N4_13, N4_14, N4_21, N4_22, N4_23, N4_24, N4_31, N4_32, N4_33, N4_34, N4_41, N4_42, N4_43, N4_44)


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
    
    error2_1_integral = 0;
    last_error2_1 = 0;
    
    error2_2_integral = 0;
    last_error2_2 = 0;
    
    error4_integral = 0;
    last_error4 = 0;

    for t_prime = t:h:T % this loop is to compute S(tau_i)
        
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
        
        %N2*N3*N4
        N2N3N4_11 = N2N3_11*N4_11 + N2N3_12*N4_21 + N2N3_13*N4_31 + N2N3_14*N4_41;
        N2N3N4_12 = N2N3_11*N4_12 + N2N3_12*N4_22 + N2N3_13*N4_32 + N2N3_14*N4_42;
        N2N3N4_13 = N2N3_11*N4_13 + N2N3_12*N4_23 + N2N3_13*N4_33 + N2N3_14*N4_43;
        N2N3N4_14 = N2N3_11*N4_14 + N2N3_12*N4_24 + N2N3_13*N4_34 + N2N3_14*N4_44;
        
        N2N3N4_21 = N2N3_21*N4_11 + N2N3_22*N4_21 + N2N3_23*N4_31 + N2N3_24*N4_41;
        N2N3N4_22 = N2N3_21*N4_12 + N2N3_22*N4_22 + N2N3_23*N4_32 + N2N3_24*N4_42;
        N2N3N4_23 = N2N3_21*N4_13 + N2N3_22*N4_23 + N2N3_23*N4_33 + N2N3_24*N4_43;
        N2N3N4_24 = N2N3_21*N4_14 + N2N3_22*N4_24 + N2N3_23*N4_34 + N2N3_24*N4_44;
        
        N2N3N4_31 = N2N3_31*N4_11 + N2N3_32*N4_21 + N2N3_33*N4_31 + N2N3_34*N4_41;
        N2N3N4_32 = N2N3_31*N4_12 + N2N3_32*N4_22 + N2N3_33*N4_32 + N2N3_34*N4_42;
        N2N3N4_33 = N2N3_31*N4_13 + N2N3_32*N4_23 + N2N3_33*N4_33 + N2N3_34*N4_43;
        N2N3N4_34 = N2N3_31*N4_14 + N2N3_32*N4_24 + N2N3_33*N4_34 + N2N3_34*N4_44;
        
        N2N3N4_41 = N2N3_41*N4_11 + N2N3_42*N4_21 + N2N3_43*N4_31 + N2N3_44*N4_41;
        N2N3N4_42 = N2N3_41*N4_12 + N2N3_42*N4_22 + N2N3_43*N4_32 + N2N3_44*N4_42;
        N2N3N4_43 = N2N3_41*N4_13 + N2N3_42*N4_23 + N2N3_43*N4_33 + N2N3_44*N4_43;
        N2N3N4_44 = N2N3_41*N4_14 + N2N3_42*N4_24 + N2N3_43*N4_34 + N2N3_44*N4_44;
        
        %move the trajectory forward
        v1_t_prime_1 = v1_t_prime_1 + a1_1_t_prime_1*h + N2_11*a2_1_t_prime_1*h + N2_12*a2_1_t_prime_2*h + N2_13*a2_2_t_prime_1*h + N2_14*a2_2_t_prime_2*h ...  
        + N2N3N4_11*a4_1_t_prime_1*h + N2N3N4_12*a4_1_t_prime_2*h + N2N3N4_13*a4_2_t_prime_1*h + N2N3N4_14*a4_2_t_prime_2*h...
        + s*sqrt(h)*(N2N3_11*eps_t_prime_1 + N2N3_12*eps_t_prime_2 + N2N3_13*eps_t_prime_3 + N2N3_14*eps_t_prime_4); %move tau ahead
        
        v1_t_prime_2 = v1_t_prime_2 + a1_1_t_prime_2*h + N2_21*a2_1_t_prime_1*h + N2_22*a2_1_t_prime_2*h + N2_23*a2_2_t_prime_1*h + N2_24*a2_2_t_prime_2*h...
        + N2N3N4_21*a4_1_t_prime_1*h + N2N3N4_22*a4_1_t_prime_2*h + N2N3N4_23*a4_2_t_prime_1*h + N2N3N4_24*a4_2_t_prime_2*h...
        + s*sqrt(h)*(N2N3_21*eps_t_prime_1 + N2N3_22*eps_t_prime_2 + N2N3_23*eps_t_prime_3 + N2N3_24*eps_t_prime_4); %move tau ahead
    
        v2_t_prime_1 = v2_t_prime_1 + a1_2_t_prime_1*h + N2_31*a2_1_t_prime_1*h + N2_32*a2_1_t_prime_2*h + N2_33*a2_2_t_prime_1*h + N2_34*a2_2_t_prime_2*h ...
        + N2N3N4_31*a4_1_t_prime_1*h + N2N3N4_32*a4_1_t_prime_2*h + N2N3N4_33*a4_2_t_prime_1*h + N2N3N4_34*a4_2_t_prime_2*h... 
        + s*sqrt(h)*(N2N3_31*eps_t_prime_1 + N2N3_32*eps_t_prime_2 + N2N3_33*eps_t_prime_3 + N2N3_34*eps_t_prime_4); %move tau ahead
        
        v2_t_prime_2 = v2_t_prime_2 + a1_2_t_prime_2*h + N2_41*a2_1_t_prime_1*h + N2_42*a2_1_t_prime_2*h + N2_43*a2_2_t_prime_1*h + N2_44*a2_2_t_prime_2*h...
        + N2N3N4_41*a4_1_t_prime_1*h + N2N3N4_42*a4_1_t_prime_2*h + N2N3N4_43*a4_2_t_prime_1*h + N2N3N4_44*a4_2_t_prime_2*h...
        + s*sqrt(h)*(N2N3_41*eps_t_prime_1 + N2N3_42*eps_t_prime_2 + N2N3_43*eps_t_prime_3 + N2N3_44*eps_t_prime_4); %move tau ahead
        
        %compute a4_t_prime
        J4_1 = q1_t_prime_1 - q2_t_prime_1;
        J4_2 = q1_t_prime_2 - q2_t_prime_2;
        J4_3 = -q1_t_prime_1 + q2_t_prime_1;
        J4_4 = -q1_t_prime_2 + q2_t_prime_2;
    
        J4_dagger_1 = J4_1/(J4_1^2 + J4_2^2 + J4_3^2 + J4_4^2);
        J4_dagger_2 = J4_2/(J4_1^2 + J4_2^2 + J4_3^2 + J4_4^2);
        J4_dagger_3 = J4_3/(J4_1^2 + J4_2^2 + J4_3^2 + J4_4^2);
        J4_dagger_4 = J4_4/(J4_1^2 + J4_2^2 + J4_3^2 + J4_4^2);
    
        sigma4 = 0.5*((q1_t_prime_1 - q2_t_prime_1)^2 + (q1_t_prime_2 - q2_t_prime_2)^2);
        sigma4_desire = ell^2/2;       
    
        error4 = sigma4_desire - sigma4;
        error4_integral = error4_integral + error4*h; %for I controller
        error4_derivative = (error4 - last_error4) / h; %for D controller
    
        a4_1_t_prime_1 = J4_dagger_1*(kp4*error4 + ki4*error4_integral + kd4*error4_derivative);
        a4_1_t_prime_2 = J4_dagger_2*(kp4*error4 + ki4*error4_integral + kd4*error4_derivative);
        a4_2_t_prime_1 = J4_dagger_3*(kp4*error4 + ki4*error4_integral + kd4*error4_derivative);
        a4_2_t_prime_2 = J4_dagger_4*(kp4*error4 + ki4*error4_integral + kd4*error4_derivative);
    
        last_error4 = error4;
        
        %=========================================================================================
        %compute a1_1_t_prime (Obs A)
        is_v_low_1_obsdir = v1_t_prime_1*(cxA - q1_t_prime_1) + v1_t_prime_2*(cyA - q1_t_prime_2) > 0;
        
        error1_1 = r_extA - sqrt((q1_t_prime_1 - cxA)^2 + (q1_t_prime_2 - cyA)^2);
        error1_1_integral = error1_1_integral + error1_1*h;
        error1_1_derivative = (error1_1 - last_error1_1) / h;
        
        if((q1_t_prime_1 - cxA)^2 + (q1_t_prime_2 - cyA)^2 <= threshold*threshold && is_v_low_1_obsdir)
%         if((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2 <= threshold*threshold)
            J1_1_1 = (q1_t_prime_1 - cxA)/sqrt((q1_t_prime_1 - cxA)^2 + (q1_t_prime_2 - cyA)^2);
            J1_1_2 = (q1_t_prime_2 - cyA)/sqrt((q1_t_prime_1 - cxA)^2 + (q1_t_prime_2 - cyA)^2);
            
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
        
        %compute a1_2_t_prime (Obs A)
        is_v_low_2_obsdir = v2_t_prime_1*(cxA - q2_t_prime_1) + v2_t_prime_2*(cyA - q2_t_prime_2) > 0;
        
        error1_2 = r_extA - sqrt((q2_t_prime_1 - cxA)^2 + (q2_t_prime_2 - cyA)^2);
        error1_2_integral = error1_2_integral + error1_2*h;
        error1_2_derivative = (error1_2 - last_error1_2) / h;
        
        if((q2_t_prime_1 - cxA)^2 + (q2_t_prime_2 - cyA)^2 <= threshold*threshold && is_v_low_2_obsdir)
%         if((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2 <= threshold*threshold)
            J1_2_1 = (q2_t_prime_1 - cxA)/sqrt((q2_t_prime_1 - cxA)^2 + (q2_t_prime_2 - cyA)^2);
            J1_2_2 = (q2_t_prime_2 - cyA)/sqrt((q2_t_prime_1 - cxA)^2 + (q2_t_prime_2 - cyA)^2);
            
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
        
        %==================================================================
        %compute a2_1_t_prime (Obs B)
        is_v_low_1_obsdir = v1_t_prime_1*(cxB - q1_t_prime_1) + v1_t_prime_2*(cyB - q1_t_prime_2) > 0;
        
        error2_1 = r_extB - sqrt((q1_t_prime_1 - cxB)^2 + (q1_t_prime_2 - cyB)^2);
        error2_1_integral = error2_1_integral + error2_1*h;
        error2_1_derivative = (error2_1 - last_error2_1) / h;
        
        if((q1_t_prime_1 - cxB)^2 + (q1_t_prime_2 - cyB)^2 <= threshold*threshold && is_v_low_1_obsdir)
%         if((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2 <= threshold*threshold)
            J2_1_1 = (q1_t_prime_1 - cxB)/sqrt((q1_t_prime_1 - cxB)^2 + (q1_t_prime_2 - cyB)^2);
            J2_1_2 = (q1_t_prime_2 - cyB)/sqrt((q1_t_prime_1 - cxB)^2 + (q1_t_prime_2 - cyB)^2);
            
            J2_1_dagger_1 = J2_1_1/(J2_1_1^2 + J2_1_2^2);
            J2_1_dagger_2 = J2_1_2/(J2_1_1^2 + J2_1_2^2);

            a2_1_t_prime_1 = J2_1_dagger_1*(kp2*error2_1 + ki2*error2_1_integral + kd2*error2_1_derivative);
            a2_1_t_prime_2 = J2_1_dagger_2*(kp2*error2_1 + ki2*error2_1_integral + kd2*error2_1_derivative);

        else
            a2_1_t_prime_1 = 0;
            a2_1_t_prime_2 = 0;
            
            J2_1_1 = 0;
            J2_1_2 = 0;
            
            J2_1_dagger_1 = 0;
            J2_1_dagger_2 = 0;
        end
        
        last_error2_1 = error2_1;
        
        %compute a2_2_t_prime (Obs B)
        is_v_low_2_obsdir = v2_t_prime_1*(cxB - q2_t_prime_1) + v2_t_prime_2*(cyB - q2_t_prime_2) > 0;
        
        error2_2 = r_extB - sqrt((q2_t_prime_1 - cxB)^2 + (q2_t_prime_2 - cyB)^2);
        error2_2_integral = error2_2_integral + error2_2*h;
        error2_2_derivative = (error2_2 - last_error2_2) / h;
        
        if((q2_t_prime_1 - cxB)^2 + (q2_t_prime_2 - cyB)^2 <= threshold*threshold && is_v_low_2_obsdir)
%         if((q2_t_prime_1 - cx)^2 + (q2_t_prime_2 - cy)^2 <= threshold*threshold)
            J2_2_1 = (q2_t_prime_1 - cxB)/sqrt((q2_t_prime_1 - cxB)^2 + (q2_t_prime_2 - cyB)^2);
            J2_2_2 = (q2_t_prime_2 - cyB)/sqrt((q2_t_prime_1 - cxB)^2 + (q2_t_prime_2 - cyB)^2);
            
            J2_2_dagger_1 = J2_2_1/(J2_2_1^2 + J2_2_2^2);
            J2_2_dagger_2 = J2_2_2/(J2_2_1^2 + J2_2_2^2);
            
            a2_2_t_prime_1 = J2_2_dagger_1*(kp2*error2_2 + ki2*error2_2_integral + kd2*error2_2_derivative);
            a2_2_t_prime_2 = J2_2_dagger_2*(kp2*error2_2 + ki2*error2_2_integral + kd2*error2_2_derivative);
        else 
            a2_2_t_prime_1 = 0;
            a2_2_t_prime_2 = 0;
            
            J2_2_1 = 0;
            J2_2_2 = 0;
            
            J2_2_dagger_1 = 0;
            J2_2_dagger_2 = 0;
        end
        last_error2_2 = error2_2;
        
        N3_11 = 1 - J2_1_1 * J2_1_dagger_1;
        N3_12 = -J2_1_2 * J2_1_dagger_1;
        N3_13 = 0;
        N3_14 = 0;
        
        N3_21 = -J2_1_1 * J2_1_dagger_2;
        N3_22 = 1 - J2_1_2 * J2_1_dagger_2;
        N3_23 = 0;
        N3_24 = 0;
        
        N3_31 = 0;
        N3_32 = 0;
        N3_33 = 1 - J2_2_1 * J2_2_dagger_1;
        N3_34 = -J2_2_2 * J2_2_dagger_1;
        
        N3_41 = 0;
        N3_42 = 0;
        N3_43 = -J2_2_1 * J2_2_dagger_2;
        N3_44 = 1 - J2_2_2 * J2_2_dagger_2;
        
        
        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_3 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_4 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
      
    S_tau_all = S_tau;