function S_tau_all = simulateMC_input_acc_2agents_only_PI(eps_t_all_1, eps_t_all_2, eps_t_all_3, eps_t_all_4, ...
    q1_t_prime_1, q1_t_prime_2, q2_t_prime_1, q2_t_prime_2, v1_t_prime_1, v1_t_prime_2, v2_t_prime_1, v2_t_prime_2, ...
    t, h, T, b, e, s, cx, cy, r, pgx, pgy, ell, eta1, eta2)

    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t
    eps_t_prime_3 = eps_t_all_3; %standard normal noise at t
    eps_t_prime_4 = eps_t_all_4; %standard normal noise at t
    
    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path

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
        
        %move the trajectory forward
        v1_t_prime_1 = v1_t_prime_1...
        + s*sqrt(h)*eps_t_prime_1; %move tau ahead
        
        v1_t_prime_2 = v1_t_prime_2...
        + s*sqrt(h)*eps_t_prime_2 ; %move tau ahead
    
        v2_t_prime_1 = v2_t_prime_1... 
        + s*sqrt(h)*eps_t_prime_3; %move tau ahead
        
        v2_t_prime_2 = v2_t_prime_2...
        + s*sqrt(h)*eps_t_prime_4; %move tau ahead
        
        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_3 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_4 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
      
    S_tau_all = S_tau;