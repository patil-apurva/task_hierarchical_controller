function S_tau_all = simulateMC_input_acc_single_agent_only_PI(eps_t_all_1, eps_t_all_2, ...
    q1_t_prime_1, q1_t_prime_2, v1_t_prime_1, v1_t_prime_2,...
    t, h, T, b, s, cx, cy, r, pgx, pgy, eta1)

    eps_t_prime_1 = eps_t_all_1; %standard normal noise at t
    eps_t_prime_2 = eps_t_all_2; %standard normal noise at t
    
    S_tau = 0; %the cost-to-go of the state dependent cost of a sample path

    for t_prime = t:h:T % this loop is to compute S(tau_i)
        
        if ((q1_t_prime_1 - cx)^2 + (q1_t_prime_2 - cy)^2 <= r^2)%if yes means the state is inside the obstacle
            S_tau = S_tau + eta1; 
        end
        
        S_tau = S_tau + h*b*((q1_t_prime_1 - pgx)^2 + (q1_t_prime_2  - pgy)^2); %add the state dependent running cost
        
        %move the trajectory forward
        q1_t_prime_1 = q1_t_prime_1 + v1_t_prime_1*h; 
        q1_t_prime_2 = q1_t_prime_2 + v1_t_prime_2*h;
        
        %move the trajectory forward
        v1_t_prime_1 = v1_t_prime_1...
        + s*sqrt(h)*eps_t_prime_1; %move tau ahead
        
        v1_t_prime_2 = v1_t_prime_2...
        + s*sqrt(h)*eps_t_prime_2 ; %move tau ahead
        
        eps_t_prime_1 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
        eps_t_prime_2 = randn; %standard normal noise at new t_prime. Will be used in the next iteration 
    end
      
    S_tau_all = S_tau;