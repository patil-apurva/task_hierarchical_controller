syms theta_1_t_prime theta_2_t_prime v1_t_prime v2_t_prime

K2 = 0.5*[cos(theta_1_t_prime) -v1_t_prime*sin(theta_1_t_prime) cos(theta_2_t_prime) -v2_t_prime*sin(theta_2_t_prime); ...
                sin(theta_1_t_prime) v1_t_prime*cos(theta_1_t_prime) sin(theta_2_t_prime) v2_t_prime*cos(theta_2_t_prime)];
            
K2_dagger = K2.'/(K2*K2.');

N3 = eye(4) - K2_dagger*K2;

N3 = simplify(N3)

% [ (v1^2*v2^2 + v1^2*cos(2*theta_1 - 2*theta_2) + v1^2 + 2*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     (v1*sin(2*theta_1 - 2*theta_2)*(v2^2 - 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),                                                          
%     -(cos(theta_1 - theta_2)*(2*v1^2 + 2*v2^2))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),                                   
%     -(2*v2*sin(theta_1 - theta_2)*(v1^2 + 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1)]
% 
% 
% [                                                           (v1*sin(2*theta_1 - 2*theta_2)*(v2^2 - 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1), 
%     (v2^2*cos(2*theta_1 - 2*theta_2) - cos(2*theta_1 - 2*theta_2) + 3*v2^2 + 1)/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     (2*v1*sin(theta_1 - theta_2)*(v2^2 + 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     -(4*v1*v2*cos(theta_1 - theta_2))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1)]
% 
% 
% [                                                          -(cos(theta_1 - theta_2)*(2*v1^2 + 2*v2^2))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     (2*v1*sin(theta_1 - theta_2)*(v2^2 + 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     (v1^2*v2^2 + v2^2*cos(2*theta_1 - 2*theta_2) + 2*v1^2 + v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     -(v2*sin(2*theta_1 - 2*theta_2)*(v1^2 - 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1)]
% 
% 
% [                                                            -(2*v2*sin(theta_1 - theta_2)*(v1^2 + 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     -(4*v1*v2*cos(theta_1 - theta_2))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     -(v2*sin(2*theta_1 - 2*theta_2)*(v1^2 - 1))/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1),
%     (v1^2*cos(2*theta_1 - 2*theta_2) - cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 1)/(v1^2*v2^2 - cos(2*theta_1 - 2*theta_2) + v1^2*cos(2*theta_1 - 2*theta_2) + v2^2*cos(2*theta_1 - 2*theta_2) + 3*v1^2 + 3*v2^2 - v1^2*v2^2*cos(2*theta_1 - 2*theta_2) + 1)]

(v1_t_prime^2*v2_t_prime^2 + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2 + 2*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
(v1_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v2_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
-(cos(theta_1_t_prime - theta_2_t_prime)*(2*v1_t_prime^2 + 2*v2_t_prime^2))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
-(2*v2_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v1_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);


(v1_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v2_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
(v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v2_t_prime^2 + 1)/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
(2*v1_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v2_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1); 
-(4*v1_t_prime*v2_t_prime*cos(theta_1_t_prime - theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);


-(cos(theta_1_t_prime - theta_2_t_prime)*(2*v1_t_prime^2 + 2*v2_t_prime^2))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
(2*v1_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v2_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
(v1_t_prime^2*v2_t_prime^2 + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 2*v1_t_prime^2 + v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
-(v2_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v1_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);


-(2*v2_t_prime*sin(theta_1_t_prime - theta_2_t_prime)*(v1_t_prime^2 + 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
-(4*v1_t_prime*v2_t_prime*cos(theta_1_t_prime - theta_2_t_prime))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
-(v2_t_prime*sin(2*theta_1_t_prime - 2*theta_2_t_prime)*(v1_t_prime^2 - 1))/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);
(v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 1)/(v1_t_prime^2*v2_t_prime^2 - cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v1_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 3*v1_t_prime^2 + 3*v2_t_prime^2 - v1_t_prime^2*v2_t_prime^2*cos(2*theta_1_t_prime - 2*theta_2_t_prime) + 1);