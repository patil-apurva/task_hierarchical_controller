% close all;
clear;

%time parameters
T = 150;
h = 0.01;

% gains for acc_3 
kp3 = 0.00731;
ki3 = 0;
kd3 = 0.144;

%errors for acc_3
error3_integral = 0;
last_error3 = 0;

% gains for acc_2
kp2 = 0.000731;
ki2 = 0;
kd2 = 0.144;

%errors for acc_2
error2_integral = zeros(2,1);
last_error2 = zeros(2,1);

%gains for acc_1
kp1 = 0.01;
ki1 = 0;
kd1 = 1;

%errors for acc_1 - agent 1
error1_1_integral = 0;
last_error1_1 = 0;

%errors for acc_1 - agent 2
error1_2_integral = 0;
last_error1_2 = 0;

r = 1;
r_ext = 2;
% r_ext = r;
cx = 3;
cy = 0;
c = [cx; cy];

%goal position
pg = [6; 0];

%when to activate a1
threshold = 2;

%initial states of two agents
qx1_0 = 0;
qy1_0 = 0.5;
vx1_0 = 0;
vy1_0 = 0;

qx2_0 = 0;
qy2_0 = -0.5;
vx2_0 = 0;
vy2_0 = 0;

%desired dist between the two agents
ell = 1;

q1_0 = [qx1_0; qy1_0];
q2_0 = [qx2_0; qy2_0];
q_0 = [q1_0; q2_0];

v1_0 = [vx1_0; vy1_0];
v2_0 = [vx2_0; vy2_0];
v_0 = [v1_0; v2_0];

q1_t = q1_0;
q2_t = q2_0;
q_t = q_0;

v1_t = v1_0;
v2_t = v2_0;
v_t = v_0;

Q1_t = [q1_t];
Q2_t = [q2_t];

for t=0:h:T-h
   t
    %compute a3
    J3 = [(q1_t - q2_t).', (q2_t - q1_t).'];
    J3_dagger = J3.'/(J3*J3.');
    
    sigma3 = 1/2*(q1_t - q2_t).'*(q1_t - q2_t);
    sigma3_desire = ell^2/2;
    
    error3 = sigma3_desire - sigma3;
    error3_integral = error3_integral + error3*h; %for I controller
    error3_derivative = (error3 - last_error3) / h; %for D controller

    a3 = J3_dagger*(kp3*error3 + ki3*error3_integral + kd3*error3_derivative);
    
    last_error3 = error3;
    
    %compute a2
    J2 = 1/2*[1, 0, 1, 0; 0, 1, 0, 1];
    J2_dagger = J2.'/(J2*J2.');
    
    sigma2 = 1/2*(q1_t + q2_t);
    sigma2_desire = pg;
    
    error2 = sigma2_desire - sigma2;
    error2_integral = error2_integral + error2*h; %for I controller
    error2_derivative = (error2 - last_error2) / h; %for D controller
    
    a2 = J2_dagger*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
    
    last_error2 = error2;
    
    N3 = eye(4) - J2_dagger*J2;
    
    %compute a1 for agent 1
    
    obs_to_robot_pos = c - q1_t;
    dist_from_obs = norm(c - q1_t);
    obs_to_robot_pos_vec = obs_to_robot_pos/dist_from_obs;
    
    error1_1 = r_ext - dist_from_obs; 
    if(norm(q1_t - c)<threshold)
        J1_1 = (q1_t - c).' ./ norm(q1_t - c);
        J1_1_dagger = J1_1.' / (J1_1 * J1_1.');
        a1_1 = -kp1*error1_1*obs_to_robot_pos_vec;
    else
        a1_1 = zeros(2,1);
        J1_1 = zeros(1,2);
        J1_1_dagger = zeros(2,1);
    end
    
    %=====================================================
    %compute a1 for agent 2
    obs_to_robot_pos = c - q2_t;
    dist_from_obs = norm(c - q2_t);
    obs_to_robot_pos_vec = obs_to_robot_pos/dist_from_obs;
    
    error1_2 = r_ext - dist_from_obs; 
    
    if(norm(q2_t - c)<threshold)
        J1_2 = (q2_t - c).' ./ norm(q2_t - c);
        J1_2_dagger = J1_2.' / (J1_2 * J1_2.');
        a1_2 = -kp1*error1_2*obs_to_robot_pos_vec;
    else
        a1_2 = zeros(2,1);
        J1_2 = zeros(1,2);
        J1_2_dagger = zeros(2,1);
    end
      
    J1 = blkdiag(J1_1, J1_2);
    J1_dagger = blkdiag(J1_1_dagger, J1_2_dagger);
    
    N2 = eye(4) - J1_dagger*J1;
    
    a1 = [a1_1; a1_2];
    
    ad = a1 + N2*(a2 + N3*a3); %with task hierarchy
%     ad = a1 + a2 + a3; %without task hierarchy
    
    q_tp1 = q_t + v_t*h;
    v_tp1 = v_t + ad*h;
    
    q_t = q_tp1;
    v_t = v_tp1;
    
    q1_t = q_t(1:2);
    q2_t = q_t(3:4);
    
    v1_t = v_t(1:2);
    v2_t = v_t(3:4);
    
    Q1_t = [Q1_t, q1_t];
    Q2_t = [Q2_t, q2_t];
    
end

figure(5)
hold on;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal
plot(Q1_t(1,:), Q1_t(2,:))
plot(Q2_t(1,:), Q2_t(2,:))