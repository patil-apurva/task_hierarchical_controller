% close all;
clear
%time parameters
T = 1000;
h = 0.1;

% gains for acc_2
kp2 = 0.000731;
ki2 = 0;
kd2 = 0.144;

% kp2 = 0;
% ki2 = 0;
% kd2 = 0.1;

%errors for acc_2
error2_integral = [0; 0];
last_error2 = [0; 0];

%gains for acc_1
kp1 = 0.035;
ki1 = 0.000306;
kd1 = 1;

% kp1 = 0.01;
% ki1 = 0;
% kd1 = 1;

%errors for acc_1
error1_integral = 0;
last_error1 = 0;

%goal position
pg = [40; -3];

%obstacle radius and center
r = 8;
r_ext = 8;
cx = 20;
cy = -5;
c = [cx; cy];

%when to activate acc_1
threshold = 10;

%initial state
qx_0 = 0;
qy_0 = 0;
vx_0 = 0;
vy_0 = 0;

q_0 = [qx_0; qy_0];
v_0 = [vx_0; vy_0];

q_t = q_0;
v_t = v_0;
Q_t = [q_t];

for t=0:h:T-h
    
    qx_t = q_t(1);
    qy_t = q_t(2);
    
    vx_t = v_t(1);
    vy_t = v_t(2);
    
    %======================================================================
    %compute acc_2
    J2 = [1, 0; 0, 1];
    J2_dagger = J2.'/(J2*J2.');
    
    error2 = (pg - q_t);
    error2_integral = error2_integral + error2*h; %for I controller
    error2_derivative = (error2 - last_error2) / h; %for D controller
%     a2 = (J2.')*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
    a2 = (J2_dagger)*(kp2*error2 + ki2*error2_integral + kd2*error2_derivative);
    
    last_error2 = error2;
    
%     is_a2_ObsDir = dot(a2, (c - q_t) ) > 0;
    is_v_t_ObsDir = dot(v_t, (c - q_t) ) > 0;
    %======================================================================
    %compute acc_1
    error1 = (r_ext - norm(q_t - c));
    error1_integral = error1_integral + error1*h;
    error1_derivative = (error1 - last_error1) / h;
    
    last_error1 = error1;
    
    if(norm(q_t - c)<threshold && is_v_t_ObsDir)
        J1 = (q_t - c).' ./ norm(q_t - c);
        J1_dagger = J1.' / (J1 * J1.');
         
%         a1 = (J1.')*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative);
        a1 = J1_dagger*(kp1*error1 + ki1*error1_integral + kd1*error1_derivative);
        
%         N2 = eye(2) - (J1.')*((pinv(J1)).');
        N2 = eye(2) - J1_dagger*J1;
    else 
        a1 = zeros(2,1);
        N2 = eye(2);
    end
    
    ad = a1 + a2; %without task hierarchy
%     ad = a1 + N2*a2; %with task hierarchy
    
    q_tp1 = q_t + v_t*h;
    v_tp1 = v_t + ad*h;
     
    q_t = q_tp1;
    v_t = v_tp1;

%     dynamics = @(t,x)double_integrator(x, ad);
%     [t_tilde, X_data] = ode45(dynamics, [t, t + h], [q_t; v_t]);
    
    Q_t = [Q_t, q_t];
end

figure(3)
hold on;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal
plot(Q_t(1,:), Q_t(2,:))