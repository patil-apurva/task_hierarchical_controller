% close all;
clear
T = 100;
h = 0.1;

%PID gains
kp = 1;
ki = 0;
kd = 1;

%goal position 
pg = [40; 0];

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

%for PID
error_integral = 0; %for I controller
last_error = 0; %for D controller
error_x =  []; %error in the x direction 

for t=0:h:T-h
    
    qx_t = q_t(1);
    qy_t = q_t(2);
    
    vx_t = v_t(1);
    vy_t = v_t(2);
    
    J2 = [1, 0; 0, 1];
    
    error = (pg - q_t);
    
%     error_norm = [error_norm, norm(error)];
    error_x = [error_x, error(1)];
    
    error_integral = error_integral + error*h; %for I controller
    
    error_derivative = (error - last_error) / h; %for D controller
    
    a = (J2.')*(kp*error + ki*error_integral + kd*error_derivative); %PID controller
    
    last_error = error;
    
    q_tp1 = q_t + v_t*h;
    v_tp1 = v_t + a*h;
     
    q_t = q_tp1;
    v_t = v_tp1;
    
    Q_t = [Q_t, q_t];
end

figure(1)
hold on;
axis equal
plot(0:h:T-h, error_x)
plot(0:h:T-h, 0*ones(T/h))
xlabel('t')
ylabel('error')

figure(2)
plot(Q_t(1,:), Q_t(2,:))