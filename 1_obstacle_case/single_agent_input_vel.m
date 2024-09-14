% close all;
clear;

%time parameters
T = 15;
h = 0.01;

%gains
k1 = 10; %gain for task 1
k2 = [1, 0; 0, 1]; %gain for task 2

%obstacle radius and center
r = 8;
cx = 20;
cy = -5;
c = [cx; cy];

%goal position
pg = [40; -3];

%when to activate v1
threshold = 10;

%initial state
qx_0 = 0;
qy_0 = 0;

q_0 = [qx_0; qy_0];

q_t = q_0;

Q_t = [q_t];

for t=0:h:T-h
    qx_t = q_t(1);
    qy_t = q_t(2);
    
    %compute v2
    J2 = [1, 0; 0, 1];
    J2_dagger = J2.'/(J2*J2.');
    
    v2 = J2_dagger*k2*(pg - q_t);
    
    is_v2_ObsDir = dot(v2, (c - q_t) ) > 0;
    
    %compute v1
    if(norm(q_t - c)<threshold && is_v2_ObsDir)
        J1 = (q_t - c).' ./ norm(q_t - c);
        J1_dagger = J1.' / (J1 * J1.');

        v1 = J1_dagger*k1*(r - norm(q_t - c));
        N_J1 = eye(2) - J1_dagger*J1;
    else
        v1 = zeros(2,1);
        N_J1 = eye(2);
    end
    
%     vd = v1 + N_J1*v2; %with task hierarchy
    vd = v1 + v2; %without task hierarchy
    
    q_tp1 = q_t + vd*h;
    
    q_t = q_tp1;
    Q_t = [Q_t, q_t];
    
end

figure(5)
hold on;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal
plot(Q_t(1,:), Q_t(2,:))