close all;
clear;

%time parameters
T = 15;
h = 0.01;

%gains
k1 = 0.5; %gain for task 1
k2 = 0.5*eye(2); %gain for task 2
k3 = 0.5; %gain for task 3

%obstacle radius and center
r = 8;
cx = 40;
cy = 0;
c = [cx; cy];

%goal position
pg = [100; 0];

%when to activate v1
threshold = 10;

%initial states of two agents
qx1_0 = 0;
qy1_0 = 3;

qx2_0 = 0;
qy2_0 = -3;

%desired dist between the two agents
ell = 6;

q1_0 = [qx1_0; qy1_0];
q2_0 = [qx2_0; qy2_0];
q_0 = [q1_0; q2_0];

q1_t = q1_0;
q2_t = q2_0;
q_t = q_0;

Q1_t = [q1_t];
Q2_t = [q2_t];

for t=0:h:T-h
   
    %compute v3
    J3 = [(q1_t - q2_t).', (q2_t - q1_t).'];
    J3_dagger = J3.'/(J3*J3.');
    
    sigma3 = 1/2*(q1_t - q2_t).'*(q1_t - q2_t);
    sigma3_desire = ell^2/2;

    v3 = J3_dagger*k3*(sigma3_desire - sigma3);
    
    %compute v2
    J2 = 1/2*[1, 0, 1, 0; 0, 1, 0, 1];
    J2_dagger = J2.'/(J2*J2.');
    
    sigma2 = 1/2*(q1_t + q2_t);
    sigma2_desire = pg;
    
    v2 = J2_dagger*k2*(sigma2_desire - sigma2);
    N_J2 = eye(4) - J2_dagger*J2;
    
    low_vel = v2 + N_J2*v3;
    
    %compute v1 for agent 1
    v_low_1 = low_vel(1:2);
    is_v_low_1_obsdir = dot(v_low_1, (c-q1_t)) > 0;
    
    if(norm(q1_t - c)<threshold && is_v_low_1_obsdir)
        J1_1 = (q1_t - c).' ./ norm(q1_t - c);
        J1_1_dagger = J1_1.' / (J1_1 * J1_1.');
        v1_1 = J1_1_dagger*k1*(r - norm(q1_t - c));
    else
        v1_1 = zeros(2,1);
        J1_1 = zeros(1,2);
        J1_1_dagger = zeros(2,1);
    end
    %=====================================================
    %compute v1 for agent 2
    v_low_2 = low_vel(3:4);
    is_v_low_2_obsdir = dot(v_low_2, (c-q2_t)) > 0;
    
    if(norm(q2_t - c)<threshold && is_v_low_2_obsdir)
        J1_2 = (q2_t - c).' ./ norm(q2_t - c);
        J1_2_dagger = J1_2.' / (J1_2 * J1_2.');
        v1_2 = J1_2_dagger*k1*(r - norm(q2_t - c));
    else
        v1_2 = zeros(2,1);
        J1_2 = zeros(1,2);
        J1_2_dagger = zeros(2,1);
    end
    
    J1 = blkdiag(J1_1, J1_2);
    J1_dagger = blkdiag(J1_1_dagger, J1_2_dagger);
    
    N_J1 = eye(4) - J1_dagger*J1;
    
    v1 = [v1_1; v1_2];
    
    vd = v1 + N_J1*low_vel; %with task hierarchy
%     vd = v1 + v2 + v3; %without task hierarchy
    
    q_tp1 = q_t + vd*h;
    
    q_t = q_tp1;
    
    q1_t = q_t(1:2);
    q2_t = q_t(3:4);
    
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
% plot(x0(1), x0(2), 'p', 'MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','y')
plot(pg(1), pg(2), 'p', 'MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','g')