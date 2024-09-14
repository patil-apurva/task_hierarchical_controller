close all;
T = 10;
h = 0.1;

k1 = -0.2*[1, 0; 1, 0; 0, 1; 0, 1];
k3 = -1*[1; 1; 1; 1];

L = 1;
r = 1;

cx = -5;
cy = 0;

c = [cx; cy];

qx1_0 = -10;
qy1_0 = 0.5;
qx2_0 = -10;
qy2_0 = -0.5;

vx1_0 = 0;
vy1_0 = 0;
vx2_0 = 0;
vy2_0 = 0;

q1_0 = [qx1_0; qy1_0];
q2_0 = [qx2_0; qy2_0];
q_0 = [q1_0; q2_0];

v1_0 = [vx1_0; vy1_0];
v2_0 = [vx2_0; vy2_0];
v_0 = [v1_0; v2_0];

q_t = q_0;
v_t = v_0;
Q_t = [q_t];

for t=0:h:T-h
    
    qx1_t = q_t(1);
    qy1_t = q_t(2);
    qx2_t = q_t(3);
    qy2_t = q_t(4);
    
    q1_t = [qx1_t; qy1_t];
    q2_t = [qx2_t; qy2_t];
    
    J1 = [(qx1_t - cx)/norm(q1_t - c), (qy1_t - cy)/norm(q1_t - c), 0, 0; 0, 0, (qx2_t - cx)/norm(q2_t - c), (qy2_t - cy)/norm(q2_t - c)];
    
    J2 = 0.5*[1, 0, 1, 0; 0, 1, 0, 1];
    
    J3 = [(q1_t - q2_t).'; (q2_t - q1_t).'];
    
    N1 = eye(4);
    N2 = eye(4) - (J1.')*((pinv(J1)).');
    N3 = N2*(eye(4) - (J2.')*((pinv(J2)).'));
    
    a2 = [1; 0; 1; 0];
%     a2 = [0; 0; 0; 0];
   
    a3 = k3*(L^2 - ((q1_t - q2_t).')*(q1_t - q2_t));
%     a3 = [0; 0; 0; 0];
      
    a1 = k1*[r - norm(q1_t - c); r - norm(q2_t - c)];
%     a1 = [0; 0; 0; 0];
     
    if(norm(q1_t - c)>(r+1))
        a1(1:2) = 0;
    end
    
    if(norm(q2_t - c)>(r+1))
        a1(3:4) = 0;
    end
    
    q_tp1 = q_t + v_t*h;
    v_tp1 = v_t + N1*a1*h + N3*a3*h + N2*a2*h;
%     v_tp1 = v_t + a1*h + a2*h + a3*h;
     
    q_t = q_tp1;
    v_t = v_tp1;
    
    Q_t = [Q_t, q_t];
end

figure(5)
hold on;
p = nsidedpoly(1000, 'Center', [cx cy], 'Radius', r);
plot(p, 'FaceColor', 'r')
axis equal
plot(Q_t(1,:), Q_t(2,:))
plot(Q_t(3,:), Q_t(4,:))