function dx = double_integrator(x,u)
A = [0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0];
B = [0, 0; 0, 0; 1, 0; 0, 1];
dx = A*x+B*u;