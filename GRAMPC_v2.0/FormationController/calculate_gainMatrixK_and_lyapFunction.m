%% Calculate gain matrix K from the LQR 

% satellite frequency w (rad/s) - can be modified (currently 500km altitude)
w=sqrt((6.6743e-11 * 5.9723e24)/(6378000 + 500000)^3);    

% create state-space model with matrices A, B, C, D from the Formation Dynamics
A=[0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 0 2*w;0 -w^2 0 0 0 0;0 0 3*w^2 -2*w 0 0];
B=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
C=eye(6);
D=zeros(6,3);
sys=ss(A,B,C,D);

% derive a gain matrix K from a LQR problem, using the state-space model
% defined above and the weighting matrices Q1 and R1
p1 = 1e-12;
Q1 = p1*(C'*C);
R1 = 0.7*eye(3);
[K] = lqr(sys,Q1,R1);

%% Calculate matrix P for the Lyapunov function V(x) = x^T*P*X

% rewrite dynamics with the control feedback matrix K:
% x_dot = Ax+Bu = Ax-BKx = (A-BK)x = A_lyap * x
temp = A - B * K;
A_lyap = temp';

% positive definite matrix Q2 with weight p2, so that the elements of P
% have an order of magnitude of approx. 1 to fit the values of state x
p2 = 1e-3;
Q2 = eye(6)*p2;  

% get matrix P as a solution to the Lyapunov equation A^T*P + P*A + Q = 0
% for A and Q set A = A_lyap and Q = Q2
P = lyap(A_lyap,Q2);                        
Pdot = temp' * P + P * temp; 
