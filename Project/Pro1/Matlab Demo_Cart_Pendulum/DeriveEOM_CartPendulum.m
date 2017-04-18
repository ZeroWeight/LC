
clear all
clc;

% coordinate paramters
syms  l    % length to the pendulum center of mass
syms  m    % mass of the pendulum
syms  I    % Inertia of the pendulum
syms  M    % inertia of the cart
syms b     % friction coefficient of the cart

% generalized coordinates
syms  theta % angle of wheel to vertical
syms  x     % angle of body to wheel

% environment parameters
syms  g     % gravity accleration

% force para
syms F

% generalized coordinates derivatives
syms theta_d 
syms theta_dd
syms x_d
syms x_dd


% coordinates by generalized coordinates
% ------------------------------------------------------------
Xcom = [    x;
            0;
            0;
            x - l*sin(theta);
            l*cos(theta);
            theta     ];
         
q = [       x;
            theta ];

q_d = [     x_d;
            theta_d ];

q_dd = [    x_dd;
            theta_dd ];

% Active Force
% ------------------------------------------------------------
Fa = [  F-x_d*b;
        -M*g;
        0;
        0;
       -m*g;
        0 
      ];
  
% -------------------------------------------------------------    
J = simplify(jacobian(Xcom,q));     
D = simplify(jacobian(J*q_d,q)*q_d);
MM = diag([M, M, 0, m, m, I]);


% linearized
% ------------------------------------------------------------ 
disp('--------------------- JMJ_R(educed) ------------------------')
JMJ = simplify(J.'*MM*J);
JMJ_R = simplify(subs(JMJ, cos(theta), 1));
JMJ_R = simplify(subs(JMJ_R, sin(theta), theta))


disp('--------------------- Matrix E_R(educed) -------------------')
E = simplify(J.'*(MM*D+Fa));
E_R = simplify(subs(E, cos(theta), 1));
E_R = simplify(subs(E_R, sin(theta), theta));
E_R = simplify(subs(E_R, theta^2, 0));
E_R = simplify(subs(E_R, theta_d^2, 0))


disp('-------------------------- EoM -----------------------------')
eom = simplify(JMJ_R\E_R)



%%
disp('---------------------- for Classical Control -----------------------')
syms X T S;
E = JMJ_R*q_dd - E_R;

ES = subs(E, x, X);
ES = subs(ES, x_d, X*S);
ES = subs(ES, x_dd, X*S^2);

ES = subs(ES, theta, T);
ES = subs(ES, theta_d, T*S);
ES = subs(ES, theta_dd, T*S^2);

ES = simplify(ES)

% 
% ------------------------------------------------------------- 
disp('you could derive the transfer function by hand here')
GS_2 = -jacobian(ES(2),T) / jacobian(ES(2),X)
GS = subs(ES(1), X, GS_2*T);
GS_1 = -simplify(jacobian(GS,F)/jacobian(GS,T))

qq=(M+m)*(I+m*l^2)-(m*l)^2;
GS_1H= (m*l/qq)*S^2 / (S^4+b*S^3*(I+m*l^2)/qq - (M+m)*m*g*l*S^2/qq - b*m*g*l*S/qq)


disp('---------------------- for Modern Control.m -----------------------')
A = [  [0,0,1,0];
       [0,0,0,1];
       jacobian(eom,q), jacobian(eom,q_d);
    ]

B= [    0;
        0;
        jacobian(eom,F)]
             
C= [1 0 0 0;
    0 1 0 0]
D = [0;0]