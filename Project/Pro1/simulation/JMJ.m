clear
clc

syms l % length ti the body center of mass
syms m_B % mass of the body
syms I_B % Inertia of the body
syms R % radius of the wheel
syms m_W % mass of the wheel
syms I_W % Inertia of the wheel

% generalized coordinates
syms phi % angle of pendulum 
syms theta % angle of body to origin

% environment parameters
syms g % gravity accleration
syms mu % friction coefficient of the wheel rotation

% input parameter
syms tau % driving torque

% force parameter
syms F_W % static friction

% differential
syms theta_d theta_dd
syms phi_d phi_dd

% coordinates by generalized coordinates
% ------------------------------------------------------------
Xcom = [theta*R;
        0;
        theta;
        theta*R + l*sin(phi);
        l*cos(phi);
        phi ];

    
q = [theta;
     phi];

q_d = [theta_d;
       phi_d ];

q_dd = [theta_dd;
        phi_dd ];

% active force
% ------------------------------------------------------------
Fa = [ 0;
       -m_W * g;
       tau - mu * theta_d;
       0;
       -m_B * g;
       0];

% JMJ
% ------------------------------------------------------------
J = simplify(jacobian(Xcom, q));
D = simplify(jacobian(J * q_d, q) * q_d);
MM = diag([m_W, m_W, I_W, m_B, m_B, I_B]);
    
% linearized
% ------------------------------------------------------------ 
disp('--------------------- JMJ_R(educed) ------------------------')
JMJ = simplify( J.' * MM * J);
JMJ_R = simplify(subs(JMJ, [cos(phi), sin(phi)], [1, phi]))

disp('--------------------- Matrix E_R(educed) -------------------')
E = simplify(J.' * (MM*D + Fa));
E_R = simplify(subs(E, [cos(phi), sin(phi), phi^2, phi_d^2], [1, phi, 0, 0]))

disp('-------------------------- EoM -----------------------------')
eom = simplify(JMJ_R\E_R)

%%
disp('---------------------- for Classical Control -----------------------')
syms X T S;
E = JMJ_R*q_dd - E_R;

ES = subs(E, [phi, phi_d, phi_dd], [X, X*S, X*S^2]);
ES = subs(ES, [theta, theta_d, theta_dd], [T, T*S, T*S^2]);
ES = simplify(ES)

% 
% ------------------------------------------------------------- 
disp('you could derive the transfer function by hand here')
GS_2 = -jacobian(ES(2),T) / jacobian(ES(2),X)
GS = subs(ES(1), X, GS_2*T);
GS_1 = -simplify(jacobian(GS, tau)/jacobian(GS,T))

% qq=(M+m)*(I+m*l^2)-(m*l)^2;
% GS_1H= (m*l/qq)*S^2 / (S^4+b*S^3*(I+m*l^2)/qq - (M+m)*m*g*l*S^2/qq - b*m*g*l*S/qq)


disp('---------------------- for Modern Control.m -----------------------')
A = [  [0,0,1,0];
       [0,0,0,1];
       jacobian(eom, q), jacobian(eom, q_d);
    ]

B = [    0;
        0;
        jacobian(eom, tau)]
             
C = [1 0 0 0;
    0 1 0 0]
D = [0;0]