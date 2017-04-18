clear
clc;

%%
m_B = 80;
I_B = 12.8;
l = 1.0;
m_W = 5;
I_W = 0.05;
R = 0.254;

g = 9.8015;
mu = 0.08;

% classical control
S = tf('s');

GS_1 = 1/(S*mu + S^2*(I_W + R^2*m_B + R^2*m_W) - (R^2*S^4*l^2*m_B^2)/((m_B*l^2 + I_B)*S^2 - g*l*m_B))

GS_2 = -(R*S^2*l*m_B)/((m_B*l^2 + I_B)*S^2 - g*l*m_B)


% modern Contorl
A = [
[ 0,                                                                                                             0,                                                                                            1, 0];
[ 0,                                                                                                             0,                                                                                            0, 1];
[ 0,                        -(R*g*l^2*m_B^2)/(I_B*I_W + I_B*R^2*m_B + I_B*R^2*m_W + I_W*l^2*m_B + R^2*l^2*m_B*m_W), -(m_B*mu*l^2 + I_B*mu)/(I_B*I_W + I_B*R^2*m_B + I_B*R^2*m_W + I_W*l^2*m_B + R^2*l^2*m_B*m_W), 0];
[ 0, (l*m_B*(I_W*g + R^2*g*m_B + R^2*g*m_W))/(I_B*I_W + I_B*R^2*m_B + I_B*R^2*m_W + I_W*l^2*m_B + R^2*l^2*m_B*m_W),           (R*l*m_B*mu)/(I_B*I_W + I_B*R^2*m_B + I_B*R^2*m_W + I_W*l^2*m_B + R^2*l^2*m_B*m_W), 0]];

B = [
                                                                                     0;
                                                                                     0;
 (m_B*l^2 + I_B)/(I_B*I_W + I_B*R^2*m_B + I_B*R^2*m_W + I_W*l^2*m_B + R^2*l^2*m_B*m_W);
      -(R*l*m_B)/(I_B*I_W + I_B*R^2*m_B + I_B*R^2*m_W + I_W*l^2*m_B + R^2*l^2*m_B*m_W)];
  
C = [1, 0, 0, 0; 0, 1, 0, 0];

D = [0;0];
 
sys_ss = ss(A, B, C, D)

sys_tf = tf(sys_ss)