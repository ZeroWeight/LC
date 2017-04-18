clear all
clc;

%%
M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;


% classical control
S = tf('s');

GS_1 = 1/(((M + m)*((m*l^2 + I)*S^2 - g*l*m))/(l*m) - S^2*l*m + (b*((m*l^2 + I)*S^2 - g*l*m))/(S*l*m))
GS_2 = ((m*l^2 + I)*S^2 - g*l*m)/(S^2*l*m)

q=(M+m)*(I+m*l^2)-(m*l)^2;
GS= (m*l/q)*S^2 / (S^4+b*S^3*(I+m*l^2)/q - (M+m)*m*g*l*S^2/q - b*m*g*l*S/q)


% modern Control
A = [ 
 0,                                       0,                                      1, 0;
 0,                                       0,                                      0, 1;
 0,       (g*l^2*m^2)/(M*m*l^2 + I*m + I*M), -(b*m*l^2 + I*b)/(M*m*l^2 + I*m + I*M), 0;
 0, (l*m*(M*g + g*m))/(M*m*l^2 + I*m + I*M),         -(b*l*m)/(M*m*l^2 + I*m + I*M), 0;
];
 
B = [ 
                                 0
                                 0
 (m*l^2 + I)/(M*m*l^2 + I*m + I*M)
       (l*m)/(M*m*l^2 + I*m + I*M)
]; 

C = [
     1     0     0     0
     0     1     0     0
];

D = [
     0
     0
];

sys_ss = ss(A,B,C,D)

sys_tf = tf(sys_ss)

%
% 问题1： 为什么我们计算的GS_1和sys_tf中的一个不一致，请说明其原因。提示：结合经典和现代现代控制理论的特点。