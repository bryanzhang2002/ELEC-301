clear; clc; close all;
f3dB = 10e3;
R = 10e3;
%% Part A
C = 1/(2*pi*f3dB*R);
A_M = 3 - sqrt(2); 

syms R1 R2
eqn1 = A_M == 1+ R2/R1;
eqn2 = R1+R2 == 10e3;
[R1, R2] = solve(eqn1,eqn2,R1,R2);
R1 = double(R1); R2 = double(R2);
H = A_M*(1/(R*C)^2)/(s^2+s*(3-A_M)/(R*C)+1/(R*C)^2);

s = tf('s');
G = -s/(R*C*s^2+3*s+1/(R*C));
figure();rlocus(G)
figure();pzmap(H); grid on

A_M = 3; syms R11 R22
eqn1 = A_M == 1+ R22/R1;
eqn2 = R11+R22 == 10e3;
[R11, R22] = solve(eqn1,eqn2,R11,R22);
R11 = double(R11); R22 = double(R22);

%% Part B
clear; clc;
R = 1e3; C = 1e-6;
f = 1/(2*pi*sqrt(6)*R*C);
f_doubled = 1/(2*pi*sqrt(6)*2*R*2*C);
f_halved= 1/(2*pi*sqrt(6)*1/2*R*1/2*C);

%% Part C
clear; clc
Vt = 25e-3;
Q1.Ic = 1.295e-3; Q1.Ib = 10.77e-6;
Rf = 100e3; Rs=5e3;

Q1.gm = Q1.Ic/Vt;
Q1.rpi = Vt/Q1.Ib;
Q1.hFE = Q1.rpi*Q1.Ic/Vt;

disp(Q1)

Q2.Ic = 2.184e-3; Q2.Ib = 15.34e-6; 

Q2.gm = Q2.Ic/Vt;
Q2.rpi = Vt/Q2.Ib;
Q2.hFE = Q2.rpi*Q2.Ic/Vt;

disp(Q2)

Rin = 706.2e-6/274.9e-9
Rout = 706.2e-6/11.29e-6

beta = -1/Rf;
Am = -127.4
Ai = Am * Rs
Af = Ai/(1+Ai*beta)
A = Af/Rs

f_L3dBf = 2.875 / (1+Ai*beta)
f_H3dBf = 91.48e3 * (1+Ai*beta)

Rif = Rin/(1+Ai*beta)
Rof = Rout / (1+Ai*beta)

Rin_f = 706.2e-6/281.7e-9

Af_1k = -0.1979 * Rs
Af_10k = -1.960 * Rs
Af_100k = -17.24 * Rs %#ok<*NOPTS> 
Af_1M = -77.80 * Rs
Af_10M = -119.8 * Rs

B_1k = -(Ai - Af_1k)/(Ai*Af_1k)
B_10k = -(Ai - Af_10k)/(Ai*Af_10k)
B_100k = -(Ai - Af_100k)/(Ai*Af_100k)
B_1M = -(Ai - Af_1M)/(Ai*Af_1M)
B_10M = -(Ai - Af_10M)/(Ai*Af_10M)

at_10k = 1 + -Ai * 0.1e-3
at_100k = 1 + -Ai * 10e-6
at_1M = 1 + -Ai * 1e-6

a_in_10k = Rin/26.3
a_in_100k = Rin/239
a_in_1M = Rin/1.31e3

a_out_10k = Rout/1.21
a_out_100k = Rout/8.34
a_out_1M = Rout/38.1

a_avg_10k = (a_in_10k+a_out_10k)/2
a_avg_100k = (a_in_100k+a_out_100k)/2
a_avg_1M = (a_in_1M+a_out_1M)/2

d_99k = sqrt((-127.49071 + 126.765187)/(-17.258379 + 17.198873))
d_101k = sqrt((-128.08 + 127.5)/(-17.298164 + 17.258379))

(d_99k + d_101k)/2