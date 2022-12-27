%% PART 1
clear; clc;
%% A
Vcc = 20; Vbe = 0.7; beta = 300; Vt = 0.025; Rs = 50; RL = 50e3;

Vc2 = 3/4*Vcc
Ve2 = 1/2*Vcc
Vb2 = Ve2 + Vbe

Vc1 = Ve2
Ve1 = 1/4*Vcc
Vb1 = Ve1 + Vbe

Rc = 2.4e3;
Ic2 = (Vcc-Vc2)/Rc
Ib2 = Ic2/beta
Ie2 = (1+beta)*Ib2
Ic1 = Ie2
Ib1 = Ic1/beta
Ie1 = (1+beta)*Ib1
I1 = Ie1/sqrt(beta)
I2 = I1-Ib2
I3 = I2-Ib1

Rb1 = (Vcc-Vb2)/I1
Rb2 = (Vb2-Vb1)/I2
Rb3 = Vb1/I3
Re = Ve1/Ie1

gm1 = Ic1/Vt;
gm2 = Ic2/Vt;
rpi1 = beta/gm1;
rpi2 = beta/gm2;
gm = gm1
rpi = rpi2

%Set all resistors to their standard values
Rc = 2.4e3; Rb1 = 75e3; Rb2 = 43e3; Rb3=51e3; Re=2.4e3;

RBB = RR(Rb2, Rb3)
R1 = 5e3-RR(RBB, rpi); R1 = 2e3;


Rce_sc = RR(Re,(rpi+RR(RBB, Rs+R1))/(1+beta))
Rcc2_sc = Rc+RL
Rcc1_sc = Rs+R1+RR(RBB, rpi)
Rcc1_oc = Rs+R1+RR(RBB, rpi+(1+beta)*Re)

syms CE
eqn1 = sqrt(1/(Rce_sc*CE)^2 - 2*(1/(Re * CE)^2)) == 1000*pi;
CE = double(solve(eqn1, CE));
CE = CE(2)
Av = -gm*RR(Rc,RL)*1*RR(RBB, rpi)/(RR(RBB, rpi)+Rs+R1)

%% B
CJE = 4.5e-12;
TF = 400e-12;
cpi = 2*CJE + TF*gm %#ok<*NOPTS> 

CJC = 3.5e-12;
VJC = 750e-3;
MJC = 330e-3;
Vcb = Vc2-Vb2;
cmu = CJC/(1+Vcb/VJC)^MJC

wlp1 = 1/(18.1*CE)
fl3db = wlp1/(2*pi)

whp1 = 1/((cpi+2*cmu)*(RR(RR(RBB, rpi), R1+Rs)))
whp2 = 1/((cpi+2*cmu)*(rpi/(1+beta)))
whp3 = 1/(cmu*RR(Rc,RL))
fh3db = whp1/(2*pi)

%% C
% x = [0 5e-3 10e-3 15e-3 20e-3 25e-3 30e-3 35e-3 40e-3 45e-3 50e-3 0.1 0.2 0.3 0.4 0.5]
% y = [591.6e-6, 60.65e-3 328.5e-3 702.4e-3 910.2e-3 1.206 1.071 670.4e-3 1.191 1.146 2.027 1.762 2.536 1.834 2.637 6.664]
% plot(x,y)

%% D
Rin = RR(RBB, rpi) + R1


function [f] = RR(R1, R2)
    f = (R1*R2)/(R1+R2);
end