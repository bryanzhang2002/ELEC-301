%% PART 2
clc; clear
%% A
Vcc = 12; Vbe = 0.7; beta = 300; alpha = beta/(beta+1); Vt = 0.025;
Rin = 50;
Ve1 = Vcc/3
Vb1 = Ve1 + Vbe
Vc2 = Vcc
Vc1 = 2/3*Vcc
Vb2 = Vc1;
Ve2 = Vb2 - Vbe

Ie1 = Vt/Rin
Ic1 = Ie1;
Ib1 = Ic1/beta
I1 = 0.1 * Ie1
I3 = I1 - Ib1

Rb1 = (Vcc-Vb1)/I1
Rb2 = Vb1/I3


syms RE2
% eqn1 = RR(RE2, (Vt*beta)/((alpha*Ve2/RE2)*(1+beta))+(Vcc-Vb2)/((Ic1+(Ve2/((beta+1)*RE2)))*(1+beta))) == 50;
Rc1 = (Vcc-Vc1)/(Ic1 + 1/(1+beta)*Ve2/RE2)
eqn2 = RR(RE2, Rc1/(1+beta) + Vt*RE2/Ve2) == 50
Re2 = double(solve(eqn2, RE2));
Re2 = Re2(2)
Ie2 = Ve2/Re2
Ib2 = Ie2/(1+beta)
Ic2 = beta * Ib2
I2 = Ic1 + Ib2
clear Rc1
Rc1 = (Vcc-Vc1)/Ic2
Re1 = Ve1/Ie1
RBB = RR(Rb1, Rb2)

gm1 = Ic1/Vt
rpi1 = beta/gm1
gm2 = Ic2/Vt
rpi2 = beta/gm2

Rcc1_sc = RR(rpi1/(1+beta),Re1)
Rcb_sc = RR(RBB,rpi1)
Rcb_oc = RR(RBB, rpi1 + (1+beta)*Re1)
Rcc2_sc = RR(Re2, (rpi2+Rc1)/(1+beta))

syms C
eqn1 = sqrt((1/(49.5*C))^2 + (1/(36.6*C))^2) == 1000*2*pi;
C = double(solve(eqn1, C));
Cc1 = C(2)
Cc2 = C(2)

wlp1 = 1/(Rcc1_sc*Cc1)
wlp2 = 1/(Rcc2_sc*Cc2)
wl3dB = sqrt(wlp1^2+wlp2^2)
fl3dB = wl3dB/(2*pi)

Cb = 1/(Rcb_oc*wl3dB/10)

%% B
Rin = 706.2e-6/12.92e-6