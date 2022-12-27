%% PART 3
Vcc = 15; Vbe = 0.7; beta = 300; alpha = beta/(beta+1); Vt = 0.025; Rc = 10e3; Rs = 50;
%% A
Io = 2e-3
Iref = Io * (1+2/beta)
Rref = (0-(-Vcc+Vbe))/Iref

%% B
gm = alpha*Io/(2*Vt)
Vcb = 3.3

CJE = 4.5e-12;
TF = 400e-12;
cpi = 2*CJE + TF*gm %#ok<*NOPTS> 

CJC = 3.5e-12;
VJC = 750e-3;
MJC = 330e-3;
Vcb = Vc2-Vb2;
cmu = CJC/(1+Vcb/VJC)^MJC

rpi = beta/gm

Ad = -2*gm*Rc*rpi/(Rs+2*rpi)
whp2 = 1/(cmu/2 * 2*Rc)
fh3dB = whp2/(2*pi)