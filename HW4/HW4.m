% Neil Waldhausen
% SIOG 239 Homework 4 (Assigment 2.3.4)
%%
% Rinex File
rinexo = 'opmt2920.19o';
[obs,t,prn,apr,hant] = read_rinexo(rinexo);

%%
% Values from Rinex File
L2
L1
P1
P2


%%

% Numerical Values for Calculations
f1 = 1.57542*10^9;
f2 = 1.2276*10^9;
c = 0.299792458*10^9;
lambda1 = c/f1;
lambda2 = c/f2;

% Calculations for Multipath Error for PRN10

MP1 = P1 - ((2/(alpha-1))+1)*L1+((2/(alpha-1)))*L2;
MP2 = P2 - (((2*alpha)/(alpha-1))+1)*L1+(((2*alpha)/(alpha-1)))*L2;
