clear all
close all
clc


[eph,alpha,beta,dow] = read_rinexn('brdc2920.19n');
%how to define t? -> see f)
% first epoch of day 292/2019, when t = 6 × 24 × 3600:
t = 6;
%how to define sv -> see f)
sv=31;
rho = getsatpos(t, sv, eph)









% 
% rho_es = get_satpos(t,1,eph);
% rho_es = get_satpos(t,30,eph)eph = read_rinexn();