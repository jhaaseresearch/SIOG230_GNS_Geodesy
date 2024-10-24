[obs,t,prn,apr,hant] = read_rinexo('opmt2920.19o');
%% Given parameters
f1 = 1.57542e9; % Hz
f2 = 1.2276e9;  % Hz
c = 0.299792458e9; % m/s
arc_duration = 8;
prn_target = 10;

plot_MP1_MP2(obs, t, prn, prn_target, hant, f1, f2, c, arc_duration)