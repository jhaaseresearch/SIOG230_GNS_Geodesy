clear all
close all
clc

%constants
f1 = 1.57542 *109;% Hz
f2 = 1.2276 * 109;% Hz
c = 0.299792458 * 109;% m/s 
lambda_1 = c/f1;
lambda_2 = c/f2;

[data,t,prn,apr,hant] = read_rinexo("opmt2920.19o");

idx = find(prn==10);
compute_multipath(f1,f2, data.L1(:,idx), data.L2(:,idx),data.P1(:,idx), data.P2(:,idx))