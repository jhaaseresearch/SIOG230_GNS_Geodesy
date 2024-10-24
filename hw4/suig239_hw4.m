clearvars
close all
clc
[obs,t,prn,apr,hant] = read_rinexo('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/opmt2920.19o');

% Output: obs = structure with observation types, e.g. obs.C1, obs.L1, etc.
%         t = time matrix (year month day hour minute second)
%         gps = vector of GPS satellite numbers
%         apr = approximate site position from rinex file header
%% 4 
f1 = 1.57542 * 10^9; %Hz
f2 = 1.2276 * 10^9; %Hz
c = 0.299792458 * 10 ^9; % m/s
lambda1 = c/f1;
lamda2 = c/f2;

P1 = obs.P1;
P2 = obs.P2 ;
alpha = (f1/f2)^2;
L1 = obs.L1;
L2 = obs.L2; 

% For PRN10 
PRN10 = prn(9,:);

MP1 = P1 - ((2/(alpha-1))+1)*L1 + (2/(alpha-1))*L2;
MP2 = P2 - (2*alpha/(alpha-1))*L1 + (2*alpha/(alpha-1)-1)*L2;

PRN10_MP1 = MP1(:, 9);
PRN10_MP2 = MP2(:, 9);

plot(PRN10_MP1)
hold on
plot(PRN10_MP2)