%{
HW4

HW 2.3.4 from geodesy lab book, page 82

Author: Tommy Stone
Class: SIOG 239
Date: 10/21/24
%}

clear all; close all; clc;

% Problems 1-3
% Adding utils created for class to this Matlab script
addpath('../matlab_utils')

% Adding products
addpath('../products')

% File added from HW
rinex_file = 'opmt2920.19o';

[obs,t,gps,apr,hant] = read_rinexo(rinex_file);

% Problems 4, setting constants
f1 = 1.57542e9; % Hz
f2 = 1.2276e9;  % Hz
c = 0.299792458e9; % m/s
lambda1 = c/f1;
lambda2 = c/f2;

% Setting alpha value
alpha = (f1/f2)^2;

% retrieving L1,L2 (phase measurement) and P1,P2 (pseudorange measurement)
L1_obs = obs.L1;
P1_obs = obs.P1;
L2_obs = obs.L2;
P2_obs = obs.P2;

% Filtering fo prn10, first using prn2 because prn10 has NaN values
prn = 2;
index_pos = find(gps == prn);

% Currently this is providing NaN files
L1 = L1_obs(:,index_pos);
L2 = L2_obs(:,index_pos);
P1 = P1_obs(:,index_pos);
P2 = P2_obs(:,index_pos);

% Calculating MP1 and MP2
p1_const = 2/(alpha - 1);
p2_const = (2*alpha)/(alpha - 1);
MP1 = P1 - (p1_const + 1).*L1 + p1_const.*L2;
MP2 = P2 - p2_const.*L1 + (p2_const - 1).*L2;

% Plotting MP1 and MP2 on same graph
%
%(\___/)
%(=^.^=) Currently, if lines 45-49 are correct then there are a bunch of
%(")_(") NaN values in the file, Need to determine what to do with NaNs
%
figure()
plot(MP1);
hold on;
plot(MP2);



