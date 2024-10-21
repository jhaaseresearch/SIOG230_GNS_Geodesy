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

% Filtering fo prn10
prn = 10;
index_pos = find(gps == prn);

% Currently this is providing NaN files
L1 = L1_obs(:,index_pos);
L2 = L2_obs(:,index_pos);
P1 = P1_obs(:,index_pos);
P2 = P2_obs(:,index_pos);


