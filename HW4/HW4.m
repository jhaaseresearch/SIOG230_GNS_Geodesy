%{
HW4

HW 2.3.4 from geodesy lab book, page 82

Author: Tommy Stone
Class: SIOG 239
Date: 10/21/24
%}

clear all; close all; clc;

% Problems 1-5
% Adding utils created for class to this Matlab script
addpath('../matlab_utils')

% Adding products
addpath('../products')

% File added from HW
rinex_file = 'opmt2920.19o';

[obs,t,gps,apr,hant] = read_rinexo(rinex_file);

% converting to datetime
times = datetime(t(:,1)+2000,t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
refTime= datetime(1900,1,1,0,0,0);

% created reference time in seconds
timeSeconds = seconds(times - refTime);

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
prn = 10;
index_pos = find(gps == prn);

% Filtering by Satellite
L1 = L1_obs(:,index_pos);
L2 = L2_obs(:,index_pos);
P1 = P1_obs(:,index_pos);
P2 = P2_obs(:,index_pos);

% Removing NaNs from file
non_nan_index = find(~isnan(L1));

% PRN10 has several discontinuous areas, finding one continuous section
non_nan_index = non_nan_index(1:600);

% Removing NaN values
L1 = L1(non_nan_index);
L2 = L2(non_nan_index);
P1 = P1(non_nan_index);
P2 = P2(non_nan_index);

% Filtering times
times = times(non_nan_index);
timeSeconds = timeSeconds(non_nan_index);

% Obtaining distance
L1 = L1*lambda1;
L2 = L2*lambda2;

% Calculating MP1 and MP2
p1_const = 2/(alpha - 1);
p2_const = (2*alpha)/(alpha - 1);
MP1 = P1 - (p1_const + 1).*L1 + p1_const.*L2;
MP2 = P2 - p2_const.*L1 + (p2_const - 1).*L2;

% Taking the mean from the plots
MP1_mean = mean(MP1);
MP2_mean = mean(MP2);

% subtrackign mean
MP1_nomean = MP1 - MP1_mean;
MP2_nomean = MP2 - MP2_mean;

% Plotting results
figure()
subplot(2,1,1)
plot(times,MP1_nomean,'');
title("MP1 for PRN" + num2str(prn));
xlabel("Time (hours)")
ylabel("MP1 (m)")
grid();
datetick;
subplot(2,1,2)
plot(times,MP2_nomean,'r');
xlabel("Time (hours)")
ylabel("MP2 (m)")
grid()
title("MP2 for PRN" + num2str(prn))
datetick;

% 6) Discuss the results
%{
In order to obtain the results from the plot I had to remove the NaN values
as well as choose a time segment that had continuous data. This appeared to
be from around 4:00-8:30 on 10/19/19. That both MP1 and MP2 had varying
values at the start of the record indicating potentially higher
interference. As we move further down the time record, roughly 4:30 onward,
we see that the values are constrainted between 1 and -1 meters indicated
that the Satellite may be over the antenna directly causing less
intereference and a better direct path

7) The noise ambiguity was seen as a constant value which we take to be the
mean of the MP1 and MP2 values which is why this was subtracted. Also, the
Ionosphere was estimated using the frequency values but MP1 and MP2 could
be better estimated with a model of the Ionosphere. We also do not have the
Troposphere in our data which can play a role. 
%}



