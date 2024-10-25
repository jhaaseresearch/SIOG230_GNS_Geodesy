% Neil Waldhausen
% SIOG 239 Homework 4 (Assigment 2.3.4)
%%
clc, clear all; close all

% Rinex File
rinex_file = 'opmt2920.19o';
[obs,t,gps,apr,hant] = read_rinexo(rinex_file);

%%
%Filter for PRN10
prn = 10;
idx = find(gps == prn);

% Initial Values from Rinex File
obsL1 = obs.L1(:,idx);
obsL2 = obs.L2(:,idx);
obsP1 = obs.P1(:,idx);
obsP2 = obs.P2(:,idx);

% Remove NaN from data
remove_nan = find(~isnan(obsL1));
nan_index = remove_nan(1:600,1);

L1 = obsL1(nan_index);
L2 = obsL2(nan_index);
P1 = obsP1(nan_index);
P2 = obsP2(nan_index);

% Numerical Values for Calculations
f1 = 1.57542e9;
f2 = 1.2276e9;
c = 0.299792458e9;

lambda1 = c/f1;
lambda2 = c/f2;

alpha = (f1/f2)^2;

Final_L1 = L1*lambda1;
Final_L2 = L2*lambda2;

% Filtering Time Values
times = datetime(t(:,1)+2000,t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
refTime= datetime(1900,1,1,0,0,0);

idx_times = times(nan_index);
idx_times_sec = seconds(times(nan_index) - refTime);

% Calculations for Multipath Error for PRN10
 MP1 = P1 - (((2/(alpha-1))+1)*Final_L1)+((2/(alpha-1))*Final_L2);
 MP2 = P2 - (((2*alpha)/(alpha-1))*Final_L1)+(((2*alpha)/(alpha-1)-1)*Final_L2);

meanMP1 = mean(MP1);
meanMP2 = mean(MP2);

MP1_dif = MP1 - meanMP1;
MP2_dif = MP2 - meanMP2;

% Create Plot for MP1 and MP2
figure
subplot(2,1,1)
plot(idx_times, MP1_dif, Color='blue')
title("MP1 for PRN 10")
xlabel("Time (Hours)")
ylabel("MP1 (meters)")
datetick;
grid on

subplot(2,1,2)
plot(idx_times, MP2_dif, Color='red')
title("MP2 for PRN 10")
xlabel("Time (Hours)")
ylabel("MP2 (meters)")
datetick;
grid on