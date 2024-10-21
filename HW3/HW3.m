clear all; close all; clc;

% add matlab_utils provided by Jennifer Haas to the path file
addpath('/Users/tgstone/Documents/UCSD/classes/SIOG239_GNSS/SIOG230_GNS_Geodesy/matlab_utils')

% Path to data file
rixen = '/Users/tgstone/Documents/UCSD/classes/SIOG239_GNSS/SIOG230_GNS_Geodesy/products/epgga9.292';
sp3file = '/Users/tgstone/Documents/UCSD/classes/SIOG239_GNSS/SIOG230_GNS_Geodesy/products/igs20756.sp3';

% confirming with Sp3 file

% Reading Rixen file
[eph,alpha,beta,dow] = read_rinexn(rixen);

[sp3,sv,excl_sat] = read_sp3(sp3file);

%SP3 file is in 
% T, X,Y,Z, dt
prn1_data = sp3.prn1;

% Set time for for Satellite
svn = 1;

% Grabbing all the times, in time of day (not week)
ts = prn1_data(:,1)+  6*86400;

% Testing function with Satellite 1 and 10
% 96 x 3
rho_es = get_satpos(ts, svn,eph);
rho_es = rho_es';


% Comparing the differences for sp3 
% 96x3
XYZ_sp3 = prn1_data(:,2:4);

%mag_sp3 = vecnorm(XYZ_sp3,2,2);
%mag_rho = vecnorm(rho_es,2,2);
XYZ_diff = XYZ_sp3 - rho_es;

% sum difference
resid = sum(XYZ_diff.^2,2);


% Taking size
[M,N] = size(XYZ_sp3);

%diff = mag_sp3/M - mag_rho/M;

figure();
plot(ts,resid,'-o');
title('Comparison of Difference magnitude for SP3 and Rho');
xlabel("time (s)")
ylabel("Norm Difference (mag_1 - mag_2) (m)")


% Comments


