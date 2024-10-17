clear all; close all; clc;

% add matlab_utils provided by Jennifer Haas to the path file
addpath('/Users/tgstone/Documents/UCSD/classes/SIOG239_GNSS/SIOG230_GNS_Geodesy/matlab_utils')

% Path to data file
rixen = '/Users/tgstone/Documents/UCSD/classes/SIOG239_GNSS/SIOG230_GNS_Geodesy/HW3/data/brdc0160.12n';

% Reading Rixen file
[eph,alpha,beta,dow] = read_rinexn(rixen);

% Set time for for Satellite
t = 518400;

% Testing function with Satellite 1 and 10
rho_es = get_satpos(t,1,eph)
rho_es = get_satpos(t,30,eph)

