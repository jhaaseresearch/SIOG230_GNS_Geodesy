clearvars;
close all;
clc;

% 2.4.3 

% constants 
c = 299792458; % speed of light 
angle_tol = 10; % cutoff angle for satellite 

% read psuedorange, eph and sp3 data 
[obs,t,prn,apr,hant] = read_rinexo('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/opmt2920.19o');

[eph, alpha, beta, dow] = read_rinexn('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/brdc2920.19n');

[sp3,sv,exc] = read_sp3('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/igs20756.sp3');

C1 = obs.C1; 

time_day = datetime(t(:,1)+2000,t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
ref_time = datetime(2019, 10, 19,0,0,0);

ts = seconds(time_day - ref_time); % dt in seconds 

% time of transmission 
ttx = ts - C1./c; 

sats_given = [2, 6, 12, 14, 24, 25, 29, 32];

time_num = length(ttx);
sats_num = length(sats_given); 

% a priori pos of reciever 
x0 = apr(1); 
y0 = apr(2); 
z0 = apr(3);


% intitializing matrix 
rho_j = zeros(1, length(sats_num));
L = zeros(1,length(sat_num));
A = zeros(length(sat_num), 4);

% finding filtered indexes 
for k=1:length(sats_given)
    idx(k) = find(prn == sats_given(k));
    ttx = ttx(:, idx);
    C1 = C1(:, idx);
    % sat pos znd clock bias 
    sp = sat_pos(i);
    sat_pos = get_satpos(ttx, sp, eph, 3)
    x = sat_pos(1);
    y = sat_pos(2); 
    z = sat_pos(3);
    dt = sapos(4); 

    % modeled observables 

    rho = qrt((x - x0)^2 + (y - y0)^2 + (z - z0)^2);
    L_k = obs.C1(epoch_index, sp) - rho + c * dt;
    
    rho_j(k)=rho; 
    L(k) = L; 

    % partial derivatives 

    A(i, 1) = (x0 - x)/rho; 
    A(i, 2) = (y0 - y)/rho; 
    A(i, 3) = (z0 - z)/rho; 
    A(i, 4) = -c * 1^-9; 

     
end 

% solving for x by least squares
all_x_pos = []
all_y_pos = []
all_y_pos = []

% covariance 

%station position 

% ellispoisal coord 

%rotation matrix 

% covariance in new frame 

% DOPS 
