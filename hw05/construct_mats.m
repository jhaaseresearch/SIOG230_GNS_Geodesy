function [A, L] = construct_mats(P, sv_arr, ephBrdc, apr, non_nan_row, tidx, gps)
% construct_mats   Construct L and design  matrix A
%     In Path: get_satpos.m
%     Input:
%       - P = obs.C1 from obs, from read_rinexo.m
%       - sv_arr = satellite PRN INDEX array for the non-NaN times
                % key to remember this is NOT SV NUMBER
%       - non_nan_row = valid times of transmission
%       - ephBrdc = ephemerides matrix created from a rinex navigation
%               file using read_rinexn.m
%       - apr = approx. positions [X0, Y0, Z0] from read_rinexo.m
%       - gps = gps numbers from read_rinexo.m

X0 = apr(1);
Y0 = apr(2);
Z0 = apr(3);
c = 0.299792458*1e9;

% Create L vector (sv x 1)
L = zeros(length(sv_arr), 1);
% Design matrix A
A = zeros(length(sv_arr), 4);

for i = 1:length(sv_arr)
    t0 = non_nan_row(i);
    % NOTE: sv_arr =/= gps
    % For this specific satellite number j = sv:
    % sat pos is computed at a t_tx, for each satellite
    sv = gps(sv_arr(i));
    satpos = get_satpos(t0,sv,ephBrdc,0);
    % output satpos = [X Y Z dt toe drel]
    X = satpos(1);
    Y = satpos(2);
    Z = satpos(3);
    dt = satpos(4);
    
    rho0 = sqrt((X-X0)^2 + (Y-Y0)^2 + (Z-Z0)^2);
    a_x = -((X-X0)/rho0);
    a_y = -((Y-Y0)/rho0);
    a_z = -((Z-Z0)/rho0);
    l = P(tidx,sv_arr(i)) - rho0 + c*dt;
    
    L(i) = l;
    A(i, 1) = a_x;
    A(i, 2) = a_y;
    A(i, 3) = a_z;
    A(i, 4) = -0.299792458; % multiply by 10^-9

    
end

