clearvars
close all
clc

eph = read_rinexn('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/brdc2910.19n');

%% make function [X,Y,Z]=get_satpos(t, sv, eph) 

% other given parameters 
% clock drift rate - row 2
% clock bias - 19 
% clock drift - 20 
% toc - 21 
% tgd - 22
% trans - 23

function [X, Y, Z] = get_satpos(t, sv, eph)
    % parameters needed 

    svprn = eph(1,:);
    M0 = eph(3,:);
    roota = eph(4,:);
    deltan = eph(5,:);
    ecc = eph(6,:);  
    omega0 = eph(7,:);
    cuc = eph(8,:);
    cus = eph(9,:);
    crc = eph(10,:);
    cic = eph(14,:);
    cis = eph(15,:);
    i0 = eph(12,:);
    idot = eph(13,:);
    Omega0 = eph(16,:);
    Omegadot = eph(17,:);
    toe = eph(18,:);

    %find correct time of toe
    t_convert = t * 24 * 3600; % assuming t is in days
    valid_toe = toe(toe <= t_convert);
    closest_toe = max(valid_toe);

    %computing basic parameters 

    tk = t_convert - closest_toe; %time elapsed since toe 
    GM = 3.986004418 * 10^(14); % m^3 s^-2 
    a = roota.^2 % compute semi-major axis
    mu = (M0 + sqrt(GM./a.^3) + deltan) * tk;
    
    % iterative solution for E 
    E = mu;  % initial guess 
    tol = 1e-11;  % tolerance for convergence
    max_iter = 100;  % Maximum number of iterations

    for i = 1:max_iter;
        E_new = mu + ecc .* sin(E);  
        if abs(E_new - E) < tol;
            break  % break loop if reached tolerance
        end
        E = E_new;  
    end

    v = atan2(sqrt(1 - ecc.^2) .* sin(E), cos(E) - ecc);

    % correcting for orbital perturbations
    omega = omega0 + cuc .* cos(2 * (omega0 + v)) + cus .* sin(2 * (omega0 + v)); % argument of perigee
    r = a .* (1 - ecc .* cos(E)) + crc .* cos(2 * (omega0 + v)) + crc .* sin(2 * (omega0 + v)); % radial distance
    incl = i0 + idot * tk + cic .* cos(2 * (omega0 + v)) + cis .* sin(2 * (omega0 + v)); % inclination
 
    % compute the right ascention 
    omegae = 7.2921151467 * 10^-5; % mean angular vel of earth (rad/s)
    
    Omega = Omega0 + (Omegadot - omegae) * tk - omegae*closest_toe;

    % oribtal frame 
    r = [r .* cos(v); r .* sin(v); zeros(size(r))];

    % rotation matrix
    
    R = [cos(Omega).*cos(omega)-sin(Omega).*sin(omega).*cos(incl), -cos(Omega).*cos(omega)-sin(Omega).*sin(omega).*cos(incl), sin(Omega).*sin(incl);
        sin(Omega).*cos(omega)+cos(Omega).*sin(omega).*cos(incl), -sin(Omega).*sin(omega)+cos(Omega).*cos(omega).*cos(incl), -cos(Omega).*sin(incl);
        sin(omega).*sin(incl), cos(omega).*sin(incl), cos(incl)]; 
 
    % applying rotation to get ECEF coords 
    pos_ecef = R*r 

    % assign output variables 

end

[X , Y , Z] = get_satpos(6, 31, eph)
print('sat pos in orbital in X Y Z: ', X ,Y, Z)