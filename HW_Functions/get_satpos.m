function rho_es = get_satpos(t,sv,eph)
%{
HW3
Class: SIOG 239
Author: Tommy Stone

Function returns the Satellite position in XYZ coordinates

Parameters
----------
t - time, in seconds of GPS week
sv - Satellite number
eph - ephemeris file

Returns
-------
rho_es - XYZ vector for Satellite positions

%}

    % Constants
    GM = 3.986004418e14; % Gravitational constant x mass of the earth (m^3/s^2)
    
    % Description of eph Data x Satellite
    % 1) svprn - Satellite PRN Number
    % 2) Clock_drift_rate - always 0
    % 3) Mo - Mean anomaly, radians [mu_0]
    % 4) Roota - sqrt(semi-major axis), sqrt(m)|sqrt(a)|
    % 5) deltan - variation of mean angular velocity, radians/s [\delta n]
    % 6) ecc - eccentricity [e]
    % 7) omega0 - arguement of perigee, raidians [\omega_0]
    % 8) cuc - correction coefficient, m or radians
    % 9) cus - correction coefficient, m or radians
    % 10) crc - correction coefficient, m or radians
    % 11) crs - correction coefficient, m or radians
    % 12) i0 - inclination, radiasn [i0]
    % 13) idot - rate of inclination, radians/s [i]
    % 14) cic - correction coefficients, m or radians
    % 15) cis - correction coefficients, m or radians
    % 16) Omega0 - right ascension, radiasn [Omega_0]
    % 17) Omegadot - rate of right ascension, radians/s [\Omega^\dot]
    % 18) toe - time of ephemeris, seconds of the current GPS week [t_oe]
    % 19) clock_bias
    % 20) clock_drift
    % 21) toc - 
    % 22) tgd - transmitter group delay
    % 23) trans - tranmission time of message
    
    % filtering data by Satellite number
    satellite_index = eph(1,:) == sv;
    epha = eph(:,satellite_index);
    
    %####### part (b)##############
    
    % Find the elapsed time tk = t - toc (eq 2.5)
    toes = epha(18,:);
    
    % For now setting Toe to the same time as t
    % This could be changed later
    tks = t - toes;                                           % eq 2.5
    
    % Determine mean anomaly, eq 2.6
    mos = epha(3,:);   % mu_0 vector
    as = epha(4,:).^2;  % a term is square root so squaring term
    a3s = as.^3;        % per equation take the third power
    dns = epha(5,:);    % setting dn values
    mus = mos + (sqrt(GM./a3s) + dns).*tks;                    % eq 2.6
    
    
    % Iteratively detemining values of E
    tol = 1e-11;  % tolerance given by Jennifer Haas
    diff = 1;     % setting initial tolerance
    ii = 1;       % displaying iteerations
    es = epha(6,:);
    Es = mus + es.*sin(mus);
    
    % iteratively finding Es                                 % eq 2.7
    while diff >= tol

        if ii >= 100
            sprinf("Reached max iterations")
            break;
        end
        
        Es_next = mus + es.*sin(Es);
    
        % taking difference
        diff = Es - Es_next;
    
        % finding largest error from the vector
        diff = abs(max(diff));
        Es = Es_next;
        ii = ii + 1;
    end
    sprintf("Iterations for convergence for E: %i",ii)
    sprintf("Difference: %e", diff)
    
    % Finding the true anomaly, v, equation 2.8
    vs = atan((sqrt(1 - es.^2).*sin(Es)./(cos(Es) - es)));  % eq. 2.8
    
    
    %######## part (c) #############
    
    %argument of perigree
    omega0s = epha(7,:);
    Cucs = epha(9,:);
    
    omegas = omega0s + Cucs.*cos(2.*(omega0s + vs)) + ...   % eq 2.9
        Cucs.*sin(2.*(omega0s + vs));
    
    % radial distance
    Crcs = epha(10,:);
    rs = as.*(1 - es.*cos(Es)) + ...                        % eq. 2.10
        Crcs.*cos(2.*(omega0s + vs)) + ...
        Crcs.*sin(2.*(omega0s + vs));  
    
    % inclination
    i0s = epha(12,:);
    idots = epha(13,:);
    Cics = epha(14,:);
    Ciss = epha(15,:);
    
    is = i0s + idots.*tks + ...                             % eq 2.11
        Cics.*cos(2.*(omega0s + vs)) + ...
        Ciss.*sin(2.*(omega0s + vs));
    
    %######### part (d) ############
    
    omega_e = 7.2921151467e-5; % (rad/s)
    Omega0s = epha(16,:);
    Omegadots = epha(17,:);
    Omegas = Omega0s + (Omegadots + omega_e).*tks - omega_e.*toes; % eq 2.12
    
    %############ part (e) ########## 
    % convert Satellite position from coordinates in orbital plane to 
    % coordinates in ECEF frame
    
    % matrix of positions
    % Initializing Rs
    N = length(rs);
    r_pos = nan(3,N);
    
    % converting to orbital frame                           % eq 2.13
    for ii = 1:N
        r = rs(ii);
        v = vs(ii);
        r_pos(:,ii) = [r.*cos(v); r.*sin(v); 0];
    end
    
    % Build matrix to apply rotation
    Rs = nan(3,3,N);
    
    for ii = 1:N                                            % equation 2.14
        % Assigning values use at each iteration
        Omega = Omegas(ii);
        omega = omegas(ii);
        ini = is(ii);
        pos11 = cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(ini);
        pos12 = -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(ini);
        pos13 = sin(Omega)*sin(ini);
    
        pos21 = sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(ini);
        pos22 = -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(ini);
        pos23 = -cos(Omega)*sin(ini);
    
        pos31 = sin(omega)*sin(ini);
        pos32 = cos(omega)*sin(ini);
        pos33 = cos(ini);
    
        R_tmp = [pos11 pos12 pos13; pos21 pos22 pos23; pos31 pos32 pos33];
        Rs(:,:,ii) = R_tmp;
    
    end
    
    % finding position
    rho_es = nan(3,N);                                 % eq 2.15
    for ii = 1:N
        rho_es(:,ii) = Rs(:,:,ii)*r_pos(:,ii);
    end
    
end    
    
    
    
    

