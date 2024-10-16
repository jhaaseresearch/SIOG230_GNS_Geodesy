function rho_e = get_satpos(t, sv, eph)
    % Computes the satellite position in ECEF coordinates given time, satellite number, and ephemeris data
    % Inputs:
    %   t   - Time in seconds (since the start of the week)
    %   sv  - Satellite number
    %   eph - Ephemeris matrix from RINEX file
    
    % Find the columns corresponding to the satellite number 'sv' in the eph matrix
    sv_column = find(eph(1, :) == sv);
    
    % Extract toe values corresponding to the satellite 'sv'
    toe_values = eph(18, sv_column);
    %% Find the closest valid toe
    % Find the toe values that are less than or equal to t and greater than 2 hours (7200 seconds)
    toe_indeces = find(toe_values <= t & toe_values >= 7200);
    
    % Check if any valid toe is found
    if isempty(toe_indeces)
        error('No valid toe found before or equal to the specified time.');
    else
        % Get the last valid toe value found before or equal to t
        toe = toe_values(toe_indeces(end));
    end
    
    %% Extract satellite-specific ephemeris values at the corresponding toe index
    M0 = eph(3, sv_column(toe_indeces(end)));
    roota = eph(4, sv_column(toe_indeces(end)));
    deltan = eph(5, sv_column(toe_indeces(end)));
    ecc = eph(6, sv_column(toe_indeces(end)));
    omega = eph(7, sv_column(toe_indeces(end)));
    Cuc = eph(8, sv_column(toe_indeces(end)));
    Cus = eph(9, sv_column(toe_indeces(end)));
    Crc = eph(10, sv_column(toe_indeces(end)));
    Crs = eph(11, sv_column(toe_indeces(end)));
    i0 = eph(12, sv_column(toe_indeces(end)));
    idot = eph(13, sv_column(toe_indeces(end)));
    Cic = eph(14, sv_column(toe_indeces(end)));
    Cis = eph(15, sv_column(toe_indeces(end)));
    Omega0 = eph(16, sv_column(toe_indeces(end)));
    Omegadot = eph(17, sv_column(toe_indeces(end)));
    
    %% Compute Time since toe (tk)
    tk = t - toe;
    
    %% Mean anomaly at t (µ)
    GM = 3.986004418e14;  % Gravitational constant times Earth's mass
    a = roota^2;  % Semi-major axis
    mu = M0 + (sqrt(GM / a^3) + deltan) * tk;
    
    %% Iteratively solve for E (eccentric anomaly) with a convergence tolerance of 0.001
    E = mu;  % Start with E = µ
    for iter = 1:100000
        E_new = mu + ecc * sin(E);  % Update equation
        if abs(E_new - E) < 1e-11  % Convergence tolerance
            break;
        end
        E = E_new;  % Update E
    end
    
    %% Compute true anomaly (v)
    v = atan2(sqrt(1 - ecc^2) * sin(E), cos(E) - ecc);
    
    %% Correct for orbital perturbations
    % Argument of perigee correction
    omega = omega + Cuc * cos(2 * (omega + v)) + Cus * sin(2 * (omega + v));
    
    % Radial distance correction
    r = a * (1 - ecc * cos(E)) + Crc * cos(2 * (omega + v)) + Crs * sin(2 * (omega + v));
    
    % Inclination correction
    i = i0 + idot * tk + Cic * cos(2 * (omega + v)) + Cis * sin(2 * (omega + v));
    
    %% Compute right ascension, accounting for Earth's rotation
    omega_e = 7.2921151467e-5;  % Earth's mean angular velocity in rad/s
    Omega = Omega0 + (Omegadot - omega_e) * tk - omega_e * toe;
    
    %% Satellite position in the orbital plane (r, v from previous steps)
    r_orbital = [r * cos(v);
                 r * sin(v);
                 0];
    
    %% Compute the rotation matrix R
    R = [cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i), -cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i), sin(Omega) * sin(i);
         sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i), -sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i), -cos(Omega) * sin(i);
         sin(omega) * sin(i), cos(omega) * sin(i), cos(i)];
    
    % Apply the rotation to the satellite position in the orbital plane
    rho_e = R * r_orbital;
    
    %% Displaying computed values for comparison (optional)

%     disp(['Time t: ', num2str(t), ' seconds']);
%      disp(['Time elapsed tk: ', num2str(tk), ' seconds']);
%      disp(['Time of ephemeris toe: ', num2str(toe), ' seconds']);
%     
% %   Display the satellite position in the orbital plane with units (meters)
%     disp(['Satellite position in orbital plane: ', num2str(r_orbital(1)), ' meters, ', ...
%       num2str(r_orbital(2)), ' meters, ', num2str(r_orbital(3)), ' meters']);
% 
% %   Display the satellite position in ECEF (XYZ) with units (meters)
%     disp(['Satellite position in ECEF (XYZ): ', num2str(rho_e(1)), ' meters, ', ...
%       num2str(rho_e(2)), ' meters, ', num2str(rho_e(3)), ' meters']);
    
% %   Display mean anomaly, eccentric anomaly, true anomaly, right ascension, and other parameters
%     disp(['M0 = ', num2str(M0)]);
%     disp(['Mean anomaly (M) = ', num2str(mu)]);
%     disp(['Eccentric anomaly (E) = ', num2str(E)]);
%     disp(['True anomaly (v) = ', num2str(v)]);
%     disp(['Right ascension (Omega) = ', num2str(Omega)]);
%     disp(['Argument of perigee (omega) = ', num2str(omega)]);
%     disp(['Radial distance (r) = ', num2str(r)]);
%     disp(['Inclination (i) = ', num2str(i)]);
%     
%     % Display the rotation matrix R
%     disp('Rotation matrix (R): ');
%     disp(R);
end
