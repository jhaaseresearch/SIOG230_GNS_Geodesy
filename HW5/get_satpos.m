function  [x,y,z, delta_t] = get_satpos(t_data,sv,eph_all)

    %Parameters for the satellite
    eph_columns = eph_all(1, :) == sv;
    eph = eph_all(:, eph_columns);

    %Delete columns with negative transmission time
    trans_neg = eph(23, :) < 0;
    eph(:, trans_neg) = [];

    %Find the closest time to t_data
    diff = abs(eph(18,:) - t_data);
    [~, ind] = min(diff);
    closest_time = eph(18,ind);

%    for i = 1:size(eph, 2) uncomment it when calculating for multiple times
    for i = ind
        Mo = eph(3,i);
        roota = eph(4,i);
        deltan = eph(5,i);
        ecc = eph(6,i);
        omega0 = eph(7,i);
        cuc = eph(8,i);
        cus = eph(9,i);
        crc = eph(10,i);
        crs = eph(11,i); 
        cic = eph(14,i);
        cis = eph(15,i);
        i0 = eph(12,i);
        idot = eph(13,i);
        Omega0 = eph(16,i);
        Omegadot = eph(17,i);
        toe = eph(18,i);
        af0 = eph(19,i); %clock_drift
        af1 = eph(20,i); %clock_bias;
        af2 = eph(2,i); %clock_drift_rate
        
        %%Compute basic parameters at request time t_data        
        %Time elapsed since toe
        t = t_data - toe;
        delta_t = af0 + af1*(t_data-toe) + af2*(t_data-toe)^2
        
        %Mean anomaly at t
        GM = 3.986004418 * 10^14;
        a =roota^2
        mu = Mo + (sqrt(GM/(a^3)) + deltan) * t;
    
        %Iterative solution for E
        e = ecc;
        E = mu; % Initial guess
        tolerance = 1e-11; 
        maxIter = 1000; 
        iter = 0; 

        while iter < maxIter
            E_old = E;
            E = mu + e * sin(E); 
            iter = iter + 1; 

            if abs(E - E_old) < tolerance
                break;
            end
        end
    
        %True anomaly
        v = atan2(sqrt(1 - (e^2)) * sin(E), cos(E) - e);
    
        %%Correct for orbital pertubations
        %Argument of perigee
        omega = omega0 + cuc * cos(2 * (omega0+v)) + cus * sin(2 * (omega0+v));
    
        %Radial Distance
        r= a * (1- e * cos(E)) + crc * cos(2 * (omega0 + v)) + crs * sin(2 * (omega0 + v));
    
        %Inclination
        inc= i0 + idot * t + cic * cos(2 * (omega0 + v)) + cis * sin(2 * (omega0 + v));
    
        %%Compute the right ascension
        we = 7.2921151467*10^(-5)
        w = Omega0 + (Omegadot - we) * t - we * toe;
    
        %%Convert satellite position from orbital frame to ECEF frame
        % Orbital position in the orbital frame
        x_orb = r * cos(v);
        y_orb = r * sin(v);
        z_orb = 0;
    
        % Rotation matrix from orbital frame to ECEF frame
        R = [
            cos(w) * cos(omega) - sin(w) * sin(omega) * cos(inc), ...
            -cos(w) * sin(omega) - sin(w) * cos(omega) * cos(inc), ...
            sin(w) * sin(inc);
            
            sin(w) * cos(omega) + cos(w) * sin(omega) * cos(inc), ...
            -sin(w) * sin(omega) + cos(w) * cos(omega) * cos(inc), ...
            -cos(w) * sin(inc);
            
            sin(omega) * sin(inc), ...
            cos(omega) * sin(inc), ...
            cos(inc)
        ];
    
         % Create the position vector in the orbital frame
        r_orbital_vector = [x_orb; y_orb; z_orb];
    
        % Convert to ECEF
        r_ecef = R * r_orbital_vector;

%        xyz(i, :) = r_ecef uncomment it when calculating for multiple times
        x = r_ecef(1)
        y = r_ecef(2)
        z = r_ecef(3)
    end
end
