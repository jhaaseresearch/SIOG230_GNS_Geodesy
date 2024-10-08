%%Converting from xyz to WGS84

%a = 6378137.0000 semi major axis
%f = 1/298.257223563 flattening
%b = a(1-f) semi minor axis
%e2 = 2f - f^2 eccentricity squared
%e1 = (a^2 - b^2)/b^2 second numerical eccentricity

%auxiliary quantities
%p = sqr(x^2 + y^2)
%theta = atan2(za, pb)
%lon = atan2(y,x)
%lat = atan2(z + (sin(theta))^3 * e1 * b, p-(cos(theta))^3 * e2 * a
%N = a/sqr(1 - (sin(lat))^2 * e2)
%elev = p * cos(lat) + z * sin (lat) - a * sqr(1-e2*sin(lat)^2)

function wgs = xyz2wgs(txyz)

    % Parameters for WGS84
    a = 6378137
    f = 1 / 298.257223563
    b = a * (1-f)
    e2 = 2*f - f^2
    e1 = (a^2 - b^2)/b^2

    N = size(txyz, 1)
    wgs = zeros(N, 4)

    for i = 1:N
        % Define the input
        time = txyz(i, 1)
        x = txyz(i, 2)
        y = txyz(i, 3)
        z = txyz(i, 4)
        
        % Longitude
        lon = atan2(y, x)
        lon_deg = rad2deg(lon) 

        % Latitude
        p = sqrt(x^2 + y^2)
        theta = atan2(z*a, p*b)      
        lat = atan2(z + (sin(theta))^3 * e1 * b, p-(cos(theta))^3 * e2 * a) 
        lat_deg = rad2deg(lat) 
        
        % Ellipsoidal height
        elev = p * cos(lat) + z * sin(lat) - a * sqrt(1-e2*(sin(lat))^2)
        
        wgs(i, :) = [time, lon_deg, lat_deg, elev];
    end
end