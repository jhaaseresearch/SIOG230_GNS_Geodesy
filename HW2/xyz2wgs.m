function output = xyz2wgs(input)
    % This function converts locations from xyz to WGS84 coordinates.
    % Input: N x 4 array of [time(sec), x(m), y(m), z(m)]
    % Output: N x 4 array [time(sec), lon(deg), lat(deg), height(m)].

    % Conventional constants for WGS84
    a = 6378137.0000; % semi-major axis
    f = 1.0 / 298.257223563; % flattening
    b = a * (1 - f); % semi-minor axis
    e2 = 2 * f - f^2; % first eccentricity squared
    e1 = (a^2 - b^2) / b^2; % second eccentricity squared

    % Number of rows
    N = size(input, 1);
    
    % Initialize output matrix
    output = zeros(N, 4);
    fprintf('Time (s)\tLatitude (deg)\tLongitude (deg)\tElevation (m)\n');
    for i = 1:N
        % Extract t, x, y, z for the i-th row
        ti = input(i, 1);
        xi = input(i, 2);
        yi = input(i, 3);
        zi = input(i, 4);
        
        % Auxiliary quantities
        p = sqrt(xi^2 + yi^2);
        theta = atan2(zi * a, p * b);
        
        % Longitude
        lamda = atan2(yi, xi);
        
        % Latitude
        phi = atan2(zi + (sin(theta))^3 * e1 * b, p - (cos(theta))^3 * e2 * a);
        
        % Elevation
        h = p * cos(phi) + zi * sin(phi) - a * sqrt(1 - e2 * sin(phi)^2);
        
        % Store time, latitude, longitude, height in output
        output(i, 1) = ti; % Time (seconds)
        output(i, 2) = rad2deg(lamda); % Longitude (degrees)
        output(i, 3) = rad2deg(phi); % Latitude (degrees)
        output(i, 4) = h; % Height (meters)
        
        % Print each result
        fprintf('%d\t\t%.4f\t\t%.4f\t\t%.4f\n', ti, rad2deg(phi), rad2deg(lamda), h);
    end
end
