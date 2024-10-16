%% SIOG 239 Homework 1 
function results = xyz2wgs(m)
    a = 6378137.0000;  % semi-major axis (meters)
    f = 1.0 / 298.257223563;  % flattening
    b = a * (1 - f);  % semi-minor axis (not used directly in this version)
    e2 = 2 * f - f^2;  % first eccentricity squared
    
    % Initialization
    N = size(m, 1);  
    results = zeros(N, 4);
    
    for i = 1:N
        % t, x, y, z for the current position
        t = m(i, 1);
        x = m(i, 2);
        y = m(i, 3);
        z = m(i, 4);
        
        % Auxiliary quantities
        p = sqrt(x^2 + y^2);  
        
        % Longitude (lambda)
        lambda = atan2(y, x);
        
        % Latitude (phi) 
        phi = atan2(z, p * (1 - e2));
        
        % Radius of curvature in the prime vertical (N)
        N_phi = a / sqrt(1 - e2 * sin(phi)^2);
        
        % Elevation (h)
        h = (p / cos(phi)) - N_phi;
        
        % Radians to degrees 
        lon = rad2deg(lambda);
        lat = rad2deg(phi);
        
     
        results(i, :) = [t, lat, lon, h];
    end
end

% checking 
N = 5;  % Number of m
x_val = 4433469.9438;
y_val = 362672.7267;
z_val = 4556211.6409;

positions = [(1:N)', repmat([x_val, y_val, z_val], N, 1)];

results = xyz2wgs(positions);

for i = 1:N
    fprintf('Time: %d, Latitude: %.4f, Longitude: %.4f, Height: %.4f\n', ...
            results(i, 1), results(i, 2), results(i, 3), results(i, 4));
end
