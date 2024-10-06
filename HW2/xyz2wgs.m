function B = xyz2wgs(A)
%{
xyz2wgs is a function created for HW2 for SIOG 239

Author: Tommy Stone
Date: 10/05/24
HW2

Parameters
----------
A - array
    Input Array of N x (time (s), x (m), y (m), z(m)). 

Returns

B - array
    Output array of N x (time (s), longitude (degree), latitude (degree),
    height (m))
%}

% Defining constants for WGS84
a = 6378137.0000;         % semi-major axis
f = 1.0/298.257223563;    % flattening
b = a*(1-f);               % semi-minor axis
e2 = 2*f - f^2;           % eccentricity squared
e1 = (a^2 - b^2)/(b^2);    % second numerical eccentricity


% Creating a Nx4 output matrix
[n,m] = size(A);

% Initializing new matrix B
B = nan(n,m);

% Time does not change so assigning values from time of A to B
B(:,1) = A(:,1);

% Iterating through each row of A to solve for the values

% Iterating through each row
for t = 1:n
    
    % Assigning x,y,z values
    tmp_xyz = A(t,2:end);
    x = tmp_xyz(1);
    y = tmp_xyz(2);
    z = tmp_xyz(3);

    p = sqrt(x^2 + y^2);      
    theta = atan2d(z*a, p*b);

    % Defining longitude
    lambda = atan2d(y,x);

    % Defining latitude
    phi = atan2d(z + ((sind(theta))^3)*e1*b,p - ((cosd(theta))^3)*e2*a);

    % Defining radius of curvature in prime vertical
    N = a/(1 - (sind(phi)^2)*e2);

    % Elevation
    h = p*cosd(phi) + z*sind(phi) - a*sqrt(1 - e2*(sind(phi)^2));

    % assigning to B to return values
    B(t,2) = lambda;
    B(t,3) = phi;
    B(t,4) = h;
end







end

