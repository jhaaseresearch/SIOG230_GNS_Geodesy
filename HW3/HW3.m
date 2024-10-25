% Neil Waldhausen
% SIOG 239 Homework 2
% Converting RINEX data file to ECEF coordinates

% Load ephemeris data
eph = read_rinexn("brdc2920.19n");

%% Read data from RINEX file
svw = eph(:, 31);   % Chosen Satellite
[i, j] = find(eph(1,:), 1);  % Ensure this works as intended

% Extract satellite parameters
svprn = svw(1, j);
clock_drift_rate = svw(2, j);
M0 = svw(3, j);
roota = svw(4, j);
deltan = svw(5, j);
ecc = svw(6, j);
omega = svw(7, j);
cuc = svw(8, j);
cus = svw(9, j);
crc = svw(10, j);
crs = svw(11, j);
i0 = svw(12, j);
idot = svw(13, j);
cic = svw(14, j);
cis = svw(15, j);
Omega0 = svw(16, j);
Omegadot = svw(17, j);
toe = svw(18, j);
clock_bias = svw(19, j);
clock_drift = svw(20, j);
toc = svw(21, j);
tgd = svw(22, j);
trans = svw(23, j);

% Find tk
tk = svw(18, j);

%% Calculations %%
% Mean Anomaly
gm = 3.986004418e14;  % Gravitational constant
mu = M0 + (sqrt(gm / (roota^2)) + deltan) * tk;

% Iterative Solution for Eccentric Anomaly
tol = 1e-11;
E_start = 0;  % Initial guess for E
E = mu;  % First approximation of E

% Iterate to solve for E
while true
    E_prev = E;
    E = mu + ecc * sin(E_prev);
    if abs(E - E_prev) < tol
        break;
    end
end

% True Anomaly
v = atan2(sqrt(1 - ecc^2) * sin(E), cos(E) - ecc);

% Correct for Orbital Perturbations
% Argument of Perigee
omega3 = omega + cuc * cos(2 * (omega + v)) + cus * sin(2 * (omega + v));

% Radial Distance
rad = (roota^2) * (1 - ecc * cos(E)) + crc * cos(2 * (omega + v)) + crs * sin(2 * (omega + v));

% Inclination
inc = i0 + idot * tk + cic * cos(2 * (omega + v)) + cis * sin(2 * (omega + v));

% Right Ascension
we = 7.2921151467e-5;  % Earth's rotation rate
omega4 = Omega0 + (Omegadot - we) * tk - we * toe;

% Position in Orbital Plane
rv = [rad * cos(v); rad * sin(v); 0];

% Rotation Matrix to ECEF Coordinates
R = [
    cos(omega4) * cos(omega3) - sin(omega4) * sin(omega3) * cos(inc), ...
    -cos(omega4) * sin(omega3) - sin(omega4) * cos(omega3) * cos(inc), ...
    sin(omega4) * sin(inc);
    
    sin(omega4) * cos(omega3) + cos(omega4) * sin(omega3) * cos(inc), ...
    -sin(omega4) * sin(omega3) + cos(omega4) * cos(omega3) * cos(inc), ...
    -cos(omega4) * sin(inc);
    
    sin(omega3) * sin(inc), ...
    cos(omega3) * sin(inc), ...
    cos(inc)
];

% Convert to ECEF coordinates
rho = R * rv;

% Output ECEF coordinates
X = rho(1);
Y = rho(2);
Z = rho(3);
