% Neil Waldhausen
% SIOG 239 Homeowork 2
% Converting RINEX data file to ECEF coordinates

%[sp3,sv,excl_sat] = read_sp3('igs20756.sp3');
%function [X,Y,Z] = get_satpos(t,sv,eph)

eph = read_rinexn("brdc2920.19n");

%% Read data from rinex file
svw = eph(:, 31);   % Chosen Satellite
[i j] = find(eph(1,:),1);

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

% find tk
tk = svw(18, j);

%% Calculations %%
% Mean Anomoly
gm = 3.986004418*10^14;
mu = M0 + (sqrt(gm/(roota^2))+deltan)*tk;

%Iterative Solution for E
tol = 1e-11;
E_start = 0;
E = mu + ecc .* sin(E_start);

k = 1;  % Index for iteration
for i = k
    E_start = E;
    E = mu + ecc .* sin(E_start);
    if abs(E - E_start) < tol;
        break
    end
    k = k+1;
end

%True Anomoly
v = atan((sqrt(1-2.718^2))*sin(E))/(cos(E)-2.718);

% Correct for Orbital Pertubations %
% Argument of Perigee
omega3 = Omega0+cuc*cos(2*(Omega0+v))+cus*sin(2*(Omega0+v));

% Radial Distance
rad = (roota^2)*(1-2.718*cos(E))+crc*cos(2*(Omega0+v))+crs*sin(2*(Omega0+v));

% Inclination
inc = i0 + tk + cic*cos(2*(Omega0+v))+cis*sin(2*(Omega0+v));

% Right Ascension
we = 7.2921151467*10^-5;
omega4 = Omega0+(Omegadot-we)*tk-we*toe;

% Rotation matrix 
rv = [rad.*cos(v); rad.*sin(v); 0];
R = [cos(omega4).*cos(omega3)-sin(omega4).*sin(omega3).*cos(inc)...
    -cos(omega4).*sin(omega3)-sin(omega4).*cos(omega3).*cos(inc)...
    sin(omega4).*sin(inc); sin(omega4).*cos(omega3)+cos(omega4).*sin(omega3).*cos(inc)...
    -sin(omega4).*sin(omega3)+cos(omega4).*cos(omega3).*cos(inc)...
    -cos(omega4).*sin(inc); sin(omega3).*sin(inc) cos(omega3).*sin(inc)...
    cos(inc)];

rho = R*rad;

X = rho(1);
Y = rho(2);
Z = rho(3);





   

