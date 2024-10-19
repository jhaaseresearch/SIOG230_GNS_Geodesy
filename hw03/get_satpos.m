function [X,Y,Z] = get_satpos(t, sv, eph)

% get_satpos	transform broadcast orbit -> ECEF
%
% Input: t  = time of interest (secs of current GPS week)
%        sv = PRN of satellite
%        eph = ephimeredes output by read_rinexn.m
%
% Output: X, Y, Z

%%% Choose correct satellite PRN number
sv_arr = eph(1,:);
idx_prn = find(sv_arr == sv);
eph = eph(:,idx_prn);

% Get tk = time elapsed since toe
toe_arr = eph(18,:);
% Need to find correct toe = time of ephemeris
% toe should be as close to input t as possible
[~, Idx] = min(abs(t - toe_arr));
toe = eph(18,Idx);
tk = t - toe;

% Read in variables at our chosen index = Idx
svprn = eph(1,Idx);
clock_drift_rate = eph(2,Idx);
M0 = eph(3,Idx);
roota = eph(4,Idx);
deltan = eph(5,Idx);
ecc = eph(6,Idx);
omega0 = eph(7,Idx);
cuc = eph(8,Idx);
cus = eph(9,Idx);
crc = eph(10,Idx);
crs = eph(11,Idx);
i0 = eph(12,Idx);
idot = eph(13,Idx);
cic = eph(14,Idx);
cis = eph(15,Idx);
Omega0 = eph(16,Idx);
Omegadot = eph(17,Idx);
toe = eph(18,Idx);
clock_bias = eph(19,Idx);
clock_drift = eph(20,Idx);
toc = eph(21,Idx);
tgd = eph(22,Idx);
trans = eph(23,Idx);


% Compute Mean Anomaly mu at t
GM = 3.986004418*10^14;
mu = M0 + (sqrt(GM)./(roota.^3) + deltan) .* tk;


% % Compute Eccentric Anomaly E %%% WRONG - WHY?
% tol = 1e-11;
% E0 = mu;
% E = mu + ecc .* sin(E0);
% 
% while max(E - E0) > tol
%     E0 = E;
%     E = mu + ecc .* sin(E0);
% end

%% Iteratively solve for E (eccentric anomaly) with a convergence tolerance of 0.001
E = mu;  % Start with E = µ
for iter = 1:100000
    E_new = mu + ecc * sin(E);  % Update equation
    if abs(E_new - E) < 1e-11  % Convergence tolerance
        break;
    end
    E = E_new;  % Update E
end


% Compute True Anomaly v
v = atan2(sqrt(1-ecc.^2).*sin(E), cos(E)-ecc);

% Correct for orbital perturbations
% Argument of perigee:
omega = omega0 ...
    + cuc .* cos(2*(omega0 + v)) ...
    + cus .* sin(2*(omega0 + v));

% Radial Distance:
r = roota.^2 .* (1-ecc.*cos(E)) ...
    + crc .* cos(2*(omega0 +v)) ...
    + crs .* sin(2*(omega0 + v));

% Inclination:
inc = i0 + idot .* tk ...
    + cic .* cos(2*(omega0 + v)) ...
    + cis .* sin(2*(omega0 + v));

% Compute right ascension
omegae = 7.2921151467 * 1e-5; % mean angular vel of Earth, rad/s

Omega = Omega0 + (Omegadot - omegae)*tk - omegae*toe;

% Convert satellite position from orbital plane -> ECEF coords
% Satellite position in orbital frame:
rvec = [r.*cos(v); r.*sin(v); 0];
% Rot matrix orbital frame -> ECEF frame:
R = [cos(Omega).*cos(omega) - sin(Omega).*sin(omega).*cos(inc), ...
    -cos(Omega).*sin(omega) - sin(Omega).*cos(omega).*cos(inc), ...
    sin(Omega).*sin(inc); ...
    sin(Omega).*cos(omega) + cos(Omega).*sin(omega).*cos(inc), ...
    -sin(Omega).*sin(omega) + cos(Omega).*cos(omega).*cos(inc), ...
    -cos(Omega).*sin(inc); ...
    sin(omega).*sin(inc), ...
    cos(omega).*sin(inc), ...
    cos(inc)];

% Perform Rotation:
rhoe = R*rvec;
X = rhoe(1);
Y = rhoe(2);
Z = rhoe(3);