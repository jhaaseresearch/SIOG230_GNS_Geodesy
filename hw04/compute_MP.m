function [MP1, MP2, t_hrs] = compute_MP(obs, t, gps, sv, f1, f2, range)

% READ_RINEXO Compute MP1, MP2
%
% Input: obs = obs output from read_rinexo.m
%        t = t output from read_rinexo.m
%        gps = gps output from read_rinexo.m
%        sv = PRN number of interest
%        f1 = frequency of L1 [Hz = cycles/period]
%        f2 = frequency of L2
%        range = range of epochs over which to calculate MP1, MP2
%               
% Output: MP1, MP2 = arrays
%         t_hrs = time in hours from obsfile start
%

t_hrs = t(range,4) + t(range,5)/60 + t(range,6)/3600;
% Find MP for PRN and range input
idx_prn = find(gps == sv);

% Phase data (L1, L2) is in units of cycles of carrier phase
L1 = obs.L1(range,idx_prn);
L2 = obs.L2(range,idx_prn);
% Pseudorange data in meters
C1 = obs.C1(range,idx_prn);
P1 = obs.P1(range,idx_prn);
P2 = obs.P2(range,idx_prn);
c = 0.299792458*1e9;
lambda1 = c/f1;
lambda2 = c/f2;
alpha = (f1/f2)^2;

% L1, L2 --> meters
L1 = L1*lambda1;
L2 = L2*lambda2;

MP1 = P1 - (2/(alpha - 1) + 1) .* L1 + (2/(alpha - 1)) .* L2;
MP2 = P2 - (2*alpha/(alpha-1)).*L1 + (2*alpha/(alpha-1) - 1) .* L2;


% Remove mean from MP1, MP2
MP1 = MP1 - mean(MP1(~isnan(MP1)));
MP2 = MP2 - mean(MP2(~isnan(MP2)));


