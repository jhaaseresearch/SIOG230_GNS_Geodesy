% Neil Waldhausen
% SIOG 239 Geodesy HW5 (Assignment 2.4.3)
% Compute position of satellites
clc, clear all, close all

% Reading in Data files
eph = read_rinexn("brdc2920.19n");                  % Broadcast Ephemeris Data
rinex_file = 'opmt2920.19o';
[obs,t,gps,apr,hant] = read_rinexo(rinex_file);     % Observation Data

sp3_file = 'igs20756.sp3';
[sp3,sv,excl_sat] = read_sp3(sp3_file);

% Define Constants
c = 0.299792458e9; % Speed of light in meters per second
cutoff = 10; % Elevation cutoff in degrees

% Calculating times in seconds
time_day = datetime(t(:,1)+2000,t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
time_ref = datetime(2019,10,19,0,0,0); % Reference time
time_sec = seconds(time_day - time_ref); % Convert time to seconds of the day
time_week = time_sec + 6*3600*24; % Convert time to seconds of the week

C1 = obs.C1; % Pseudorange C1
sats = width(C1); % Number of satellites
time = length(C1); % Number of epochs

% Initial Values
clock_bias = 1; % Receiver clock bias in nanoseconds

Xo = apr(1,1); % A priori receiver X position
Yo = apr(1,2); % A priori receiver Y position
Zo = apr(1,3); % A priori receiver Z position

% Time of Transmission
ttx = time_week - C1/c; % Time of transmission for each satellite

% Get satellites with available data
isat = ~isnan(C1(1,:)); % Mask for satellites with valid C1 measurements
sat = sv(isat); % Active satellites

% Initialize design matrix A and observation vector L
A = zeros(length(sat), 4);
L = zeros(length(sat), 1);
unknowns = zeros(4, 1); % Unknowns: [dX, dY, dZ, clock_bias]

% Loop over each satellite to build the system of equations
for ic = 1:length(sat)
    % Get satellite position (in ECEF) and clock bias at the time of transmission
    satpos = get_satpos(t, sv(sat(ic)), eph, sat(ic)); 

    % Compute the geometric range from receiver to satellite (rho_0)
    jP = sqrt((X - Xo)^2 + (Y - Yo)^2 + (Z - Zo)^2); 

    % Compute the satellite clock bias (dt), assuming it's given by get_satpos
    satellite_clock_bias = clock_bias * c; % Convert clock bias from seconds to meters

    % Compute the observation (L) - measured pseudorange minus modeled geometric range + clock bias
    L(ic) = C1(ic) - jP + satellite_clock_bias;

    % Variables for the design matrix A
    ax = (X - Xo) / jP; 
    ay = (Y - Yo) / jP;
    az = (Z - Zo) / jP;
    dt = 1; % Partial derivative with respect to clock bias is 1

    % Row of the design matrix A
    A(ic, :) = [ax ay az c];

    % Solve the least squares problem using the pseudo-inverse
    % The unknowns vector contains adjustments to receiver position (dX, dY, dZ) and clock bias
    unknowns = pinv(A) * L;
    
    % Compute the adjusted receiver position (in ECEF) and clock bias
    dX = unknowns(1); 
    dY = unknowns(2);
    dZ = unknowns(3); 
    clock_bias_adjusted = unknowns(4); % Clock bias in meters
    
    % Compute the corrected receiver position (in ECEF)
    X_corr = Xo + dX;
    Y_corr = Yo + dY;
    Z_corr = Zo + dZ;
end

% Display results
fprintf('Corrected Receiver Position (ECEF):\n');
fprintf('X = %f meters\n', X_corr);
fprintf('Y = %f meters\n', Y_corr);
fprintf('Z = %f meters\n', Z_corr);
fprintf('Clock bias (adjusted): %f meters\n', clock_bias_adjusted);

% Compute the covariance matrix for the unknowns
cov_matrix = pinv(A' * A);

% Extract DOP values (GDOP, PDOP, HDOP, VDOP, TDOP)
DOP = sqrt(diag(cov_matrix));

% Display DOP values
fprintf('DOP values (in meters):\n');
fprintf('GDOP = %f\n', DOP(1));
fprintf('PDOP = %f\n', DOP(2));
fprintf('HDOP = %f\n', DOP(3));
fprintf('VDOP = %f\n', DOP(4));
fprintf('TDOP = %f\n', DOP(5));

    