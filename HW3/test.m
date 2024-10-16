%1
[eph, alpha, beta, dow] = read_rinexn('Data/epgga9.292');  % Load ephemeris data
%Check the program
%At the first epoch of day 292/2019
t = dow * 24 * 3600;  % Convert time to seconds
sv = 31;  % Satellite number

% Call the function to compute satellite position
disp(['Satellite position at t = ', num2str(t), ' seconds: ']);
rho_e = get_satpos(t, sv, eph);
disp(['Satellite position in ECEF (XYZ): ', num2str(rho_e(1)), ' meters, ', ...       
        num2str(rho_e(2)), ' meters, ', num2str(rho_e(3)), ' meters']);
%% 
[eph, alpha, beta, dow] = read_rinexn('Data/epgga9.292');  % Load ephemeris data
% 15 minute before the end of day 292/2019,
t = dow * 24 * 3600 + 24 * 3600 - 15 * 60;  % Convert time to seconds
sv = 31;  % Satellite number

% Call the function to compute satellite position
disp(['Satellite position at t = ', num2str(t), ' seconds: ']);
rho_e = get_satpos(t, sv, eph);
disp(['Satellite position in ECEF (XYZ): ', num2str(rho_e(1)), ' meters, ', ...       
        num2str(rho_e(2)), ' meters, ', num2str(rho_e(3)), ' meters']);
%% 
% 2
[eph, alpha, beta, dow] = read_rinexn('Data/epgga9.292');  % Load ephemeris data
%Compute ECEF coordinates for satellite 31 every 15 minutes
% Interval of 15 minutes in seconds
interval = 15 * 60;  % 900 seconds

% Day of the week (dow) and corresponding time range
start_time = dow * 24 * 3600;         % Start of the day in seconds
end_time = dow * 24 * 3600 + 24 * 3600 - 15 * 60;     % End of the day in seconds

% Initialize time vector
t = start_time:interval:end_time;  % Time in seconds every 15 minutes

% Satellite number
sv = 31;

% Initialize matrix to store satellite positions for each time step
satellite_positions = zeros(length(t), 4);  % Preallocate for speed

% Iterate over each time step in 't'
time_passed = 0;
for i = 1:length(t)
    current_time = t(i);  % Current time in seconds
    disp(['Satellite position at t = ', num2str(current_time), ' seconds: ']);
    disp(['Time (seconds since start of the day): ', num2str(time_passed)]);
    % Compute satellite position at the current time
    rho_e = get_satpos(current_time, sv, eph);
    
    % Store the satellite position in the matrix at time since the start of the
    % day
    satellite_positions(i, :) = [time_passed, rho_e'];
    time_passed = time_passed + 15 * 60;
    % Display the current time and satellite position in ECEF (XYZ)
    % coordinates (meters)
    disp(['Satellite position in ECEF (XYZ): ', num2str(rho_e(1)), ' meters, ', ...       
        num2str(rho_e(2)), ' meters, ', num2str(rho_e(3)), ' meters']);
end
%% 
%compare with the coordinates given in the corresponding sp3 file (precise IGS orbits).
% extracts the XYZ position of a given satellite sv from a .sp3 file 
[sp3,sv_all,excl_sat] = read_sp3('Data/igs20756.sp3');
% Check if the satellite exists in the SP3 data
sv_field = ['prn' num2str(sv)];
if isfield(sp3, sv_field)
    sp3_sv_data = getfield(sp3, sv_field);  % Extract data for satellite 31
else
    error(['Satellite ' num2str(sv) ' not found in SP3 data.']);
end

% Extract SP3 positions for satellite 31
satellite_id = 31;
sp3_field_name = sprintf('prn%d', satellite_id);  % Construct field name
sp3_positions = sp3.(sp3_field_name);  % Get SP3 positions [T, X, Y, Z, dT]

% Initialize array for residuals
num_positions = size(sp3_positions, 1);
residuals = zeros(num_positions, 1);  % Preallocate for residual norms

% Compute residuals
for i = 1:num_positions
    % Time from SP3
    T_sp3 = sp3_positions(i, 1);  % Time (seconds of the day)
    
    % Find corresponding broadcast position
    % Locate the index in satellite_positions that matches T_sp3
    [~, index] = min(abs(satellite_positions(:, 1) - T_sp3));
    
    % Extract broadcast position
    broadcast_position = satellite_positions(index, 2:4);  % Get [X, Y, Z]
    
    % Extract SP3 position
    sp3_position = sp3_positions(i, 2:4);  % Get [X, Y, Z]
    
    % Compute the 3-D residual
    residuals(i) = sqrt( ...
        (sp3_position(1) - broadcast_position(1))^2 + ...
        (sp3_position(2) - broadcast_position(2))^2 + ...
        (sp3_position(3) - broadcast_position(3))^2 ...
    );
end

% Create time vector for plotting
time_vector = sp3_positions(:, 1);  % Using SP3 time for x-axis

% Plotting the norm of the residuals
figure;
plot(time_vector/3600, residuals, '-o');
title('3-D Residual Position Norm for Satellite 31');
xlabel('Time (Hours since start of day)');
ylabel('Residual Norm (meters)');
grid on;
% Save the figure
saveas(gcf, 'residuals_plot.png'); 
%% 
% 3
% the orbit difference between  sp3 and broadcast XYZ positionis is the lowest at the time of toe, every 2 hours, 
% and due to luck of data and interpolation the error increases with time
% until the next toe 2 hours later.

