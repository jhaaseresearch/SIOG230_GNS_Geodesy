% APPROX POSITION XYZ
c = 0.299792458e9; % m/s
[obs,t,prn,apr,hant] = read_rinexo('opmt2920.19o');  % Read RINEX data
C1 = obs.C1; 
[eph, alpha, beta, dow] = read_rinexn('brdc2920.19n'); % Read ephemeris data
t_sec = (0:length(t) - 1) * 30;  % Time in seconds from epoch
% A priori position of the receiver
x0 = apr(1);
y0 = apr(2);
z0 = apr(3);
%% 
% Number of epochs (length of the observation data)
num_epochs = length(t_sec);
% Loop over each epoch
for epoch_index = 1:num_epochs
    disp(epoch_index)
    % Receiver time for the current epoch
    t = dow * 24 * 3600 + t_sec(epoch_index);

    % Find satellites with valid C1 data for the current epoch
    nonNaN_satelites = find(~isnan(obs.C1(epoch_index, :)));
    
    % Initialize arrays to store results for the current epoch
    rho_0_all = zeros(1, length(nonNaN_satelites));  % Distance for each epoch and satellite
    L_all = zeros(1, length(nonNaN_satelites));  % Observed minus computed distance
    H_all = zeros(length(nonNaN_satelites), 4);  % Design matrix for this epoch
    
    % Loop over each satellite for the current epoch
    for i = 1:length(nonNaN_satelites)
        sv = prn(nonNaN_satelites(i)); % Current satellite PRN number

        % Transmission time calculation
        t_tx = t - obs.C1(epoch_index, nonNaN_satelites(i)) / c;

        % Get satellite position and clock bias
        satpos = get_satpos(t_tx, sv, eph, 3);
        xs = satpos(1);
        ys = satpos(2);
        zs = satpos(3);
        dt = satpos(4);

        % Compute the modeled observables
        rho_0 = sqrt((xs - x0)^2 + (ys - y0)^2 + (zs - z0)^2);
        L = obs.C1(epoch_index, nonNaN_satelites(i)) - rho_0 + c * dt;

        % Store computed values
        rho_0_all(i) = rho_0;
        L_all(i) = L;

        % Calculate partial derivatives for the design matrix H
        H_all(i, 1) = (x0 - xs) / rho_0;  % d(rho)/dx
        H_all(i, 2) = (y0 - ys) / rho_0;  % d(rho)/dy
        H_all(i, 3) = (z0 - zs) / rho_0;  % d(rho)/dz
        H_all(i, 4) = -c * 1e-9;  % Partial derivative of the clock bias term
    end
    all_x_positions = [];
    all_y_positions = [];
    all_z_positions = [];
    % Solve for position and clock bias for this epoch using least squares
    % Construct the design matrix H and observed vector L for valid satellites
    valid_satellites = find(~isnan(L_all));  % Indices of valid satellites
    
    % Extract the H matrix and L vector for valid satellites
    H = H_all(valid_satellites, :);  
    L = L_all(valid_satellites)';  % Row vector for L
    
    % Solve for the unknowns (position corrections and clock bias)
    X = lscov(H, L);  % Solve least squares
    
    % Update receiver's position and clock bias
    x0 = x0 + X(1);
    y0 = y0 + X(2);
    z0 = z0 + X(3);
    delta_bias = X(4) / 1000;  % Convert clock bias adjustment to seconds
    clock_biases(epoch_index) = delta_bias;

    % Store the updated positions for this epoch
    x_positions(epoch_index) = x0;
    y_positions(epoch_index) = y0;
    z_positions(epoch_index) = z0;
    
    % Store the updated positions for this epoch
    all_x_positions = [all_x_positions, x_positions];
    all_y_positions = [all_y_positions, y_positions];
    all_z_positions = [all_z_positions, z_positions];
    % Compute covariance matrix for the final receiver position
    Cov_X = inv(H' * H);  % Covariance matrix (simplified)
    
    % Construct the rotation matrix R (assumed constant over epochs)
    input = [x_positions(end), y_positions(end), z_positions(end), clock_biases(end)];
    
    output = xyz2wgs(input);
    lambda = output(2);  % Longitude (degrees)
    phi = output(3);     % Latitude (degrees)
    R = [
        -sin(phi) * cos(lambda), -sin(phi) * sin(lambda), cos(phi);
        -sin(lambda), cos(lambda), 0;
        cos(phi) * cos(lambda), cos(phi) * sin(lambda), sin(phi)
    ];
    
    % Compute the matrix CL = R * Cov_X * R'
    CL = R * Cov_X(1:3, 1:3) * R';  % Only position-related covariance
    
    % Extract variances and covariances
    sigma_n2 = CL(1, 1);  % Variance in the north direction
    sigma_e2 = CL(2, 2);  % Variance in the east direction
    sigma_u2 = CL(3, 3);  % Variance in the up direction
    sigma_t = Cov_X(4, 4);  % Time uncertainty
    
    % Compute the DOP factors for this epoch
    VDOP = sqrt(sigma_u2);  % Vertical DOP
    HDOP = sqrt(sigma_n2 + sigma_e2);  % Horizontal DOP
    PDOP = sqrt(sigma_n2 + sigma_e2 + sigma_u2);  % Position DOP
    TDOP = sigma_t;  % Time DOP
    GDOP = sqrt(sigma_n2 + sigma_e2 + sigma_u2 + sigma_t^2);  % Geometric DOP
    
    % Display results for this epoch
    fprintf('Epoch %d:\n', epoch_index);
    fprintf('Receiver Position: x = %.4f m, y = %.4f m, z = %.4f m\n', x0, y0, z0);
%     fprintf('VDOP: %.4f\n', VDOP);
%     fprintf('HDOP: %.4f\n', HDOP);
%     fprintf('PDOP: %.4f\n', PDOP);
%     fprintf('TDOP: %.4f\n', TDOP);
%     fprintf('GDOP: %.4f\n', GDOP);
% Display receiver and satellite positions for this epoch
end
figure;
hold on;

% Plot receiver position
plot3(all_x_positions, all_y_positions, all_z_positions, 'bo', 'MarkerSize', 5, 'LineWidth', 3);

title('Receiver Position');
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
legend('Receiver');
grid on;
hold off;
