function plot_MP1_MP2(obs, t, prn, prn_target, hant, f1, f2, c, arc_duration)
    % Calculate wavelengths
    lambda1 = c / f1; % Wavelength corresponding to f1
    lambda2 = c / f2; % Wavelength corresponding to f2
    
    % Compute ionospheric delay ratio
    alpha = (f1 / f2)^2;
    
    % Find the index of the target satellite (PRN)
    index_target = find(prn == prn_target); 
    
    % Extract observations for the target satellite
    % Phase measurement (cycles)
    L1 = obs.L1(:, index_target); 
    L2 = obs.L2(:, index_target); 
    % Pseudorange measurement (meters)
    P1 = obs.P1(:, index_target); 
    P2 = obs.P2(:, index_target); 
    
    % Compute MP1
    MP1 = P1 - ((2 / (alpha - 1)) + 1) * (lambda1 * L1) + (2 / (alpha - 1)) * (lambda2 * L2);
    % Compute MP2
    MP2 = P2 - (2 * alpha / (alpha - 1)) * (lambda1 * L1) + (2 * alpha / (alpha - 1) - 1) * (lambda2 * L2);
%% Calculate and plot MP1 and MP2 with arcs
    % Assuming 30-second intervals, compute samples per arc
    samples_per_arc = 8; % arc_duration hours
    
    % Total number of arcs
    num_arcs = floor(length(t) / samples_per_arc);
    
    % Preallocate for arcs plotting
    MP1_arc_all = nan(length(MP1), 1); % Fill with NaNs initially
    MP2_arc_all = nan(length(MP2), 1); % Fill with NaNs initially
    
    % Loop over each arc
    for arc = 1:num_arcs
        % Define the start and end of the arc
        arc_start = (arc - 1) * samples_per_arc + 1;
        arc_end = min(arc * samples_per_arc, length(MP1)); % Handle last partial arc if exists
        
        % Extract MP1 and MP2 for the current arc
        MP1_arc = MP1(arc_start:arc_end);
        MP2_arc = MP2(arc_start:arc_end);
        
        % Compute the mean for the current arc
        mean_MP1_arc = mean(MP1_arc(~isnan(MP1_arc)));
        mean_MP2_arc = mean(MP2_arc(~isnan(MP2_arc)));
        
        % Subtract the mean from the current arc
        MP1_arc_all(arc_start:arc_end) = MP1_arc - mean_MP1_arc;
        MP2_arc_all(arc_start:arc_end) = MP2_arc - mean_MP2_arc;
    end
    % Time in hours
    time_hours = ((0:length(t) - 1) * 30) / 3600;
    % Plot the results with arcs
    figure;
    subplot(1, 2, 1);
    plot(time_hours, MP1_arc_all, 'b-', 'LineWidth', 1.5);
    xlabel('Time (hours)');
    ylabel('MP1 [m]');
    title('MP1 vs Time (with arcs)');
    grid on;

    subplot(1, 2, 2);
    plot(time_hours, MP2_arc_all, 'r-', 'LineWidth', 1.5);
    xlabel('Time (hours)');
    ylabel('MP2 [m]');
    title('MP2 vs Time (with arcs)');
    grid on;
    
    % Title for the arc-based subplots
    sgtitle(['MP1 and MP2 with ', num2str(arc_duration), '-hour arcs for PRN ', num2str(prn_target)]);
end
