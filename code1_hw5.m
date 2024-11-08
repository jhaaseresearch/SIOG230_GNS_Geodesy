%
%
%
eph_all =read_rinexn("brdc2920.19n")
[obs,t,gps,apr,hant] = read_rinexo('opmt2920.19o')
% apr are the given coordinates of the satellites
%
svprn = eph_all(1,:);
clock_drift_rate = eph_all(2,:);
M0 = eph_all(3,:);
roota = eph_all(4,:);
deltan = eph_all(5,:);
ecc = eph_all(6,:);
omega = eph_all(7,:);
cuc = eph_all(8,:);
cus = eph_all(9,:);
crc = eph_all(10,:);
crs = eph_all(11,:);
i0 = eph_all(12,:);
idot = eph_all(13,:);
cic = eph_all(14,:);
cis = eph_all(15,:);
Omega0 = eph_all(16,:);
Omegadot = eph_all(17,:);
toe = eph_all(18,:);
clock_bias = eph_all(19,:);
clock_drift = eph_all(20,:);
toc = eph_all(21,:);
tgd = eph_all(22,:);
trans = eph_all(23,:);
%
obs_c1 = obs.C1
t0 = 518400 % time in seconds
c = 0.299792458*10^(9)% m/s

validColumnsFirstRow = ~isnan(obs_c1(1, :));
validColumnIndices = find(validColumnsFirstRow);
obs_c1_index_Rvalues = validColumnIndices'


t_D = datetime(2000+t(:,1),t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
t_seconds = t_D-t_D(1)
first_t_gpsTime = 518400
last_t_gpsTime = 604770
time_steps_30s = first_t_gpsTime:30:last_t_gpsTime
% First coordinates from rinex file
x0 = apr(1); y0 = apr(2); z0 = apr(3);
l = length(time_steps_30s)
% Loop for every time
for i=1:l
    x0 = x0 ; y0=y0 ; z0=z0

    validColumnsFirstRow = ~isnan(obs_c1(i, :));
    % Get the indices of the columns where the first row is not NaN
    validColumnIndices = find(validColumnsFirstRow);
    obs_c1_index_Rvalues = validColumnIndices'
    L_vector = ones(length(obs_c1_index_Rvalues), 1);
    A_matrix = ones(length(obs_c1_index_Rvalues), 4);
    k = length(validColumnIndices)
    solutions = ones(length(validColumnIndices), 4);
    
    % Loop over each satellite and store the solution
    for j=1:k
        % Get the solution for satellite sv at time t
        satellite_id = validColumnIndices(j);
        t_cor = (time_steps_30s(i) - (obs_c1_index_Rvalues/c));
        coord = get_satpos(t_cor(i),satellite_id,eph_all,0)
        solutions(j, :) = coord(1:4);
        x_sapos = solutions(:,1);
        y_sapos = solutions(:,2);
        z_sapos = solutions(:,3);
        dt_sapos = solutions(:,4);
        x0 = x0 .* (ones(length(x_sapos),1));
        y0 = y0 .* (ones(length(y_sapos),1));
        z0 = z0 .* (ones(length(z_sapos),1));
        c_vector = c .* (ones(length(x_sapos),1))
        rho_sat = sqrt((x0-x_sapos).^2 + (y0-y_sapos).^2 + (z0-z_sapos).^2)
        
        % Constructing L vector
        L_vector = [obs_c1_index_Rvalues - rho_sat + (-c_vector.*dt_sapos)]
        
  
        % Elements of A matrix
        ax = (x0-x_sapos) .* (1./rho_sat);
        ay = (y0-y_sapos) .* (1./rho_sat);
        az = (z0-z_sapos) .* (1./rho_sat);
        
        A_matrix = [ax ay az c_vector];
         
        x_sol = pinv(A_matrix)*L_vector
    end

        x0 = x_sol(1);
        y0 = x_sol(2);
        z0 = x_sol(3);
        dt = x_sol(4);
end

