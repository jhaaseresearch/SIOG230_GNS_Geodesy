% Example input
% For a single time step
t = 1; % Time step in seconds
x = 4433469.9438; % X coordinate in meters
y = 362672.7267; % Y coordinate in meters
z = 4556211.6409; % Z coordinate in meters
input = [t, x, y, z];
output = xyz2wgs(input);
%% 
% For multiple time steps
N = 10;
t = (1:N)'; % Time step in seconds (column vector)
x = 4433469.9438 * ones(N, 1); % X vector, repeated N times
y = 362672.7267 * ones(N, 1);  % Y vector, repeated N times
z = 4556211.6409 * ones(N, 1); % Z vector, repeated N times
input = [t, x, y, z];
output = xyz2wgs(input);