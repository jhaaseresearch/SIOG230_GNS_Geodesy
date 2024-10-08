clear all
close all
clc

%input data

time = 1;
X = 4433469.9438;
Y = 362672.7267;
Z = 4556211.6409;
N = 4;
time_output = zeros(N,1);
lon = zeros(N,1);
lat = zeros(N,1);
h = zeros(N,1);

for t = 1:N
    [time_output(t), lon(t), lat(t), h(t)] = xyz2wgs(t, X,Y,Z);
end

output = [time_output, lon, lat, h]

%one single run
%[time_output, lon, lat, h] = xyz2wgs(time, X,Y,Z);
% fprintf('time output [s]:\n %g\n', time_output);
% fprintf('longitude [deg] output:\n %g\n', lon);
% fprintf('latitude [deg]:\n %g\n', lat);
% fprintf('height [m]:\n %g\n', h);