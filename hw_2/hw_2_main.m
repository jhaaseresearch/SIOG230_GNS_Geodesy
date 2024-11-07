clear all
close all
clc

%input data

time = 1;
X = 4433469.9438;
Y = 362672.7267;
Z = 4556211.6409;
coord_time_ECEF_example = [time X Y Z];

%number of entries
N = 5;

%preallocate
coord_time_wgs84 = zeros(N,4);
coord_time_ECEF = zeros(N,4);

%generate input matrix
for t = 1:N
    coord_time_ECEF(t,:) = coord_time_ECEF_example;
    coord_time_ECEF(t,1) = t;
end

%check for other coordinates
% Munich: Longitude	11.576124, Latitude	48.137154, Height 520m
% XYZ generated by https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm
x_munich =  4177971;
y_munich =  855800;
z_munich =  4727455;
t = t+1;

%add additional entry to input matrix
coord_time_ECEF(t,:) = [t x_munich y_munich z_munich];

%calculate WGS84 coordinates for all ECEF entries, keeping the time value
coord_time_wgs84 = xyz2wgs(coord_time_ECEF)
