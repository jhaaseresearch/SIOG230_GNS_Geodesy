% Code which tests HW2 function xyz2wgs 
% Below code confirms that xyz2wgs.m is working properly
% Author: Tommy Stone
% Date: 10/06/24
% Quarter: Fall 2024

clear all; close all; clc;

% Write initial values for testing 
x = 4433469.9438;
y = 362672.7267;
z = 4556211.6409;
t = 1;
A = [t,x,y,z];
B = xyz2wgs(A);

% return test values
lat = 45.8791;
lon = 4.6766;
elevation = 432.4222;

% printing out comparisons
sprintf("x = " + num2str(x) + "m, Longitude from xyz2wgs " + num2str(B(2)) + " degrees")
sprintf("y = " + num2str(y) + "m, Latitude from xyz2wgs " + num2str(B(3)) + " degrees")
sprintf("elevation = " + num2str(z) + " m, Height from xyz2wgs " + num2str(B(4)) + " m")

% creating a new matrix of values
N = 30;
t2 = 1:N';
x2 = linspace(4433469.9438,4500000,N)';
y2 = linspace(y,y + 10000,N)';
z2 = linspace(z, z + 1000,N)';

% Assigning to matrix
A2 = nan(N,4);
A2(:,1) = t2;
A2(:,2) = x2;
A2(:,3) = y2;
A2(:,4) = z2;

% confirming this works with a matrix
B2 = xyz2wgs(A2);





