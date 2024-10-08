% The homework assignment 2 is to write a function xyz2wgs.m to convert from xyz to wgs84 following the equations in the geodesy_lab book section 1.1.2 p. 3.
% please write the function so that the input is a N x 4 array of time(sec), x(m), y(m), z(m) and the output is an array of time, long(deg), lat(deg), height(m).
% Set the time = 1 in your input data.
% To check your code, verify that:
% Lat = 45.8791, Lon = 4.6766, Ele = 432.4222 m
% corresponds to:
% X = 4433469.9438 m, Y = 362672.7267 m, Z = 4556211.6409 m
% 
% To check that your code works for a list of N positions at N different time points repeat the x, y, z values above N times and set the time values to be 1:N


close
clear
clc


x = [ 4433469.9438 4433469.9438 4433469.9438]
y = [362672.7267 362672.7267 362672.7267]
z = [4556211.6409 4556211.6409 4556211.6409]

  
a = 6378137.0000; %meters
f = 1.0/298.257223563;
b = a*(1-f);

e1 = (a^2 - b^2)/(b^2);
e2 = 2*f - f^2;

p = sqrt(x.^2 + y.^2);
theta = atan2(z.*a,p.*b);


Lon = (atan2(y,x)) * (180/pi) ;
Lat = atan2(z+((sin(theta)).^3).*e1.*b, p - (cos(theta)).^3 .* e2 .* a) * (180/pi);
Elev = p.*cos(Lat*(pi/180)) + z.*sin(Lat*(pi/180)) - a.*sqrt(1.0 - e2.*(sin(Lat *(pi/180)).^2));

disp("--------------------------------")

disp(Lat);
disp(Lon);
disp(Elev);




