%%   SIOG 239 - Homework 2    %%
format longG
labelm = 'm';
labeldeg = 'deg';

% array1 = [t; x; y; z];
array1 = [1.0; 4433469.9438; 362672.7267; 4556211.6409];

% Input Variables in x,y,z coordinates
x = array1(2, 1);   % meters
y = array1(3, 1);    % meters
z = array1(4, 1);   % meters
t = array1(1, 1);

disp('XYZ Coordinates')
disp(array1)

% Converstion Constants
a = 6378137.0000;
f = 1.0/298.257223563;
b = a*(1-f);
e2 = (2*f)-(f^2);
e1 = (a^2-b^2)/(b^2);

% Auxillary Quantites
p = sqrt(x^2+y^2);
theta = atan2(z*a, p*b);

i = z+(sin(theta))^3*e1*b;
j = p-(cos(theta))^3*e2*a;


%Solving for Lat/Long from x, y, z
long = atan2(y, x);         % radians
lat = atan2(i, j);          % radians
longdeg = long*(180/pi);    % degrees
latdeg = lat*(180/pi);      % degrees
height = p*cos(lat)+z*sin(lat)-a*sqrt(1.0-e2*sin(lat)^2);   % meters

array2 = [t; longdeg; latdeg; height];
disp('Lat/Long Coordinates')
disp(array2)





