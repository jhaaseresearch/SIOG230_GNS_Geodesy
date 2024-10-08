function [time_output, long, lat, height] = xyz2wgs(time_input, x,y,z)
    time_output = time_input;
    
    %constants
    a = 6378137.0000;
    f = 1.0/298.257223563;
    b = a * (1 - f );
    e2 = 2*f - f^2;
    e1 = (a^2 - b^2)/b^2;
    p = sqrt(x^2 + y^2);
    theta = atan2(z*a, p*b);
    long = atan2(y,x);
    lat = (atan2(z + (sin(theta))^3 * e1 * b, p - (cos(theta))^3 * e2 * a));
    height = p * cos(lat) + z * sin(lat) - a * sqrt(1.0 - (e2 * sin(lat)^2));
    lat = lat * 180/pi;
    long = long * 180/pi;
end

%