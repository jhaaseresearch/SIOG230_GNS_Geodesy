function [time,long,lat,height] = xyz2wgs(input)
    times = input(:,1);
    a = 6378137.0000;
    f = 1.0/298.257223563;
    b = a * (1 -f);
    e2 = (2*f) - f^2 ;
    e1 = (a^2 - b^2)/(b^2);
    X = input(:,2)
    Y = input(:,3)
    Z = input(:,4);
    p = sqrt(X.^2 + Y.^2);
    sigma = atan2(Z.*a, p*b);

    long = (180/pi)* atan2(Y,X)
    lat = (180/pi) * atan2(Z + (sin(sigma).^3) * e1 * b, p - (cos(sigma).^3)* e2 *a)
    time = times;
    height = p .* cos(lat * (pi/180)) + Z .* sin(lat*(pi/180)) - a * sqrt(1.0 -e2*sin(lat*(pi/180)).^2)

end