function [t, lon, lat, h] = xyz2wgs(t, x, y, z)
    a = 6378137;
    f = 1/298.257223563;
    b = a*(1-f);
    e1 = (a^2 - b^2)/(b^2);
    e2 = 2*f - f^2;
    p = sqrt(x.^2 + y.^2);
    theta = atan2(z*a, p*b);
    
    lon = atan2(y,x);
    
    Yinput = z+(sin(theta).^3) *e1*b;
    Xinput = p - (cos(theta).^3)*e2*a;
    lat = atan2(Yinput, Xinput);
    lon = lon*180/pi;
    lat = lat*180/pi;
    h = p.*cos(lat)+z.*sin(lat)-a*sqrt(1-e2*(sin(lat).^2));
end