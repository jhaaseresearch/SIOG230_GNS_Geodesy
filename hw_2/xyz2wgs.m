function coord_wgs = xyz2wgs(coord_ECEF)
    % this function generates coordinate containing timestamp (s),
    % longitude (deg), latitude (deg) and height (m) in WGS84 for
    % a given coordinate in the format timestamp (s), x (m), y (m), z(m)
    % ECEF

    time_input = coord_ECEF(:,1);
    x = coord_ECEF(:,2);
    y = coord_ECEF(:,3);
    z = coord_ECEF(:,4);
    time_output = time_input;
    
    %constants
    a = 6378137.0000;
    f = 1.0/298.257223563;
    b = a .* (1 - f );
    e2 = 2.*f - f.^2;
    e1 = (a.^2 - b.^2)./b.^2;
    p = sqrt(x.^2 + y.^2);
    theta = atan2(z.*a, p.*b);
    long = atan2(y,x);
    lat = (atan2(z + (sin(theta)).^3 .* e1 .* b, p - (cos(theta)).^3 .* e2 .* a));
    height = p .* cos(lat) + z .* sin(lat) - a .* sqrt(1.0 - (e2 .* sin(lat).^2));
    lat = lat .* 180/pi;
    long = long .* 180/pi;
    coord_wgs = [time_output, long, lat, height];
end

%[time_output, long, lat, height]
%time_input, x,y,z