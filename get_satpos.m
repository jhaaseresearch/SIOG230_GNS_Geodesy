function [X,Y,Z] = get_satpos(t1,t2,sv,filename)

    
    %t1 = time step in the given ephemerides file for the specified
         %satillite
    %t2 = transmition time
    %sv = satellite number
    %filename = name of ephemerides file given with .txt exstention

    
    
    %inputs t eph sv
    %t = 600000
    %sv = 1
   

    
    GM  = 3.986005e14;  % Universal gravitational constant times the mass of the Earth, [m^3/s^2]
    gpsPi  = 3.1415926535898;
    omega_e = 7.2921151467E-05;
    %file = readlines('brdc2920.19n.txt');
    file = readlines(filename);

    if sv > 4

    eph = file( [((sv-2) * 104) + (t1-1)*8 + 14] : [ ((sv-2) * 104) + (t1-1) * 8 + 21]);

        else
    
        eph = file( [((sv-1) * 104) + (t1-1)*8 + 14] : [ ((sv-1) * 104) + (t1-1) * 8 + 21]);

    end
 

    Line1 = char(file(14));
    Line2 = char(file(15));
    Line3= char(file(16));
    Line4= char(file(17));
    Line5=char(file(18));
    Line6=char(file(19));
    
    
    day = Line1(10:11);
    hou = Line1(13:14);
    min = Line1(16:17);
    sec = Line1(19:21);
    toe = Line4(4:22);
    roota = Line3(61:79);
    mu_0 = Line2(61:79);
    DeltaN = Line2(42:60);
    e = Line3(23:41);
    C_rs = Line2(23:41);
    C_uc = Line3(4:22);
    omega_0 = Line5(42:60);
    Omega_0 = Line4(42:60);
    Omega_dot = Line5(61:79);
    C_us = Line3(42:60);
    C_rc = Line5(23:41);
    C_ic = Line4(23:41);
    C_is = Line4(61:79);
    i0 = Line5(4:22);
    idot = Line6(4:22);
    
    
    toe = strrep(toe, 'D', 'E');
    roota = strrep(roota, 'D', 'E');
    mu_0 = strrep(mu_0, 'D', 'E');
    DeltaN = strrep(DeltaN, 'D', 'E');
    e = strrep(e, 'D', 'E');
    C_rs = strrep(C_rs, 'D', 'E');
    C_uc = strrep(C_uc, 'D', 'E');
    omega_0 = strrep(omega_0, 'D', 'E');
    Omega_0 = strrep(Omega_0, 'D', 'E');
    Omega_dot = strrep(Omega_dot, 'D', 'E');
    C_us = strrep(C_us, 'D', 'E');
    C_rc = strrep(C_rc, 'D', 'E');
    C_ic = strrep(C_ic, 'D', 'E');
    C_is = strrep(C_is, 'D', 'E');
    i0 = strrep(i0, 'D', 'E');
    idot = strrep(idot, 'D', 'E');
    
    
    toe = str2double(toe);
    roota = str2double(roota);
    mu_0 = str2double(mu_0);
    DeltaN = str2double(DeltaN);
    e = str2double(e);
    
    C_rs = str2double(C_rs);
    C_uc = str2double(C_uc);
    omega_0 = str2double(omega_0);
    Omega_0 = str2double(Omega_0);
    Omega_dot = str2double(Omega_dot);
    C_us = str2double(C_us);
    C_rc = str2double(C_rc);
    C_ic = str2double(C_ic);
    C_is = str2double(C_is);
    i0 = str2double(i0);
    idot = str2double(idot);
    
    day_s = str2double(day) * 86400;
    hour_s = str2double(hou) * 3600;
    min_s = str2double(min) * 60;
    sec = str2double(sec) ;
    
    %t =  sec + min_s + hour_s + toe ;
    
    t_k = t2 - toe;
    
    a = roota^2 ;
    
    
    mu = mu_0 + ( sqrt(GM / a^3) + DeltaN) * t_k;
    
    E  = mu;
    
        for ii = 1:10
            E_old   = E;
            E       = mu + e * sin(E);
            dE      = rem(E - E_old, 2*gpsPi);
    
            if abs(dE) < 1.e-12
                % Necessary precision is reached, exit from the loop
                break;
            end
        end
    
         nu = atan2(sqrt(1 - e^2) * sin(E), cos(E)- e);
    
    phi = nu + omega_0;
    
    
        u = omega_0 + ...
            C_uc * cos(2*phi) + ...
            C_us * sin(2*phi);
        %radius
        r = a * (1 - e*cos(E)) + ...
            C_rc * cos(2*phi) + ...
            C_rs * sin(2*phi);
        %inclination
        i = i0 + idot * t_k + ...
            C_ic * cos(2*phi) + ...
            C_is * sin(2*phi);
            %Compute the angle between the ascending node and the Greenwich meridian
        Omega = Omega_0 + (Omega_dot - omega_e)*t_k - ...
            omega_e * toe;
    
    
    
    
       X = cos(nu)*r * cos(Omega) - sin(nu)*r * cos(i)*sin(Omega);
       Y = cos(nu)*r * sin(Omega) + sin(nu)*r * cos(i)*cos(Omega);
       Z = sin(nu)*r * sin(i);



end





