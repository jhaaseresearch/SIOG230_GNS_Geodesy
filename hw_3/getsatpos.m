function [rho] = getsatpos(t,sv, eph)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    

    %constants:
    GM = 3.98600441831014;% m^3/(s^2);
    omegaE = 7.2921151467e-5;% rad/s

    %papameters from navigation file:
    %closest toe:
    toe = eph(18,:);
    t = t*3600*24;
    del_toe_t = abs(toe-t);
    [min_val, idx] = min(del_toe_t);
    toe_closest = toe(idx);
    t_k = t - toe_closest;
    
    svprn = eph(1,:);
    M0 = eph(3,idx);
    roota = eph(4,idx);
    deltan = eph(5,idx);
    ecc = eph(6,idx);  
    omega0 = eph(7,idx);
    cuc = eph(8,idx);
    cus = eph(9,idx);
    crc = eph(10,idx);
    crs = eph(11,idx);
    cic = eph(14,idx);
    cis = eph(15,idx);
    i0 = eph(12,idx);
    idot = eph(13,idx);
    Omega0 = eph(16,idx);
    Omegadot = eph(17,idx);
    a = roota.^2;
    
   
    %mean anomaly at t
    mu = M0 + (sqrt(GM)./(roota.^3) + deltan) .* t_k;
    E = mu;
    tol = 1e-11;
    for k = 1:20
        E_new = mu + ecc * sin(E);
        if(abs(E-E_new) <= tol)
            break;
        else
            E = E_new;
        end
    end
    
    %true anomaly v
    v = atan2((sqrt(1-ecc.^2).*sin(E)), (cos(E)-ecc));%maybe have to flip y and x!!
    omega = omega0 + cuc.*cos(2.*(omega0 + v)) + cus.*sin(2.*(omega0 + v));
    r = a.*(1 - ecc.*cos(E)) + crc.* cos(2.*(omega0 +v) + crs .* sin(2.*(omega0 +v)));
    i = i0 + idot.* t_k + cic.*cos(2.*(omega0 + v)) + cis.*sin(2.*(omega0 + v));
    

    Omega = Omega0 + (Omegadot - omegaE).*t_k - (omegaE.*toe_closest);

    r_coord = [r.*cos(v); r.*sin(v); 0];
    R = [(cos(Omega).*cos(omega) - sin(Omega).*sin(omega).*cos(i))...
        (-cos(Omega).*sin(omega)-sin(Omega).*cos(omega).*cos(i))...
        (sin(Omega).*sin(i));...
        (sin(Omega).*cos(omega)+cos(Omega).*sin(omega).*cos(i))...
        (-sin(Omega).*sin(omega)+cos(Omega).*cos(omega).*cos(i))...
        (-cos(Omega).*sin(i));...
        (sin(omega).*sin(i)) (cos(omega).*sin(i)) (cos(i))];
    
    rho = R.*r;
    yuw=3+3;


end