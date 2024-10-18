function [rho] = get_satpos(t, sv, eph)
     
    % choosing satellite number
    svar = eph(1,:);
    idx_prn=find(svar == sv);
    eph=eph(:,idx_prn);
    
    %finding correct 
    toe_ar = eph(18,:);
    [Min, Idx] = min(abs(t-toe_ar));
    % parameters needed
    svprn = eph(1,Idx);
    M0 = eph(3,Idx);
    roota = eph(4,Idx);
    deltan = eph(5,Idx);
    ecc = eph(6,Idx);  
    omega0 = eph(7,Idx);
    cuc = eph(8,Idx);
    cus = eph(9,Idx);
    crc = eph(10,Idx);
    crs = eph(11,Idx);
    cic = eph(14,Idx);
    cis = eph(15,Idx);
    i0 = eph(12,Idx);
    idot = eph(13,Idx);
    Omega0 = eph(16,Idx);
    Omegadot = eph(17,Idx);
    toe = eph(18,Idx);

    

    %computing basic parameters 
    %find correct time of toe
    % t_convert = t * 24 * 3600; % assuming t is in days
    valid_toe = toe(toe <= t);
    closest_toe = max(valid_toe);

    tk = t - closest_toe; %time elapsed since toe 
    GM = 3.986004418 * 10^(14); % m^3 s^-2 
    a = roota.^2 %semi-major axis
    mu = (M0 + sqrt(GM./a.^3) + deltan) * tk;
    
    % iterative solution for E 
    E = mu;  % initial guess 
    tol = 1e-11;  % tolerance for convergence
    max_iter = 100;  % Maximum number of iterations

    for i = 1:max_iter;
        E_new = mu + ecc .* sin(E);  
        if abs(E_new - E) < tol;
            break  
        end
        E = E_new;  
    end

    v = atan2(sqrt(1 - ecc.^2) .* sin(E), cos(E) - ecc);

    % correcting for orbital perturbations
    omega = omega0 + cuc .* cos(2 * (omega0 +v)) + cus .* sin(2 * (omega0 + v)); % argument of perigee
    r = a .* (1 - ecc .* cos(E)) + crc .* cos(2 * (omega0 + v)) + crs .* sin(2 * (omega0 + v)); % radial distance
    incl = i0 + idot * tk + cic .* cos(2 * (omega0 + v)) + cis .* sin(2 * (omega0 + v)); % inclination
 
    % compute the right ascention 
    omegae = 7.2921151467 * 10^-5; % mean angular vel of earth (rad/s)
    
    Omega = Omega0 + (Omegadot - omegae) * tk - omegae*closest_toe;

    % oribtal frame 
    r_orbital = [r .* cos(v); r .* sin(v); 0];
    % fprintf('Size of r_orbital is %d by %d\n', size(r_orbital, 1), size(r_orbital, 2));


    % rotation matrix
    
    R = [cos(Omega).*cos(omega)-sin(Omega).*sin(omega).*cos(incl), -cos(Omega).*cos(omega)-sin(Omega).*sin(omega).*cos(incl), sin(Omega).*sin(incl);
        sin(Omega).*cos(omega)+cos(Omega).*sin(omega).*cos(incl), -sin(Omega).*sin(omega)+cos(Omega).*cos(omega).*cos(incl), -cos(Omega).*sin(incl);
        sin(omega).*sin(incl), cos(omega).*sin(incl), cos(incl)]; 
    % fprintf('Dimensions of R: %d x %d\n', size(R, 1), size(R, 2));

    % applying rotation to get ECEF coords 
   rho = R*r_orbital 

   % X = rho(1);
   % Y = rho(2);
   % Z = rho(3);
   fprintf('Variables t, tk, M0, M, E, v, Omega, omega, r, i, R', t, tk, M0, mu, E, v, Omega, omega, r, R);

end 
%% make function [X,Y,Z]=get_satpos(t, sv, eph) 

% other given parameters 
% clock drift rate - row 2
% clock bias - 19 
% clock drift - 20 
% toc - 21 
% tgd - 22
% trans - 23
