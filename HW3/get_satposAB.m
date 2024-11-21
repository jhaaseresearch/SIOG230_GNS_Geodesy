% This is my function for get satellite position!
%Inputs: the time I wants, the satellite, and the ephemeris:
%Outputs: X,Y,Z of satellite
function [X,Y,Z] = get_satpos(time,sv, eph)
    toe = eph(18,:); % time of eph
    tk = time - toe  ;
% Get useful eph values:
%Establishing variables:
    M0 = eph(3, :); % Mean anomaly in radians
    roota = eph(4, :); % sqrt(semi-major axis)
    a = roota .^2;
    deltan = eph(5, :); % variation of mean angular velocity radians/s
    ecc = eph(6, :); % eccentricity
    omega = eph(7, :); % perigee radians
    cuc = eph(8, :); %correction coef
    cus = eph(9, :); %correction coef
    crc = eph(10, :); %correction coef
    crs = eph(11, :); %correction coef
    i0 = eph(12, :); %inclination, radians
    idot = eph(13, :); %rate of inclination, radians/s
    cic = eph(14, :); %correction coef
    cis = eph(15, :); %correction coef
    Omega0 = eph(16, :); %right ascention radians
    Omegadot = eph(17,:);%rate of right ascention radians/s

%Now setting up to iter. to find E:
    tol = 1e-11; %Tolerance Jennifer Haas used
    Mew =M0 +(sqrt(((3.986004418*10^14))./(roota.^3))+ deltan) .* tk;
    E =  Mew;
    maxiter=100
    for kg = 1 : maxiter
        Efut = Mew + ecc.*sin(E);
        if max(Efut - E) < tol
            break
        end
        E = Efut;
    end


%True anomaly v:
    v = atan((sqrt(1-(ecc.^2)).* sin(E))./(cos(E)-ecc));
%Orbital Pertubation Corrections:
    omegacorrect = omega + cuc.*cos(2*(omega+v)) +cus .*sin(2*(omega+v));
    r = (a.*(1-ecc.*cos(E))) + (crc.*cos(2*(omega+v))) + (crs .*sin(2*(omega+v)));
    i = i0 + idot.*tk + cic.*cos(2*(omega+v)) +cis .*sin(2*(omega+v));

% Right Ascension
    OMEGA = Omega0 + ((Omegadot - 7.2921151467*1e-5).*tk) - (7.2921151467*1e-5 .* toe);

% Converting Satellite position to ECEF frame:
    rvec = zeros(3, length(r));
    rvec(1,:) = (r .* cos(v));
    rvec(2,:)= r .* sin(v);

    R11 =(cos(OMEGA).*cos(omegacorrect))-(sin(OMEGA).*sin(omegacorrect).*cos(i));
    R12 = (-cos(OMEGA).*sin(omegacorrect))-(sin(OMEGA).*cos(omegacorrect).*cos(i));
    R13 = sin(OMEGA) .* sin(i);
    R21 = (sin(OMEGA).*cos(omegacorrect))+(cos(OMEGA).*sin(omegacorrect).*cos(i));
    R22 = (-sin(OMEGA).*sin(omegacorrect))+(cos(OMEGA).*cos(omegacorrect).*cos(i));
    R23 = -cos(OMEGA).*sin(i);
    R31 = sin(omegacorrect).*sin(i);
    R32 = cos(omegacorrect).*sin(i);
    R33 = cos(i);
    for kkh = 1: length(r)
        R = [R11(kkh) R12(kkh) R13(kkh)
            R21(kkh) R22(kkh) R23(kkh)
            R31(kkh) R32(kkh) R33(kkh)];
        disp(size(R))
        rnow = rvec(1:3, kkh)
        rhoe = R .*rnow;
       X = rhoe(1);
       Y = rhoe(2);
       Z = rhoe(3);
    end
end