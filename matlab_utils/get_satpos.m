function satpos = get_satpos(t,sv,eph,flag);

% SATPOS   Computes X,Y,Z ECEF coordinates from ephemerides data.
%          Also returns the satellite clock bias dt.
%
%     Input:
%       - t = time of data
%             **MUST** be given in seconds of the current GPS week
%       - sv = satellite PRN number
%       - eph = ephemerides matrix created from a rinex navigation
%               file using read_rinexn.m
%       - flag:
%         0 --> use closest toe
%         1 --> use last toe (as GPS receiver would do)
%         2 --> use first toe of day for the whole session
%         3 --> use next toe (most accurate)
%      Output:
%         - satpos = [X Y Z dt toe drel]
%           XYZ = ECEF coordinates (m)
%           dt = satellite clock bias (s)
%           toe = time of ephemerides (sec of GPS week)
%           drel = relativistic correction (m)
%           tgd = transmitter group delay (s)
%
%      Usage:   satpos = get_satpos(t,sv,eph,flag);

% Earth's universal gravitational parameter, m^3/s^2
GM = 3.986004418e14;

% earth rotation rate, rad/s
Omegae_dot = 7.2921151467e-5;

% speed of light
c = 0.299792458e9;

% find correct sv in eph matrix
col = find(eph(1,:)==sv);
eph = eph(:,col);

% make sure there is data for that PRN
if (isempty(col))
  disp(['WARNING: no data for PRN' num2str(sv)]);
  satpos = [];
  return;
end

% check that message transmission time is positive (should actually be close to toe)
I = find(eph(23,:) < 0);
if (~isempty(I))
  disp(['WARNING: found negative transmission time for PRN ' num2str(eph(1,I)) ' => removed']);
  J = find(eph(23,:) >= 0);
  eph = eph(:,J);
end

% find correct ephemerides for that sv
if (flag == 0)
  % the closest toe to the time of data (t):
  [y,col] = min(abs(t-eph(18,:)));
elseif (flag == 1)
  % use last toe (as GPS receiver would do):
  [y,col] = max(find((t-eph(18,:))>=0));
elseif (flag == 2)
  % use the first toe for the entire day:
  col = 1;
elseif (flag == 3)
  % use next toe (most accurate!):
  [col,y] = min(find((t-eph(18,:))<=0));
  if (isempty(col)) col = size(eph,2); end;
end

% make sure there is data for that time
if (isempty(col))
  %disp(['WARNING: no data for PRN' num2str(sv) ' for time ' num2str(t)]);
  satpos = [];
  return;
end

% assign the correct column
eph = eph(:,col);

% read ephemerides
svprn   =   eph(1);
af2     =   eph(2);
M0      =   eph(3);
roota   =   eph(4);
deltan  =   eph(5);
ecc     =   eph(6);
omega0  =   eph(7);
cuc     =   eph(8);
cus     =   eph(9);
crc     =  eph(10);
crs     =  eph(11);
i0      =  eph(12);
idot    =  eph(13);
cic     =  eph(14);
cis     =  eph(15);
Omega0  =  eph(16);
Omegadot=  eph(17);
toe     =  eph(18);
af0     =  eph(19);
af1     =  eph(20);
toc     =  eph(21);
tgd     =  eph(22);

% Now compute satellite position
A = roota*roota;

% time elapsed since toe
tk = t - toe;

% mean anomaly at t
M = M0 + (sqrt(GM/A^3)+deltan)*tk;

% iterative solution for E
E_old = M;
dE = 1;
%while (dE > 1e-12)
while (dE > 1e-11)
  E = M + ecc*sin(E_old);
  dE = abs(E-E_old);
  E_old = E;
end

% true anomaly
v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);

% same as: (Hoffmann-Wellenhof)
%v = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));

% other solution for the true anomaly (Tsui, 2000)
%v1 = acos((cos(E)-ecc)/(1-ecc*cos(E)));
%v2 = asin(((1-ecc^2)^.5)*sin(E)/(1-ecc*cos(E)));
%v = v1 * sign(v2);

% longitude of ascending node, accounting for Earth rotation Omegae_dot
Omega = Omega0 + (Omegadot-Omegae_dot)*tk - Omegae_dot*toe;

% intermediate quantities needed to compute corrections
phi = omega0 + v;
phi = rem(phi,2*pi);

% argument of perigee
omega = omega0 + cuc*cos(2*phi) + cus*sin(2*phi);

% radial distance
r = A*(1-ecc*cos(E)) + crc*cos(2*phi) + crs*sin(2*phi);

% inclination
i = i0 + idot*tk + cic*cos(2*phi) + cis*sin(2*phi);

% rotation matrix from orbital to Earth fixed
R = [cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i) -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i)  sin(Omega)*sin(i);
     sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i) -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i) -cos(Omega)*sin(i);
     sin(omega)*sin(i)                                   cos(omega)*sin(i)                                   cos(i)];

% satellite position vector w.r.t. Earth's center in inertial frame
satpos_in = [r*cos(v)
             r*sin(v)
             0];

% rotate into ECEF frame
satpos = R * satpos_in;

% output to check with students at t=0h and t=22h of the current day
format long
if (sv==31 && (t == 6*24*3600)) % PRN31 at 00:00:00
%if (sv==31 && (t == 6*24*3600 + 12*3600)) % PRN31 at 12:00:00
%if (sv==31 && (t == 6*24*3600 || t == 6*24*3600 + 24*3600 - 15*60))
%if (sv==2 && (t == 6*24*3600 || t == 6*24*3600 + 24*3600 - 15*60))
  display(['==========================================']);
  display(['t = ' num2str(t) ', tk = ' num2str(tk) ',  toe = ' num2str(toe)]);
  display(['sat pos in orbital plane = ' num2str(satpos_in')]);
  display(['sat pos ECEF = ' num2str(satpos')]);
  display(['M0 = ' num2str(M0) ', M = ' num2str(M) ', E = ' num2str(E) ', v = ' num2str(v) ', OMEGA = ' num2str(Omega)]);
ecc
M0
M
E
v
R
Omega
omega
i
  display(['omega = ' num2str(omega) ', r = ' num2str(r) ', i = ' num2str(i)]);
  display(['R = [' num2str(R(1,:))]);
  display(['     ' num2str(R(2,:))]);
  display(['     ' num2str(R(3,:)) ']']);
  display(['==========================================']);
end

% main relativistic correction
drel = sqrt(GM) * roota * ecc * sin(E) * 2 / c;

% compute and append satellite clock bias in seconds
dt = af0 + af1*tk + af2*tk^2;

% output
satpos = [satpos; dt; toe; drel; tgd];

