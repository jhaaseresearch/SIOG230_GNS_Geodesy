% Neil Waldhausen
% SIOG 239 Homeowork 2
% Converting RINEX file to ECEF coordinates

function [eph,alpha,beta,dow] = read_rinexn(rinexn)

% READ_RINEXN	Reads a rinex navigation file.
%
% Input: rinexn  = rinex navigation file name
%
% Output: eph = matrix with 23 rows and as many columns
%               as there are ephemerides.
%         alpha = ionospheric coefficients (4 element vector)
%         beta  = ionospheric coefficients (4 element vector)
%         dow   = day of week
%
% Usage:  eph = read_rinexn(rinexn)

% open navigation file
disp(['-------------']);
fid = fopen(rinexn,'r');
disp(['Reading ' rinexn]);

% read header
alpha = [];
beta = [];
line = fgetl(fid);
while ~(strcmp(line(61:63),'END'))
  line = fgetl(fid);
  if (strfind(line,'ION ALPHA'))
    alpha = [str2num(line(4:14)) str2num(line(16:26)) str2num(line(28:38)) str2num(line(40:50))];
  end
  if (strfind(line,'ION BETA'))
    beta = [str2num(line(4:14)) str2num(line(16:26)) str2num(line(28:38)) str2num(line(40:50))];
  end
end

% read body
j = 1;
while 1
   line = fgetl(fid);
   if ~isstr(line), break, end
   tmp = sscanf(line, '%d %d %d %d %d %d %f');

   % find day of week
   if (j==1)
     [Y,M,DoM,DoY,GPSW,dow,SoGPSW,JD,DecY] = gpsdate(tmp(2)+2000,tmp(3),tmp(4),tmp(5),tmp(6),tmp(7));
   end

   svprn = tmp(1);
   toc = tmp(5) + tmp(6)/60 + tmp(7)/3600;
   clock_bias = str2num(line(23:41));
   clock_drift = str2num(line(42:60));
   clock_drift_rate = str2num(line(61:79));

   line = fgetl(fid);
   iode = str2num(line(4:22));
   crs = str2num(line(23:41));
   deltan = str2num(line(42:60));
   M0 = str2num(line(61:79));

   line = fgetl(fid);
   cuc = str2num(line(4:22));
   ecc = str2num(line(23:41));
   cus = str2num(line(42:60));
   roota = str2num(line(61:79));

   line = fgetl(fid);
   toe = str2num(line(4:22));
   cic = str2num(line(23:41));
   Omega0 = str2num(line(42:60));
   cis = str2num(line(61:79));

   line = fgetl(fid);
   i0 = str2num(line(4:22));
   crc = str2num(line(23:41));
   omega = str2num(line(42:60));
   Omegadot = str2num(line(61:79));

   line = fgetl(fid);
   idot = str2num(line(4:22));
   codesL2 = str2num(line(23:41));
   gpsweek = str2num(line(42:60));
   L2Pflag = str2num(line(61:79));

   line = fgetl(fid);
   sv_accuracy = str2num(line(4:22));
   sv_health = str2num(line(23:41));
   tgd = str2num(line(42:60));
   iodc = str2num(line(61:79));

   line = fgetl(fid);
   trans = str2num(line(4:22));

   eph(1,j) = svprn;
   eph(2,j) = clock_drift_rate;
   eph(3,j) = M0;
   eph(4,j) = roota;
   eph(5,j) = deltan;
   eph(6,j) = ecc;
   eph(7,j) = omega;
   eph(8,j) = cuc;
   eph(9,j) = cus;
   eph(10,j) = crc;
   eph(11,j) = crs;
   eph(12,j) = i0;
   eph(13,j) = idot;
   eph(14,j) = cic;
   eph(15,j) = cis;
   eph(16,j) = Omega0;
   eph(17,j) = Omegadot;
   eph(18,j) = toe;
   eph(19,j) = clock_bias;
   eph(20,j) = clock_drift;
   eph(21,j) = toc;
   eph(22,j) = tgd;
   eph(23,j) = trans;

   j=j+1;
end

fclose(fid);

function [sp3,sv,excl_sat] = read_sp3(sp3file);

% READ_SP3  Reads sp3 file and returns position vectors.
%
% Input: sp3 = orbit file name
%
% Output: sp3 = structure with for each field (= prn):
%               T = time (in seconds of the day)
%               X,Y,Z (in meters) for each field (= prn),
%               dT = satellite clock (in seconds)
%         sv = list of prn numbers present in sp3 file (vector)
%         excl_sat = list of bad satellites
%
% Usage: [sp3,sv,excl_sat] = read_sp3(sp3file);

% initialize sp3 structure
sp3 = [];
sv  = [];
excl_sat = [];