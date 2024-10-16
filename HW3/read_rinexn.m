function [eph, alpha, beta, dow] = read_rinexn(rinexn)

% READ_RINEXN	Reads a RINEX navigation file.
%
% Input: rinexn  = RINEX navigation file name
%
% Output: eph = matrix with 23 rows and as many columns
%               as there are ephemerides.
%         alpha = ionospheric coefficients (4 element vector)
%         beta  = ionospheric coefficients (4 element vector)
%         dow   = day of week

% open navigation file
disp(['-------------']);
fid = fopen(rinexn, 'r');
disp(['Reading ' rinexn]);

% initialize output variables
alpha = [];
beta = [];
eph = [];

% read header
line = fgetl(fid);
while ~(strcmp(line(61:63), 'END'))
    line = fgetl(fid);
    
    if (contains(line, 'ION ALPHA'))
        line = strrep(line, 'D', 'E');
        alpha = [str2double(line(4:14)) str2double(line(16:26)) str2double(line(28:38)) str2double(line(40:50))];
    end
    
    if (contains(line, 'ION BETA'))
        line = strrep(line, 'D', 'E');
        beta = [str2double(line(4:14)) str2double(line(16:26)) str2double(line(28:38)) str2double(line(40:50))];
    end
end

% read ephemeris body
j = 1;
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    
    % Convert "D" to "E" in scientific notation for correct parsing
    line = strrep(line, 'D', 'E');
    tmp = sscanf(line, '%d %d %d %d %d %d %f');
    
    % find day of week
    if (j == 1)
        [Y, M, DoM, DoY, GPSW, dow, SoGPSW, JD, DecY] = gpsdate(tmp(2) + 2000, tmp(3), tmp(4), tmp(5), tmp(6), tmp(7));
    end
    
    % parse ephemeris data for each satellite
    svprn = tmp(1);
    toc = tmp(5) + tmp(6) / 60 + tmp(7) / 3600;
    
    % read the clock and orbit parameters, converting "D" to "E"
    clock_bias = str2double(strrep(line(23:41), 'D', 'E'));
    clock_drift = str2double(strrep(line(42:60), 'D', 'E'));
    clock_drift_rate = str2double(strrep(line(61:79), 'D', 'E'));
    
    % read other parameters from the following lines
    line = strrep(fgetl(fid), 'D', 'E');
    iode = str2double(line(4:22));
    crs = str2double(line(23:41));
    deltan = str2double(line(42:60));
    M0 = str2double(line(61:79));
    
    line = strrep(fgetl(fid), 'D', 'E');
    cuc = str2double(line(4:22));
    ecc = str2double(line(23:41));
    cus = str2double(line(42:60));
    roota = str2double(line(61:79));
    
    line = strrep(fgetl(fid), 'D', 'E');
    toe = str2double(line(4:22));
    cic = str2double(line(23:41));
    Omega0 = str2double(line(42:60));
    cis = str2double(line(61:79));
    
    line = strrep(fgetl(fid), 'D', 'E');
    i0 = str2double(line(4:22));
    crc = str2double(line(23:41));
    omega = str2double(line(42:60));
    Omegadot = str2double(line(61:79));
    
    line = strrep(fgetl(fid), 'D', 'E');
    idot = str2double(line(4:22));
    codesL2 = str2double(line(23:41));
    gpsweek = str2double(line(42:60));
    L2Pflag = str2double(line(61:79));
    
    line = strrep(fgetl(fid), 'D', 'E');
    sv_accuracy = str2double(line(4:22));
    sv_health = str2double(line(23:41));
    tgd = str2double(line(42:60));
    iodc = str2double(line(61:79));
    
    line = strrep(fgetl(fid), 'D', 'E');
    trans = str2double(line(4:22));
    
    % store ephemeris data
    eph(1, j) = svprn;
    eph(2, j) = clock_drift_rate;
    eph(3, j) = M0;
    eph(4, j) = roota;
    eph(5, j) = deltan;
    eph(6, j) = ecc;
    eph(7, j) = omega;
    eph(8, j) = cuc;
    eph(9, j) = cus;
    eph(10, j) = crc;
    eph(11, j) = crs;
    eph(12, j) = i0;
    eph(13, j) = idot;
    eph(14, j) = cic;
    eph(15, j) = cis;
    eph(16, j) = Omega0;
    eph(17, j) = Omegadot;
    eph(18, j) = toe;
    eph(19, j) = clock_bias;
    eph(20, j) = clock_drift;
    eph(21, j) = toc;
    eph(22, j) = tgd;
    eph(23, j) = trans;
    
    j = j + 1;
end

fclose(fid);
end
