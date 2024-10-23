function [obs,t,gps,apr,hant] = read_rinexo(rinexo);

% READ_RINEXO Read RINEX observation file. Slower than readrinex but
%             does not need compiling. Reads rinex files with Glonass
%             and Galileo data as well.
%
% Input: rinexo = rinex file name
%
% Output: obs = structure with observation types, e.g. obs.C1, obs.L1, etc.
%         t = time matrix (year month day hour minute second)
%         gps = vector of GPS satellite numbers
%         apr = approximate site position from rinex file header
%
% Usage: [obs,t,prn,apr,hant] = read_rinexo(rinexo);

% defaults
lineflag = 1;

% initialize outputs
obs  = [];
t    = [];
gps  = [];
glo  = [];
gal  = [];
hant = [0 0 0];

% open rinex file
disp(['-------------']);
fid = fopen(rinexo,'r');
disp(['Reading ' rinexo]);

% read header information
while(~feof(fid))
  line = fgetl(fid);
  if strfind(line,'APPROX POSITION XYZ')  apr = str2num(line(1:60)); end;
  if strfind(line,'INTERVAL')             samplint = str2num(line(5:10)); end;
  if strfind(line,'ANTENNA: DELTA H/E/N') hant = str2num(line(1:60)); end;
  if strfind(line,'# / TYPES OF OBSERV')
    nobs_type = str2num(line(6));
    if (nobs_type>5) lineflag = 2; end
    [obs_type] = strread(line(8:60),'%s','delimiter',' ');
  end
  if strfind(line,'END OF HEADER') break; end;
end
if ~exist('apr') warning('cannot find approximate position'); end;
disp(['  approximate position (m): ' num2str(apr(1)) ' ' num2str(apr(2)) ' ' num2str(apr(3))]);
if ~exist('nobs_type') error('cannot find observation types'); end;
disp(['  ' num2str(nobs_type) ' observable types: ' cell2mat(obs_type')]);
if ~exist('samplint') error('cannot find sampling interval'); end;
%samplint = 30;
disp(['  sampling interval: ' num2str(samplint) ' seconds']);

% determine number of epochs and prns
n_ep = 0;
while(~feof(fid))
  % find start of data block
  line = fgetl(fid);
  if (length(line)>=33 & ~isempty(regexp(line(33),'[GRE]')))
    n_ep = n_ep + 1;
    % get total number of satellites
    nsat = str2num(line(31:32));
    % extract satellites numbers
    sat_list = [];
    sat_tmp = regexp(line(33:length(line)),'\w{0,3}', 'match');
    sat_list = [sat_list sat_tmp];
    gps_tmp = cell2mat(regexp(line(33:length(line)),'G..','match'));
    for i=2:3:length(gps_tmp) gps = [gps ; str2num(gps_tmp(i:i+1))]; end;
    glo_tmp = cell2mat(regexp(line(33:length(line)),'R..','match'));
    for i=2:3:length(glo_tmp) glo = [glo ; str2num(glo_tmp(i:i+1))]; end;
    gal_tmp = cell2mat(regexp(line(33:length(line)),'E..','match'));
    for i=2:3:length(gal_tmp) gal = [gal ; str2num(gal_tmp(i:i+1))]; end;
    % if more lines with satellite numbers
    while (length(sat_list) < nsat)
      line = fgetl(fid);
      sat_tmp = regexp(line(33:length(line)),'\w{0,3}', 'match');
      sat_list = [sat_list sat_tmp];
      gps_tmp = cell2mat(regexp(line(33:length(line)),'G..','match'));
      for i=2:3:length(gps_tmp) gps = [gps ; str2num(gps_tmp(i:i+1))]; end;
      glo_tmp = cell2mat(regexp(line(33:length(line)),'R..','match'));
      for i=2:3:length(glo_tmp) glo = [glo ; str2num(glo_tmp(i:i+1))]; end;
      gal_tmp = cell2mat(regexp(line(33:length(line)),'E..','match'));
      for i=2:3:length(gal_tmp) gal = [gal ; str2num(gal_tmp(i:i+1))]; end;
    end
  end
end
disp(['  ' num2str(n_ep) ' epochs']);
gps = unique(gps); n_gps = length(gps); disp(['  ' num2str(n_gps) ' GPS satellites']);
glo = unique(glo); n_glo = length(glo); disp(['  ' num2str(n_glo) ' GLONASS satellites']);
gal = unique(gal); n_gal = length(gal); disp(['  ' num2str(n_gal) ' GALILEO satellites']);
frewind(fid);

% initialize time and observable matrices
T  = zeros(n_ep,1);
for i=1:nobs_type
  eval(['O' num2str(i) ' = NaN(n_ep,n_gps);']);
end

% read body
j = 0;
disp(['  reading file body -- this part is slow, sorry...']);
while(~feof(fid))
  line = fgetl(fid);
  % found start of data block
  if (length(line)>=33 & ~isempty(regexp(line(33),'[GRE]')))
    j = j+1;
    % get time of observation
    [yr mo dd hr mn ss] = strread(line(2:26),'%d %d %d %d %d %f','delimiter',' ');
    t = [t ; yr mo dd hr mn ss];
    % get total number of satellites
    nsat = str2num(line(31:32));
    % extract all satelllites on line 1
    sat_list = regexp(line(33:length(line)),'\w{0,3}', 'match');
    % if there are more line with sat numbers, read them
    while (length(sat_list) < nsat)
      line = fgetl(fid);
      sat_tmp = regexp(line(33:length(line)),'\w{0,3}', 'match');
      sat_list = [sat_list sat_tmp];
    end
    % for each satellite, read observation data -- possibly two lines
    for i = 1:nsat
      % is this a GPS satellite?
      if ~isempty(regexp(sat_list{i},'G..'))
        % find column for this sat
        sv = str2num(sat_list{i}(2:3));
        col = find(gps == sv);
        % this part is slow
        line = fgetl(fid);
        if ((length(line))>=14 & ~isempty(str2num(line(1:14))))  O1(j,col) = str2num(line(1:14)); end;
        if ((length(line))>=30 & ~isempty(str2num(line(17:30)))) O2(j,col) = str2num(line(17:30)); end;
        if ((length(line))>=46 & ~isempty(str2num(line(33:46)))) O3(j,col) = str2num(line(33:46)); end;
        if ((length(line))>=62 & ~isempty(str2num(line(49:62)))) O4(j,col) = str2num(line(49:62)); end;
        if ((length(line))>=78 & ~isempty(str2num(line(65:78)))) O5(j,col) = str2num(line(65:78)); end;
        if (lineflag == 2)
          line = fgetl(fid);
          if ((length(line))>=14 & ~isempty(str2num(line(1:14))))  O6(j,col) = str2num(line(1:14)); end;
          if ((length(line))>=30 & ~isempty(str2num(line(17:30)))) O7(j,col) = str2num(line(17:30)); end;
          if ((length(line))>=46 & ~isempty(str2num(line(33:46)))) O8(j,col) = str2num(line(33:46)); end;
          if ((length(line))>=62 & ~isempty(str2num(line(49:62)))) O9(j,col) = str2num(line(49:62)); end;
          if ((length(line))>=62 & ~isempty(str2num(line(65:78)))) O10(j,col) = str2num(line(65:78)); end;
        end
      % if not a GPS satellite, skip line
      else
        line = fgetl(fid);
      end
    end
  end
end
fclose(fid);

% assign observations to proper observable name
for i=1:nobs_type
  eval (['obs_tmp = O' num2str(i) ';']);
  field = [obs_type{i}];
  obs = setfield(obs,field,obs_tmp);
end

