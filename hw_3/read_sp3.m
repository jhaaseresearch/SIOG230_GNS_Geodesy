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

% open, read, and close sp3 file
fid = fopen(sp3file,'r');
while 1
  line = fgetl(fid);
  if ~isstr(line), break, end
  c = line(1:2);
  switch c
    case '+ '
      deli = line(10);
      if (~isempty(line(5:6))) n_sv = str2num(line(5:6)); end;
      [svtmp] = strread(line(11:length(line)),'%f','delimiter',deli);
      sv = [sv;svtmp(find(svtmp))];
    case '* '
      [t] = sscanf(line, '%c %d %d %d %d %d %f');
      Ts = t(5)*3600+t(6)*60+t(7);
    case 'PG'
      [p] = sscanf(line(3:end), '%2d %lf %lf %lf %f');
      % convert sp3 to meters
      Xs = p(2)*1000; Ys = p(3)*1000; Zs = p(4)*1000;
      % convert sat clock bias to seconds
      dT = p(5)*1e-6;
      % satellite number
      prn = line(3:4);
      % if satellite ok
      if (p(5) ~= 999999.999999)
        % build one matrix per satellite
        eval(['if ~exist(''sp3_' prn ''') sp3_' prn ' = []; end;']);
        eval(['sp3_' prn '= [sp3_' prn '; Ts Xs Ys Zs dT];']);
      else
        excl_sat = [excl_sat ; prn];
      end
  end
end
fclose(fid);

% list of bad satellites
excl_sat = unique(str2num(excl_sat));

% assign data to structure
sat = [];
for i=1:length(sv)
  if isempty(find(excl_sat==sv(i)))
    prn = num2str(sprintf('%.2d',sv(i)));
    eval(['sat_tmp = sp3_' prn ';']);
    field = ['prn' num2str(sv(i))];
    sp3 = setfield(sp3,field,sat_tmp);
    sat = [sat ; sv(i)];
  end
end
sv = sat;
