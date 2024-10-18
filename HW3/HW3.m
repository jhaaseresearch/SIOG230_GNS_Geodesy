%Homework 03


%Compute ECEF coordinates, sv=31, and compare with coordinates in the
% sp3 file

%For the day October 19

%Minimum time - Seconds of week SoGPSW_min
[Y_min,M_min,DoM_min,DoY_min,GPSW_min,DoGPSW_min,SoGPSW_min,JD_min,DecY_min] = gpsdate(2019,10,19,0,0,0)

%Maximum time - Seconds of week SoGPSW_max
[Y_max,Mmax,DoMmax,DoYmax,GPSWmax,DoGPSWmax,SoGPSW_max,JDmax,DecYmax] = gpsdate(2019,10,19,23,45,0)

%Open Rinex file
eph_all =read_rinexn("brdc2920.19n")


%Open SP3 file
[sp3,sv,excl_sat] = read_sp3("igs20756.sp3")
%For sv = 31 
T_s=sp3.prn31(:,1)
X_s=sp3.prn31(:,2)
Y_s=sp3.prn31(:,3)
Z_s=sp3.prn31(:,4)

% 15 minutes are 900 seconds and they are in the T_s vector
T_vector= SoGPSW_min + T_s

%For each time in the T_vector, calculate the XYZ-get_satpos and the
%difference between XYZ-get_satpos and sp3

N=size(T_vector, 1) 

Diff_X= zeros(N, 1)
Diff_Y= zeros(N, 1)
Diff_Z= zeros(N, 1)
Vector_Diff = zeros(N, 1)
Range_Diff = zeros(N, 1)

for i=1:N
    [x,y,z]=get_satpos(T_vector(i,1),31,eph_all)
    Diff_X(i)= x - X_s(i,1) 
    Diff_Y(i)= y - Y_s(i,1) 
    Diff_Z(i)= z - Z_s(i,1)
    Vector_Diff(i)=sqrt(Diff_X(i)^2 + Diff_Y(i)^2 +Diff_Z(i)^2)
    Range_Diff(i)=sqrt(x^2+y^2+z^2) - sqrt(X_s(i,1)^2+Y_s(i,1)^2+Z_s(i,1)^2)
end



