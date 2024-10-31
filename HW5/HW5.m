%Homework 05

%Open Rinex file
eph_all =read_rinexn("brdc2920.19n")

[obs,t,prn,apr,hant] =read_rinexo("opmt2920.19o")

[sp3,sv,excl_sat] = read_sp3("igs20756.sp3")


%Define constants (c), a priori GPS receiver position (X0,Y0,Z0) and satellite elevation cutoff
%angle (10◦ for instance);

c=0.299792458*10^9
el_cutoff= 10

%Read the pseudorange data jP and its corresponding epoch of observation, which you convert into seconds of the current day. 
%Again, I would start by hard–coding the pseudorange data from the rinex file for one epoch. You can later replace this by 
%the output of function read_rinexo

ind_time = find(t(:, 1) == 19 & ...
                t(:, 2) == 10 & ...
                t(:, 3) == 19 & ...
                t(:, 4) == 0 & ...
                t(:, 5) == 0 & ...
                t(:, 6) == 0);

date_time=datetime(t(ind_time,1), t(ind_time,2), t(ind_time,3), t(ind_time,4), t(ind_time,5), t(ind_time,6))
seconds_time = seconds(date_time - dateshift(date_time, 'start', 'day'))


%Compute the time of transmission of the data, which is not the same as the time of observation 
% in the receiver because the signal has traveled the distance jP between the satellite and the 
% receiver. In other words, pseudorange data tagged time = t in the rinex observation file were 
% sent a bit earlier, at ttx = t − j P /c. That is the time at which you will need to compute the 
% satellite position.

%For C1 
jP_all = obs.C1(ind_time, :);
nonNaN = ~isnan(jP_all);
sv_ind = find(nonNaN);
jP=jP_all(sv_ind)'
ttx = seconds_time - jP / c;

%Use your function get_satpos to calculate (from broadcast ephemerides) the satellite positions 
% in the ECEF frame and satellite clock errors (see equation 2.16) at time of transmission ttx. 
% Make sure they are expressed in meters and seconds.

%Choose the satellite and calculate the xyz
svprn=prn(sv_ind)
M=size(svprn, 1) 
eph_sel=struct(); 

for sv_i=1:M
    [Y_t,M_t,DoM_t,DoY_t,GPSW_t,DoGPSW_t,SoGPSW_t,JD_t,DecY_t] = gpsdate(2019,10,19,0,15,0)
    T_vector= SoGPSW_t

    %For each time in the T_vector, calculate the XYZ-get_satpos   
    N=size(T_vector, 1) 
    txyz = zeros(N,4)
    
    for i=1:N
        [x,y,z]=get_satpos(T_vector(i,1),svprn(sv_i),eph_all)
        txyz(i, :)=[T_vector(i,1),x,y,z]
    end
    
    svprn_number=svprn(sv_i)
    eph_sel.(['prn', num2str(svprn_number)]) = txyz;    
end


%XYZ_0 = zeros(M,4)
%Selecting X0, Y0, and Z0 
%for prn_i=1:M
%    i_time=find(sp3.(['prn', num2str(svprn(prn_i))])(:,1) == seconds_time);
%    X0=sp3.(['prn', num2str(svprn(prn_i))])(i_time,2)
%    Y0=sp3.(['prn', num2str(svprn(prn_i))])(i_time,3)
%    Z0=sp3.(['prn', num2str(svprn(prn_i))])(i_time,4)
%    XYZ_0(prn_i,:)=[svprn(prn_i), X0, Y0, Z0]
%end

%XO, YO, Z0 for satellite 2
i_time=find(sp3.(['prn', num2str(svprn(prn_i))])(:,1) == seconds_time);
X0=sp3.prn2(i_time,2)
Y0=sp3.prn2(i_time,3)
Z0=sp3.prn2(i_time,4)

%Selecting X, Y and Z 
XYZ = zeros(M,4)
for prn_i=1:M
%    i_time=find(sp3.(['prn', num2str(svprn(prn_i))])(:,1) == seconds_time);
    X=eph_sel.(['prn', num2str(svprn(prn_i))])(1,2)
    Y=eph_sel.(['prn', num2str(svprn(prn_i))])(1,3)
    Z=eph_sel.(['prn', num2str(svprn(prn_i))])(1,4)
    XYZ(prn_i,:)=[svprn(prn_i), X, Y, Z]
end


%Compute the modeled observables jρ0;

Diff_X= zeros(M, 1)
Diff_Y= zeros(M, 1)
Diff_Z= zeros(M, 1)
rho_0 = zeros(M,1)
%for i_sv=1:M 
%    Diff_X(i_sv)= XYZ(i_sv,2) - XYZ_0(i_sv,2) 
%    Diff_Y(i_sv)= XYZ(i_sv,3) - XYZ_0(i_sv,3) 
%    Diff_Z(i_sv)= XYZ(i_sv,4) - XYZ_0(i_sv,4)
%    rho_0(i_sv)=sqrt(Diff_X(i_sv)^2 + Diff_Y(i_sv)^2 +Diff_Z(i_sv)^2)
%end

for i_sv=1:M 
    Diff_X(i_sv)= XYZ(i_sv,2) - X0
    Diff_Y(i_sv)= XYZ(i_sv,3) - Y0
    Diff_Z(i_sv)= XYZ(i_sv,4) - Z0
    rho_0(i_sv)=sqrt(Diff_X(i_sv)^2 + Diff_Y(i_sv)^2 +Diff_Z(i_sv)^2)
end



%Compute the observation vector by subtracting known quantities from the pseudorange observation:
% L = j P − iρ0 + c j δ. Warning: satellite clock biases must be added to the measured pseudoranges;

%Assuming no clock biases because I did not added it to the getsapos yet
c_s=ones(M,1)*c*1e-9
%For C1 
delta_t=-0.00032874 


L = jP - rho_0 + c * delta_t

%Compute the partial derivatives and build the design matrix A. Trick: multiply c by 10−9 in 
% the design matrix in order to avoid numerical instabilities in the inversion. 
% The receiver clock bias will be output in nanoseconds.

A = zeros(length(jP), 4);
jaX = zeros(M, 1);
jaY = zeros(M, 1);
jaZ = zeros(M, 1);

for j = 1:length(jP)
    jaX(j) = -(Diff_X(j)) / rho_0(j);
    jaY(j) = -(Diff_Y(j)) / rho_0(j);
    jaZ(j) = -(Diff_Z(j)) / rho_0(j);
end

A = [jaX, jaY, jaZ, -c_s];
%Invert the design matrix and find the vector of unknowns, or solve the least squares problem 
% directly. In MATLAB this can be done using functions pinv or lscov;

X = pinv(A) * L;  

%Compute the covariance of the unknowns in ECEF frame;
C = pinv(A' * A);

X_final = [X0; Y0; Z0] + X(1:3);  % Adjustments




