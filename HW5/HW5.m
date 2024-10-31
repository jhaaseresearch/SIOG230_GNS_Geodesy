%Homework 05

%Open Rinex file
eph_all =read_rinexn("brdc2920.19n")

[obs,t,prn,apr,hant] =read_rinexo("opmt2920.19o")

[sp3,sv,excl_sat] = read_sp3("igs20756.sp3")


%Define constants (c), a priori GPS receiver position (X0,Y0,Z0) and satellite elevation cutoff
%angle (10◦ for instance);

c=0.299792458*10^9
el_cutoff= 10


%After I will create a loop here with time
%For each time in the T_vector, calculate the XYZ-get_satpos   
%N=size(T_vector, 1) 
%txyz = zeros(N,5


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
[Y_t,M_t,DoM_t,DoY_t,GPSW_t,DoGPSW_t,SoGPSW_t,JD_t,DecY_t] = gpsdate(2019,10,19,0,0,0)

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

%Use your function get_satpos to calculate (from broadcast ephemerides) the satellite positions 
% in the ECEF frame and satellite clock errors (see equation 2.16) at time of transmission ttx. 
% Make sure they are expressed in meters and seconds.
ttx = seconds_time - jP/c;
%Choose the satellite and calculate the xyz
svprn=prn(sv_ind)
M=size(svprn, 1) 
%if M > 4 continue: 
t_data= SoGPSW_t
c_s=ones(M,1)*c*1e-9
X0_adj=apr(1,1) %Initial guess
Y0_adj=apr(1,2) %Initial guess
Z0_adj=apr(1,3); %Initial guess
tolerance = 1; 
maxIter = 100; 
iter = 0;
while iter < maxIter
    Diff_X= zeros(M, 1)
    Diff_Y= zeros(M, 1)
    Diff_Z= zeros(M, 1)
    rho_0 = zeros(M,1)
    A = zeros(M, 4);
    jaX = zeros(M, 1);
    jaY = zeros(M, 1);
    jaZ = zeros(M, 1);
    X0_old = X0_adj;
    Y0_old = Y0_adj;
    Z0_old = Z0_adj;
    iter = iter + 1;
    for i=1:M
        [x,y,z, delta_t]=get_satpos(t_data,svprn(i),eph_all)
        Diff_X(i)= x - X0
        Diff_Y(i)= y - Y0
        Diff_Z(i)= z - Z0
        rho_0(i)=sqrt(Diff_X(i)^2 + Diff_Y(i)^2 +Diff_Z(i)^2)
        jaX(i) = -(Diff_X(i)) / rho_0(i);
        jaY(i) = -(Diff_Y(i)) / rho_0(i);
        jaZ(i) = -(Diff_Z(i)) / rho_0(i);
        L = jP(i) - rho_0(i) + c_s * delta_t
        A = [jaX, jaY, jaZ, -1*c_s];
    end
    X = pinv(A) * L;  
    %Compute the covariance of the unknowns in ECEF frame;
    C = pinv(A' * A);
    X0_adj=X0 + X(1,1);  % Adjustments
    Y0_adj=Y0 + X(2,1);  % Adjustments
    Z0_adj=Z0 + X(3,1);  % Adjustments

    if abs(X0_adj-X0_old) && abs(Y0_adj-Y0_old) && abs(Z0_adj-Z0_old) < tolerance
        break;
    end
end






