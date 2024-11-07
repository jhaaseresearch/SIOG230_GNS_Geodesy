%Homework 05

%Open Rinex file
eph_all =read_rinexn("brdc2920.19n")

[obs,t,prn,apr,hant] =read_rinexo("opmt2920.19o")

%Define constants (c), a priori GPS receiver position (X0,Y0,Z0) and satellite elevation cutoff
%angle (10◦ for instance);

c=0.299792458*10^9
el_cutoff= 10

%Each 15 minutes in seconds after the first date in the file
%Minimum time - Seconds of week SoGPSW_min
[Y_min,M_min,DoM_min,DoY_min,GPSW_min,DoGPSW_min,SoGPSW_min,JD_min,DecY_min] = gpsdate(t(1,1)+2000,t(1,2),t(1,3),t(1,4),t(1,5),t(1,6))

%Maximum time - Seconds of week SoGPSW_max
[Y_max,Mmax,DoMmax,DoYmax,GPSWmax,DoGPSWmax,SoGPSW_max,JDmax,DecYmax] = gpsdate(t(end,1)+2000,t(end,2),t(end,3),t(end,4),t(end,5),t(end,6))

T_vector=[SoGPSW_min:900:SoGPSW_max]'

start_time=datetime(t(1,1)+2000,t(1,2),t(1,3),t(1,4),t(1,5),t(1,6))
end_time=datetime(t(end,1)+2000,t(end,2),t(end,3),t(end,4),t(end,5),t(end,6))
timeVec = start_time:minutes(15):end_time;

%For each time in the T_vector, calculate the XYZ-get_satpos   
N=size(T_vector,1) 
txyz = zeros(N,4)

for time_i=1:N
    
    year_d=year(timeVec(time_i))-2000
    ind_time = find(t(:, 1) == year_d  & ...
                t(:, 2) == month(timeVec(time_i)) & ...
                t(:, 3) == day(timeVec(time_i)) & ...
                t(:, 4) == hour(timeVec(time_i)) & ...
                t(:, 5) == minute(timeVec(time_i)) & ...
                t(:, 6) == second(timeVec(time_i)));    
    
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
    
    %Choose the satellite and calculate the xyz
    svprn=prn(sv_ind)
    M=size(svprn, 1) 
    %if M > 4 continue: 
    t_data= T_vector(time_i) * ones(M,1) - jP/c
    c_s=ones(M,1)*c*1e-9
    X0_adj=apr(1,1) %Initial guess
    Y0_adj=apr(1,2) %Initial guess
    Z0_adj=apr(1,3); %Initial guess
    tolerance = 18; 
    maxIter = 1; 
    iter = 0;
    while iter < maxIter
        Diff_X= zeros(M, 1)
        Diff_Y= zeros(M, 1)
        Diff_Z= zeros(M, 1)
        rho_0 = zeros(M,1)
        A = zeros(M, 4);
        L = zeros(M, 1);
        jaX = zeros(M, 1);
        jaY = zeros(M, 1);
        jaZ = zeros(M, 1);
        iter = iter + 1;
        for i=1:M
            [x,y,z, delta_t]=get_satpos(t_data(i),svprn(i),eph_all)
            Diff_X(i)= x - X0_adj
            Diff_Y(i)= y - Y0_adj
            Diff_Z(i)= z - Z0_adj
            rho_0(i)=sqrt(Diff_X(i)^2 + Diff_Y(i)^2 + Diff_Z(i)^2)
            jaX(i) = -(Diff_X(i)) / rho_0(i);
            jaY(i) = -(Diff_Y(i)) / rho_0(i);
            jaZ(i) = -(Diff_Z(i)) / rho_0(i);
            L(i) = jP(i) - rho_0(i) + c * delta_t
            A(i,:) = [jaX(i), jaY(i), jaZ(i), -1*c_s(i)];
        end
        X = pinv(A) * L;  
        %Compute the covariance of the unknowns in ECEF frame;
        C = pinv(A' * A);
        X0_old = X0_adj;
        Y0_old = Y0_adj;
        Z0_old = Z0_adj;
        X0_adj=X0_old + X(1,1);  % Adjustments
        Y0_adj=Y0_old+ X(2,1);  % Adjustments
        Z0_adj=Z0_old + X(3,1);  % Adjustments 
        if  all(abs([X0_adj - X0_old, Y0_adj - Y0_old, Z0_adj - Z0_old]) < tolerance)
            break;
        end
    end
    txyz(time_i,:) = [T_vector(time_i), X0_adj, Y0_adj, Z0_adj];
end
    
    
    
    

