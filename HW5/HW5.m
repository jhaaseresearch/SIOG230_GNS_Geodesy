%{
 HW5

Author: Tommy Stone
Class: SIOG 239G
Date: Fall 2024

%(\___/)
%(=^.^=) A matrix is correct but L is a little off, not sure why
%(")_(")
%}

clear all; close all; clc;

% Adding utils created for class to this Matlab script
addpath('../matlab_utils')

% Adding products
addpath('../products')

% Add previous HW function
% addpath('../HW_Functions')

% Adding constants
c = 299792458; % speed of light, (m/s)
theta_cutoff = 10; % degrees is cutoff angel for Satellite

% File added from HW
rixen_file = 'opmt2920.19o';
brd_file = 'brdc2920.19n';

%% 2. Read the pseudorange data 
[obs,ts,gps,apr,hant] = read_rinexo(rixen_file);

% confirming with sp3 file
[sp3,sv,excl_sat] = read_sp3('igs20756.sp3');

% filtering eph data
[eph, alpha, beta, dow] = read_rinexn(brd_file);

% Assign variable for C1
C1 = obs.C1;

% convert times to seconds
timeDay = datetime(ts(:,1)+2000,ts(:,2),ts(:,3),ts(:,4),ts(:,5),ts(:,6));
refTime= datetime(2019,10,19,0,0,0);

% created reference time in seconds
timeSeconds = seconds(timeDay - refTime);

%%
%3. Compute time of transmission 
% This step will need to be done for P1, P2, and C1
ttx = timeSeconds - C1./c;

% Satellites list, this was taken from the book. Using the first four
prns = [2,6,12,14,24,25,29,32];
% when using all Satellites, prns = gps;
%prns = gps;

% Find indexes
for p = 1:length(prns)
    idx_prns(p) = find(gps == prns(p));
end

% filtering for Satellite positions
% over time we can get rid of this and use all Satellite positions
ttx = ttx(:,idx_prns);
C1 = C1(:,idx_prns);

%% 
% Obtain number of time records
T = length(ttx);

% initializing storage for Cxs
Cxs = nan([T,4,4]);
% Initializing R and Cxyz matrix
Cxyz = nan([T,3,3]);
R = nan([T,3,3]);
Cl = nan([T,3,3]);

% initializing dop values
vdops = nan([T,1]);
hdops = nan([T,1]);
pdops = nan([T,1]);
tdops = nan([T,1]);
gdops = nan([T,1]);


% Creating matrix to hold variables at each time step
% values
% delta,x0,y0,z0,error,num satellites,iterationsvalid
INFO = nan(T,7);

% iterate overtime to make calculations
start = 1;
for ii = 1:T
    ttx_temp = ttx(ii,:);
    
    % Find positions without nan values in ttx_temp
    % This should have at least 4 non Nan values
    idx_nn = find(~isnan(ttx_temp));
    
    % number of valid values
    num_real = length(idx_nn);

    % establishing initial guess coordinates
    x0 = apr(1);
    y0 = apr(2);
    z0 = apr(3);

    % tolerance should be within 1m
    tol = 0.5;
    error = 5;
    itr = 0;
    
    % If value is greater or equal to 4 then we can solve, if not we will
    % skip the record if not
    if num_real >= 4
        

        % determine the Satellites to use
        svs = prns(idx_nn);

        % initialize matrix A and vector l
        A = nan([num_real,4]);

        l = nan([num_real,1]);
        
        % get the time component for each
        t_real = ttx(ii,idx_nn);
        C1_real = C1(ii,idx_nn);

        % using least squares to solve. 
        while (tol < error && itr <= 20)

            itr = itr + 1;

            % iterate through each Satellite 
            % to build matrix A and l
            for idx_sv = 1:length(svs)
                sv = svs(idx_sv);
                t = t_real(idx_sv);

                t = t + 3600*24*6;
                
                % testing to see if first iteration works
                %if ii == 1
                % This code is for testing initially with sp3
                % we do not need this now that the code is working
                % I am leaving just in case though
                if 0 == 1
                        if sv == 2
                            satpos = sp3.prn2(ii,2:end);
                        elseif sv == 6
                            satpos = sp3.prn6(ii,2:end);
                        elseif sv == 12
                            satpos = sp3.prn12(ii,2:end);
                        elseif sv == 14
                            satpos = sp3.prn14(ii,2:end);
                        elseif sv == 24
                            satpos = sp3.prn24(ii,2:end);
                        elseif sv == 25
                            satpos = sp3.prn25(ii,2:end);
                        elseif sv == 29
                            satpos = sp3.prn29(ii,2:end);
                        elseif sv == 32
                            satpos = sp3.prn32(ii,2:end);
                        end
                else
                % [X Y Z dt toe drel]
                    satpos = get_satpos(t,sv,eph,3);
                end
                xj = satpos(1);
                yj = satpos(2);
                zj = satpos(3);
                deltaj = satpos(4);       % Satellite clock bias

                % calling object with sp3 file
                
    
                % Determining constants required to build matrix for each
                % Satellite
                rho0j = sqrt((xj - x0).^2 +(yj - y0).^2 +(zj - z0).^2);
    
                C1j = C1_real(idx_sv);
    
                % Calculate l for Satellite j
                % Multiplier is to avoid instabilities
                lj = C1j - rho0j + c*deltaj;
    
                % Calculate the constants for the matrix
                axj = -(xj - x0)/rho0j;
                ayj = -(yj - y0)/rho0j;
                azj = -(zj - z0)/rho0j;
    
                % Saving values for A and l
                % multiplying c by 10^-9 to avoid instability
                A(idx_sv,:) = [axj ayj azj -c*10^(-9)];
                l(idx_sv) = lj;
            end

            % Computing least square solution
            AT = A';
            Cx = pinv(AT*A);      % covaraince matrix of unknowns            
            dX = Cx*AT*l;
    
            % Determine tolerance
            x0_new = x0 + dX(1);
            y0_new = y0 + dX(2);
            z0_new = z0 + dX(3);
            delta = dX(4);
            error = sqrt((x0 - x0_new)^2 + (y0 - y0_new)^2 + (z0 - z0_new)^2);
    
            % Assigning new values for gradient descent
            x0 = x0_new;
            y0 = y0_new;
            z0 = z0_new;

        end


    % saving nan values in area where we cannot use least squares to solve
    else
        Cx = nan([4,4]);
        x0 = nan;
        y0 = nan;
        z0 = nan;
        delta = nan;
    end
    % matrix A and l have been built, solving for x
    if ii == 1
        % Storing values for comparison with HW at epoch 0:00
        A_store = A;
        l_store = l;
    end

    % Storing the gradient for later use 
    Cxs(ii,:,:) = Cx;

    Cxyz(ii,:,:) = Cx(1:3,1:3);

    % storing data for information
    INFO(ii,:) = [delta x0 y0 z0 error num_real itr];
end

% Taking data out of matrix
deltas = INFO(:,1);
X = INFO(:,2:4);

sprintf("Value for first A after convergence")
disp(A_store)
sprintf("Value for first L after convergence")
disp(l_store)

% Now solving for the DOPs latitude and longitude of each point
%
% r = sqrt(x0^2 + y0^2 + z0^2)
% theta = tan-1(y/x)
% phi = arccos(z/r)

% Solving for R, lambda and phi at each coordinate
for ii = 1:T
    nan_entries = sum(isnan(X(ii,:)));

    if nan_entries == 0
        r = norm(X(ii,:));
        lambda = atand(X(ii,2)/X(ii,1));  % in degrees
        phi = acosd(X(ii,3)/r);        % in degrees

        r11 = -sind(phi)*cosd(lambda);
        r12 = -sind(phi)*sind(lambda);
        r13 = cosd(phi);

        r21 = -sind(lambda);
        r22 = cosd(lambda);
        r23 = 0;

        r31 = cosd(phi)*cosd(lambda);
        r32 = cosd(phi)*sind(lambda);
        r33 = sin(phi);

        R(ii,:,:) = [r11 r12 r13;...
                    r21 r22 r23;...
                    r31 r32 r33];
    else
        r = nan;
        lambda = nan;
        phi = nan;
    end

    R_temp = reshape(R(ii,:,:),[3,3]);
    C_xyz_temp = reshape(Cxyz(ii,:,:),[3,3]);
    R_temp_T = reshape(R(ii,:,:),[3,3]);
    R_temp_T = R_temp_T';

    Cl_temp = R_temp*C_xyz_temp*R_temp_T;

    Cl(ii,:,:) = Cl_temp;

    % Defining DOPs values
    Cl_diag = diag(Cl_temp);
    sigma_n = Cl_diag(1);      % Each value is squared
    sigma_e = Cl_diag(2);
    sigma_u = Cl_diag(3);
    sigma_t = Cxs(ii,4,4);


    % Determining values for DOPs
    vdops(ii,1) = sqrt(sigma_u);
    hdops(ii,1) = sqrt(sigma_n + sigma_e);
    pdops(ii,1) = sqrt(sigma_n + sigma_e + sigma_u);
    tdops(ii,1) = sqrt(sigma_t);
    gdops(ii,1) = sqrt(sigma_n + sigma_e + sigma_u + sigma_t);
end



