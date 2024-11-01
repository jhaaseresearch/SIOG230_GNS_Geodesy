%{
 HW5

Author: Tommy Stone
Class: SIOG 239G
Date: Fall 2024
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

% Unsure about where the data is being read for for the Satellite
% clock bias. Currently getting an error
%(\___/)
%(=^.^=) current issue here
%(")_(")
brd_file = 'igs20756.sp3';

%% 2. Read the pseudorange data 
[obs,ts,gps,apr,hant] = read_rinexo(rixen_file);

% filtering eph data
[eph, alpha, beta, dow] = read_rinexn(brd_file);

C1 = obs.C1;

% initial guess for positions
x0 = apr(1);
y0 = apr(2);
z0 = apr(3);

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

% iterate overtime to make calculations
start = 1;
for ii = start:start + 2%T
    ttx_temp = ttx(ii,:);
    
    % Find positions without nan values in ttx_temp
    % This should have at least 4 non Nan values
    idx_nn = find(~isnan(ttx_temp));
    
    % number of valid values
    num_real = length(idx_nn);
    
    % If value is greater or equal to 4 then we can solve, if not we will
    % skip the record if not
    if num_real >= 4
        % determine the Satellites to use
        svs = prns(idx_nn);

        % initialize matrix A and vector l
        A = nan([num_real,4]);

        l = nan([num_real,1]);
        
        % get the time component
        t_real = ttx(ii,idx_nn);

        C1_real = C1(ii,idx_nn);

        % tolerance should be within 1m
        tol = 1;
        error = 5;
        itr = 1;

        % using least squares to solve. 
        while (tol < error && itr <= 20)

            itr = itr + 1;

            % iterate through each Satellite to get the position
            for idx_sv = 1:length(svs)
                sv = svs(idx_sv);
                t = t_real(idx_sv);
                
                % [X Y Z dt toe drel]
                satpos = get_satpos(t,sv,eph,3);
                xj = satpos(1);
                yj = satpos(2);
                zj = satpos(3);
                deltaj = satpos(4);       % Satellite clock bias
    
                % Determining constants required to build matrix for each
                % Satellite
                rho0j = sqrt((xj - x0).^2 +(yj - y0).^2 +(zj - z0).^2);
    
                C1j = C1_real(idx_sv);
    
                % Calculate l for Satellite j
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
    
            % matrix A and l have been built, solving for x
            AT = A';
            X = pinv(AT*A)*AT*l;
    
            % Determine tolerance
            error = sqrt((x0 - X(1))^2 + (y0 - X(1))^2 + (z0 - X(3))^2);
    
            % Assigning new values
            x0 = X(1);
            y0 = X(2);
            z0 = X(3);

            disp(X)
        end

    end
end

