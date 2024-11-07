close all
clear all
clc

%Define constants (c)
c = 2.99792458e8;

% define satellite elevation cutoffangle (10◦ for instance):
%not yet

%Reading ephemeris (broadcast data): EPH
[eph,alpha,beta,dow] = read_rinexn("brdc2920.19n");

%Reading observation data: OBS
[obs,t,gps,apr,hant] = read_rinexo("opmt2920.19o");

%find the time in seconds:
time_day = datetime(t(:,1)+2000,t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
time_ref = datetime(2019,10,19,0,0,0);
time_s = seconds(time_day - time_ref);


%The a priori position of the receiver in ECEF frame (in meters) is 
% provided in the rinex file header (APPROX POSITION XYZ) + t_0
rx_clock_bias_0 = 1;
xyz_del_0 = [apr rx_clock_bias_0];
number_unknowns = length(xyz_del_0);
x_prev = [0 0 0];

%first step only C1
C1 = obs.C1;
number_of_timepoints = length(C1);
number_of_satellites = width(C1);

result_matrix = zeros(number_of_timepoints, 4);
%later: compare solutions using C1, P1, and P2

for timepoint = 1:number_of_timepoints
    row_C1 = C1(timepoint,:);
    svs_at_timepoint = find(~isnan(row_C1));
    number_svs_at_timepoint = length(svs_at_timepoint);
    if(number_svs_at_timepoint < 5)
        disp('too few satellites at this timepoint to solve for x')
        result_matrix(timepoint, : ) = nan(1, 4);
        disp('Continue with next timepoint!')
    else
        if ((xyz_del_0(1:3) - x_prev) < 1)
            disp('Accuracy of position converged to less than 1m!')
            disp(x_0);
            disp('Continue with next timepoint!')
            %hop to next timepoint
            break
        else
            disp('Keep improving x for this timepoint!')
            A_1 = zeros(number_svs_at_timepoint, number_unknowns-1);
            c_vector = c.*ones(number_svs_at_timepoint, 1);
            A = [A_1 c_vector];
            l = zeros(number_svs_at_timepoint, 1);
            row_of_A = 1;
            for sv = svs_at_timepoint
                %pseudorange data tagged time = t in the rinex observation file were sent 
                %a bit earlier, at ttx = t − j P /c: 
                % "time(s) - pseudorange("that timepoint", j^th satellite)/c"
                ttx = time_s(timepoint) - row_C1(sv)/c;

                %get_satpos (from Jennifer: with sat clock bias dt,toe,drel,tgd)
                satpos = get_satpos(ttx, sv, eph,3);
                %calculate residual jX - X_0
                residual = satpos(1:3) - xyz_del_0(1:3)';
                %calculate rho_0
                rho_0 = sqrt(sum((residual).^2));
                disp('test')
                
                % %calculate elements of A for this satellite:
                % %jaX
                A(row_of_A, 1:3) = (- residual./rho_0)';
                l(row_of_A) =   row_C1(sv) - rho_0 + (c * satpos(4));
                % %jl
                % A(row_of_A, 4) = ;
                row_of_A = row_of_A + 1;

            end

            %solve for array of four unknowns
            new = pinv(A)*l;
            disp('tst')

        end

    end
    
end

