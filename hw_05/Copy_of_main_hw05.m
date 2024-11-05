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
x_prev = [0 0 0 nan];

%first step only C1
C1 = obs.C1;
number_of_timepoints = length(C1);
number_of_satellites = width(C1);

result_matrix = zeros(number_of_timepoints, 4);
%later: compare solutions using C1, P1, and P2

for timepoint = 1:1%number_of_timepoints
    row_C1 = C1(timepoint,:);
    svs_at_timepoint = find(~isnan(row_C1));
    number_svs_at_timepoint = length(svs_at_timepoint);
    if(number_svs_at_timepoint < 5)
        disp('too few satellites at this timepoint to solve for rx position and clock bias')
        result_matrix(timepoint, : ) = nan(1, 4);
        disp('Continue with next timepoint!')
    else
        iteration = 1;
        while(~(abs(xyz_del_0(1:3) - x_prev(1:3)) < 1))
            % disp('Accuracy of position converged to less than 1m!')
            % disp(x_0);
            % disp('Continue with next timepoint!')
            % %hop to next timepoint
            % break
        % else
            disp('Keep improving rx position and clock bias for this timepoint!')
            A_1 = zeros(number_svs_at_timepoint, number_unknowns-1);
            %prevent numerical instabilities: use c *10^-9: get rx clock
            %bias in ns! - > instabilities got even worse
            c_vector = 1e-9*c.*ones(number_svs_at_timepoint, 1); 
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
                
                % %calculate elements of A for this satellite:
                % %jaX
                A(row_of_A, 1:3) = (- residual./rho_0)';
                l(row_of_A) =   row_C1(sv) - rho_0 + (c * satpos(4));
                % %jl
                % A(row_of_A, 4) = ;
                row_of_A = row_of_A + 1;

            end

            %solve for array of four unknowns
            x_prev = xyz_del_0;
            xyz_del_0 = (pinv(A)*l)';
            
            if (iteration == 50)
                disp("The solution for the receiver position doesn't converge at this timepoint")
                disp("difference between this and previous position:")
                disp(abs(xyz_del_0(1:3) - x_prev(1:3)))
                result_matrix(timepoint, : ) = nan(1, 4);
                break
            end
            iteration = iteration + 1;
        end
        if(abs(xyz_del_0(1:3) - x_prev(1:3)) < 1)
            disp("solution converged to")
            disp(xyz_del_0)
            disp("difference between this and previous position:")
            disp(abs(xyz_del_0(1:3) - x_prev(1:3)))
            result_matrix(timepoint, : ) = xyz_del_0;
        end
    end
    
end

