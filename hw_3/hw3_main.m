clear all
close all
clc

task = 2.5;

[eph,alpha,beta,dow] = read_rinexn('brdc2920.19n');

sv=31;
if task==1
    %  % t = 6;
    %  % t = t*3600*24; %518400
    t = 6*3600*24 + (24*3600 - 15*60);%603900
    rho = getsatpos(t, sv, eph);


elseif task ==2
    %2)
    % first with for loop than vectorized!
    % times = linspace(0,24*3600-1,96)';
    % t = [6;7];
    t = 6*3600*24;
    rho = getsatpos(t, sv, eph)
    
    % begin_epoch_interest = 518400;
    % end_epoch_interest = 603900;
    % time_del_steps = 15*60;
    % entries = (end_epoch_interest-begin_epoch_interest)/time_del_steps;
    % rho = zeros(3, entries);
    % i = 1;
    % for time = begin_epoch_interest:time_del_steps:end_epoch_interest 
    %     rho(:,i) = getsatpos(time, sv, eph);
    %     i = i+1;
    % end
    % save('rho.mat','rho')


elseif task == 2.5
    load('rho.mat');
    [sp3,sv_sp3,bad_sat] = read_sp3('igs20756.sp3');
    reference_data = sp3.prn31;
    rho = rho';
    residual = sqrt((rho(:,1)-reference_data(:,2)).^2 ...
        +(rho(:,2)-reference_data(:,3)).^2 ...
        + (rho(:,3)-reference_data(:,4)).^2);
    plot(residual, 'LineWidth', 2);
    old_labels = xticklabels;
    num_label = size(xticklabels);
    num_label = num_label(1);
    new_labels = strings(1,num_label);
    for ii= 1:num_label
        new_label = round(0.25*str2num(old_labels{ii,1}),1);
        new_labels(ii) = num2str(new_label);

    end
    % new_labels = xticklabels*0.25
    xticklabels(new_labels);
    title("Comparison of ECEF coordinates by navigation message and from sp3 file for SV31");
    xlabel("time of the day [hours]");
    ylabel("Norm of 3D residual position [m]");
    grid on
    saveas(gcf, 'comparison_ECEF_rho_sp3.png'); 
    absolute_diff_X = abs(rho(:,1)-reference_data(:,2));
    absolute_diff_Y = abs(rho(:,2)-reference_data(:,3));
    absolute_diff_Z = abs(rho(:,3)-reference_data(:,4));
    % figure
    % plot(absolute_diff_X)
    % figure
    % plot(absolute_diff_Y)
    % figure
    % plot(absolute_diff_Z)
end

% Taks 3)
% The error, difference in a range of 350m is 
% a facto 100 greater than the expected one 
% (than it should be). There must be a mistake in 
% my calculation of rho (probably in one variable)
% that causes these additional 100m offset in all directions.