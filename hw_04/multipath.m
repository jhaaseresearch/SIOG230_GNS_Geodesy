clear all
close all
clc

%constants
f1 = 1.57542 *109;% Hz
f2 = 1.2276 * 109;% Hz
c = 0.299792458 * 109;% m/s 
lambda_1 = c/f1;
lambda_2 = c/f2;

[data,t,prn,apr,hant] = read_rinexo("opmt2920.19o");

sv = 10;
idx = find(prn==sv);
L1 = data.L1(:,idx);
L2 = data.L2(:,idx);
P1 = data.P1(:,idx);
P2 = data.P2(:,idx);
idx_n_nan=find(~isnan(L1));
k=2;
start_idx = 1;
while k < length(idx_n_nan)
    while idx_n_nan(k) == idx_n_nan(k-1)+1
        k=k+1;
        if k == length(idx_n_nan)
            break
        end
    end
    disp(k-1)
    end_idx = k-1;
    group_start = idx_n_nan(start_idx);
    group_end = idx_n_nan(end_idx);
    [MP1, MP2] = compute_multipath(f1,f2, lambda_1.*L1(group_start:group_end), lambda_2.*L2(group_start: group_end),P1(group_start:group_end), P2(group_start:group_end));
    subplot(2,1,1);
    plot(linspace(group_start, group_end, group_end - group_start+1),MP1 -mean(MP1));hold on
    
    subplot(2,1,2);
    plot(linspace(group_start, group_end, group_end - group_start+1),MP2 -mean(MP2));hold on
    
    group_start=[];
    group_2_start=[];
    MP1 = [];
    MP2 = [];
    k = k+1;
    start_idx = k;
end

%%L2 doesn't have the same NaNs as L1..why don't I get results for all the
%other groups (but group 1) for MP2??
% idx_n_nan=find(~isnan(L2));
% k=2;
% start_idx = 1;
% while k < length(idx_n_nan)
%     while idx_n_nan(k) == idx_n_nan(k-1)+1
%         k=k+1;
%         if k == length(idx_n_nan)
%             break
%         end
%     end
%     disp(k-1)
%     end_idx = k-1;
%     group_start = idx_n_nan(start_idx);
%     group_end = idx_n_nan(end_idx);
%     [MP1, MP2] = compute_multipath(f1,f2, lambda_1.*L1(group_start:group_end), lambda_2.*L2(group_start: group_end),P1(group_start:group_end), P2(group_start:group_end));
%     % subplot(2,1,1);
%     % plot(linspace(group_start, group_end, group_end - group_start+1),MP1 -mean(MP1));hold on
%     % 
%     subplot(2,1,2);
%     plot(linspace(group_start, group_end, group_end - group_start+1),MP2 -mean(MP2));hold on
% 
%     group_start=[];
%     group_2_start=[];
%     MP1 = [];
%     MP2 = [];
%     k = k+1;
%     start_idx = k;
% end

subplot(2,1,1);
title("MP1 for PRN" + num2str(sv));
xlim([0 idx_n_nan(end)]+200)
xlabel("Time")
ylabel("MP1 (m)")
xticklabels({"04:09:30", "08:19:30","12:29:30","16:39:30","20:49:30"});
grid();
subplot(2,1,2);
xlim([0 idx_n_nan(end)]+200)
xlabel("Time")
ylabel("MP2 (m)")
xticklabels({"04:09:30", "08:19:30","12:29:30","16:39:30","20:49:30"});
grid()
title("MP2 for PRN" + num2str(sv))
savefig("multipath_results.fig")

% labels=xticklabels;
% hour_min = zeros(length(labels), 2);
% for i = 1:length(labels)
%     i
%     hour_min(i,:) = t(1+str2num(cell2mat(labels(i))), 4:5)
% end

%% 6) Discuss the results
%
% For PRN 10, I get results for four groups of "consecutive not NaNs" that 
% I could calculate MP1 for. I substracted the mean of each group
% individually before plotting the multipath result (subplot 1).
% Times that involved NaNs on L1 are skipped. I use L1 to find NaNs (the NaNs
% in L1 and L2 are not the same, but when using L2 to find them I only get
% a plot for the first group for MP1 and MP2.) This way I can only plot one
% group of MP2 (corresponding to the first time intervall of "consecutive 
% not NaNs") from 04:09:00 UTC to 09:21:00 UTC on 10/19/2019 (blue in both
% subplots). Comparing this group, it can be seen that  MP1 and MP2 vary in
% values at the start of the recording. This could be a due to more interference.
% Over time, the range of MP gets smaller (in between -1 and 1m), and then 
% become greater again (towards the end of the recording of the "blue" group).
% This means that for both (MP1 and MP2) the multipath gets less and then
% bigger again. This corresponds to the "typical" cycle of "overflight"
% by a GNSS satellite. While it is rising and setting its signal is more affected
% by multipath effects in contrast to the time when it is "straight" above
% the receiver in the sky (LOS), which corresponds to the middle of the 
% "blue group". This cycle of one overflight is seen in MP1 and MP2.
% Subplot 1 furthermore shows the diurnal cylce of a GNSS satellite which
% is "in view" of a receiver (of fixed location) twice a day. Considering
% the "three plots in the evening" as one group, the overall "cycle for MP
% in one overflight" can also be seen (MP big at the start and the end, 
% smaller in between). But it can also be seen that within this cycle there
% was a complete loss of signel (gap between yellow and purple plot). It 
% should be noticed that before the signal is lost MP1 keeps increasing 
% (end of yellow plot). Comparing the "end of the yellow plot" with the
% "end of the blue plot", a similar value of MP1 (~2m) can be seen.
% The re-acquirement of the signal (begin of purple plot) shows a very
% narrow range of MP1, which corresponds to a satellite "which is still
% in rather high elevation". If the satellite were already lower, it would
% be less likely that the re-acquirement work, since the signal would 
% already be more affected by multipath.

%% 7) 
% Additional noise terms that I did not correct for (assumed gaussian noise
% when "just substracting the mean") were caused be the Ionosphere and 
% Troposphere.