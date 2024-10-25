% clear all
% close all
% clc
% 
% %constants
% f1 = 1.57542 *109;% Hz
% f2 = 1.2276 * 109;% Hz
% c = 0.299792458 * 109;% m/s 
% lambda_1 = c/f1;
% lambda_2 = c/f2;
% 
% [data,t,prn,apr,hant] = read_rinexo("opmt2920.19o");

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
xlabel("Time")
ylabel("MP1 (m)")
grid();
subplot(2,1,2);
xlim([0 idx_n_nan(end)])
xlabel("Time")
ylabel("MP2 (m)")
grid()
title("MP2 for PRN" + num2str(sv))
savefig("multipath_results.fig")
