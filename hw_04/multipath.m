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

idx = find(prn==10);
L1 = data.L1(:,idx);
L2 = data.L2(:,idx);
P1 = data.P1(:,idx);
P2 = data.P2(:,idx);
idx_n_nan=find(~isnan(L1))
k=2;
while idx_n_nan(k) == idx_n_nan(k-1)+1
    k=k+1;
end
disp(k-1)
group_1_start = idx_n_nan(1);
group_1_end = idx_n_nan(k-1);

[MP1, MP2] = compute_multipath(f1,f2, lambda_1.*L1(group_1_start:group_1_end), lambda_2.*L2(group_1_start: group_1_end),P1(group_1_start:group_1_end), P2(group_1_start:group_1_end));
subplot(2,1,1);
plot(MP1 -mean(MP1));
subplot(2,1,2);
plot(MP2 -mean(MP2));