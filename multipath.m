% Homework 4
% Compute multipath
% Rubi Garcia      OCT 2024
%clc; clear all

%Variables
f1 = 1.57542*10^9%Hz 
f2 = 1.2276*10^9 %Hz
c = 0.299792458*10^9 %m/s
lambda1 = c/f1
lambda2 = c/f2
alpha = (f1/f2)^2
 
%[obs,t,gps,apr,hant] = read_rinexo('./../hw3/products/opmt2920.19o')

time_n = datetime(2000+t(:,1),t(:,2),t(:,3),t(:,4),t(:,5),t(:,6))

%For sv 10
sv = 10
index = find(gps==sv)

L1 = obs.L1(:,index)*lambda1;
L1_ind = find(L1>0)
L1_ind2 = find(L1<0)
L1_indexes = [L1_ind2;L1_ind]
L1 = L1(~isnan(L1));
L2 = obs.L2(:,index)*lambda2;
L2 = L2(~isnan(L2))
L2 = L2(1:851)
C1 = obs.C1(L1_indexes,index);
P1 = obs.P1(L1_indexes,index);
P2 = obs.P2(L1_indexes,index);
S1 = obs.S1(L1_indexes,index);
S2 = obs.S2(L1_indexes,index);

%For sv 3
% sv = 3
% index = find(gps==sv)
% 
% L1 = obs.L1(:,index)*lambda1;
% L1 = L1(~isnan(L1));
% L2 = obs.L2(:,index)*lambda2;
% L2 = L2(~isnan(L2))
% C1 = obs.C1(L1_indexes,index);
% P1 = obs.P1(L1_indexes,index);
% P2 = obs.P2(L1_indexes,index);
% S1 = obs.S1(L1_indexes,index);
% S2 = obs.S2(L1_indexes,index);

% %Compute MP1 Multipath for the L1 longwavelength
MP1 = P1 - ((2/(alpha-1)) + 1)*L1+(2/(alpha-1))*L2
MP2 = P2 - ((2*alpha/(alpha-1)))*L1+((2*alpha/(alpha-1))-1)*L2

MP1_mean = mean(MP1(~isnan(MP1)))
MP2_mean = mean(MP2(~isnan(MP2)))

% MP1_mmean = (MP1(~isnan(MP1))) - MP1_mean;
MP1_mmean = MP1 - MP1_mean;
MP2_mmean = MP2 - MP2_mean;

%Time for sv 3
time_n_prn10 = time_n(L1_indexes)


subplot(2, 1, 1);
plot(time_n_prn10, MP1_mmean, 'b', 'LineWidth', 2);
title('L1 wavelength');
xlabel('hour');
ylabel('m');
ylim([-400000,400000])
grid on;
set(gca, 'FontSize', 20)
set(gcf, 'Color', 'w')
box off

subplot(2, 1, 2);
plot(time_n_prn10, MP2_mmean, 'r', 'LineWidth', 2);
title('L2 wavelength');
xlabel('hour');
ylabel('m');
ylim([-400000,400000])
grid on;
set(gca, 'FontSize', 20)
set(gcf, 'Color', 'w')
box off

sgtitle('Multipath for G10','FontSize',20);