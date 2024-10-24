%Homework 04

%Open Rinex file
[obs,t,prn,apr,hant] =read_rinexo("opmt2920.19o")

time_vector = datetime(t(:,1), t(:,2), t(:,3), t(:,4), t(:,5), t(:,6));
%Compute and plot MP1 and MP2 for PRN3 and PRN10 
%Compute MP1 and MP2
f1=1.57542*10^9
f2=1.2276*10^9
c=0.299792458*10^9
lambda1=c/f1
lambda2=c/f2
A=40.3

alpha=(f1/f2)^2

%Select the satellite index
sv=3
prn_i = find(prn == sv)

L1=obs.L1(:,prn_i)*lambda1
L2=obs.L2(:,prn_i)*lambda2
P1=obs.P1(:,prn_i)
P2=obs.P2(:,prn_i)
MP1_all=P1 - ((2/(alpha - 1))+1) * L1 + (2/(alpha - 1)) * L2
MP2_all=P2 - ((2 * alpha)/(alpha - 1)) * L1 + (((2 * alpha)/(alpha - 1))-1) * L2

MP1_mean=mean(MP1_all(~isnan(MP1_all)))
MP2_mean=mean(MP2_all(~isnan(MP2_all)))

MP1 = MP1_all - MP1_mean
MP2 = MP2_all - MP2_mean

subplot(4, 1, 1);  
plot(time_vector, MP1, 'b-', 'LineWidth', 2);
title('MP1 PRN3');
xlabel('Time (s)');
ylabel('MP1');
grid on;

subplot(4, 1, 2);  
plot(time_vector, MP2, 'r-', 'LineWidth', 2);
title('MP2 PRN3');
xlabel('Time (s)');
ylabel('MP2');
grid on;

%Select the satellite index
sv=10
prn_i = find(prn == sv)

L1=obs.L1(:,prn_i)*lambda1
L2=obs.L2(:,prn_i)*lambda2
P1=obs.P1(:,prn_i)
P2=obs.P2(:,prn_i)
MP1_all=P1 - ((2/(alpha - 1))+1) * L1 + (2/(alpha - 1)) * L2
MP2_all=P2 - ((2 * alpha)/(alpha - 1)) * L1 + (((2 * alpha)/(alpha - 1))-1) * L2

threshold = 1; 
Group_ID = cumsum([1; diff(MP1_all > threshold) ~= 0]); 
[G, ~] = findgroups(Group_ID); 
GroupMeans = groupsummary(MP1_all, G, 'mean'); 
MP1 = MP1_all - GroupMeans(G)

threshold = 1; 
Group_ID = cumsum([1; diff(MP2_all > threshold) ~= 0]);
[G, ~] = findgroups(Group_ID); 
GroupMeans = groupsummary(MP2_all, G, 'mean');
MP2 = MP2_all - GroupMeans(G)


subplot(4, 1, 3);  
plot(time_vector, MP1, 'b-', 'LineWidth', 2);
title('MP1 PRN10');
xlabel('Time (s)');
ylabel('MP1');
grid on;

subplot(4, 1, 4);  
plot(time_vector, MP2, 'r-', 'LineWidth', 2);
title('MP2 PRN10');
xlabel('Time (s)');
ylabel('MP2');
grid on;


%6-Discuss the results

%I can observe that the MP1 and MP2 can achieve higher values like 2-4
%meters and how they can vary on time. According to the Calais (2024), 
% the noise can achieve 1-10 m for range measurements 
% and less than 5 cm for phase measurements.
% Based on the PRN3 results which is more continuous than PRN10, 
%I observed that there is an interval time where the multipath error is
%lower. This may due the satellites movement, so the multipath noise increase as the 
% satellite zenith angle increases. 

%7-Isn’t there an important additional noise term that is not modeled here and which, conse-
%quently, is present in MP1 and MP2 as derived above?

%We did not use a model to correct the troposphere. 
%Moreover, we used the alpha parameter (based on the L1 and L2 frequency) instead of 
%using a model to correct the ionospheric delay. 
%