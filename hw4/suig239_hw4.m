clearvars
close all
clc
[obs,t,prn,apr,hant] = read_rinexo('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/opmt2920.19o');

% Output: obs = structure with observation types, e.g. obs.C1, obs.L1, etc.
%         t = time matrix (year month day hour minute second)
%         gps = vector of GPS satellite numbers
%         apr = approximate site position from rinex file header
%% 4 
 % add up all time 
time = datetime(t(:,1)+2000,t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
ref_time = datetime(1900,1,1,0,0,0); 
rel_time = seconds(time - ref_time );
f1 = 1.57542 * 10^9; %Hz
f2 = 1.2276 * 10^9; %Hz
c = 0.299792458 * 10 ^9; % m/s
lambda1 = c/f1;
lambda2 = c/f2;
alpha = (f1/f2)^2;

sv = 10;
prn_index = find(prn==sv)

P1 = obs.P1(:,prn_index);
P2 = obs.P2(:,prn_index) ;
L1 = obs.L1(:,prn_index) 
L2 = obs.L2(:,prn_index);

non_nan_index = find(~isnan(L1));
non_nan_index = non_nan_index(1:600);
L1 = L1(non_nan_index)*lambda1;
L2 = L2(non_nan_index) *lambda2;
P1 = P1(non_nan_index);
P2 = P2(non_nan_index);

time = time(non_nan_index); 
time_sec = rel_time(non_nan_index);

MP1 = P1 - ((2/(alpha-1)) + 1)*L1+(2/(alpha-1))*L2;
MP2 = P2 - ((2*alpha/(alpha-1)))*L1+((2*alpha/(alpha-1))-1)*L2;

MP1mean=mean(MP1(~isnan(MP1)));
MP2mean = mean(MP2(~isnan(MP2)));

MP1_anom = MP1-MP1mean;
MP2_anom = MP2-MP2mean;

figure()
subplot(2,1,1)
plot(time,MP1_anom,'b');
title("MP1 for PRN" + num2str(sv));
xlabel("Time (hours)")
ylabel("MP1 (m)")
grid();
datetick;
subplot(2,1,2)
plot(time,MP2_anom,'r');
xlabel("Time (hours)")
ylabel("MP2 (m)")
grid()
title("MP2 for PRN" + num2str(sv))
datetick;
%%

% discussing reults  

%{

6) To get these plots I removed the NaN values to get a section of data that
had continous data (around 4:00 to 9:15). There is greater noise showing at
the beginning and end of MP1 and MP2, suggesting that there is some interference occuring around
these times. Otherwise, the values remained between -1 and 1 meters with
little interference, posssibly indicating that the satellite was directly
over the antenna. 

7) An additional noise term that could also be a factor of additional noise
is the lack of ionosphere and tropospheric models. This could cause some additional noise since we did not factor either of those into the 
equations (although the ionosphere was estimated in this case). 


%}