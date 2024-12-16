% Loading observation file via Jennifer's function:
[obstypes, t, gpssat, sitepos, hant] = read_rinexo('opmt2920.19o');
[eph,alpha,beta,dow] = read_rinexn('epgga9.292')
% Initializing needed constants:
f1 = 1.57542 * (10^9); 
f2 = 1.2276 * (10^9) ;
c = 0.299792458 * 10^9; 
lam1 = c/f1;
lam2 = c/f2;
alpha = (f1/f2)^2;



% For PRN3 to check my equations
% Get necessary data from the observation types structure:
L1 = obstypes.L1;
L1PRN3 = L1(:,3)* lam1;
L2 = obstypes.L2;
L2PRN3 = L2(:,3) * lam2;
P1 = obstypes.P1;
P1PRN3 = P1(:,3);
P2 = obstypes.P2;
P2PRN3 = P2(:,3);
%Time in seconds array:
timesecs = 0:30:(30* length(L1PRN3) -1);
% Calculate MP1 and MP2:
constmp1 = 2/(alpha - 1);
constmp2 = (2 * alpha)/ (alpha - 1);
MP1 = (P1PRN3 - ((constmp1 +1)*L1PRN3) + (constmp1 * L2PRN3));
MP1 = MP1 - mean(MP1(~isnan(MP1)));
MP2 = P2PRN3 - (constmp2 * L1PRN3) + ((constmp2 - 1)*L2PRN3);
MP2 = MP2 - mean(MP2(~isnan(MP2)));

% Plot MP1 and MP2:
figure(1);
clf;
hold on;
sgtitle('PRN 3')
subplot(2,1,1)
plot( MP1, '-', linewidth=2)
ylabel('MP1 Error (m)')
xlabel('Time (in seconds past midnight on 10/19/19)')
subplot(2,1,2)
plot(MP2, '-', linewidth = 2)
ylabel('MP2 Error (m)')
xlabel('Time (in seconds past midnight on 10/19/19)')

% For PRN10 to check my equations
% Get necessary data from the observation types structure:
xind1 = 2 * 60 * 4;
xind2 = 120 * 9;
L1PRN10 = L1(:,9)* lam1;
L2PRN10 = L2(:,9) * lam2;
P1PRN10 = P1(:,9);
P2PRN10 = P2(:,9);


% Calculate MP1 and MP2:
MP110 = (P1PRN10 - ((constmp1 +1)*L1PRN10) + (constmp1 * L2PRN10));
% Took this process from Lani
threshold = 1; 
Group_ID = cumsum([1; diff(MP110 > threshold) ~= 0]); 
[G, ~] = findgroups(Group_ID); 
GroupMeans = groupsummary(MP110, G, 'mean'); 
MP110 = MP110 - GroupMeans(G)
%
MP210 = P2PRN10 - (constmp2 * L1PRN10) + ((constmp2 - 1)*L2PRN10);
% Took this process from Lani
threshold = 1; 
Group_ID = cumsum([1; diff(MP210 > threshold) ~= 0]);
[G, ~] = findgroups(Group_ID); 
GroupMeans = groupsummary(MP210, G, 'mean');
%
MP210 = MP210 - GroupMeans(G)

% Plot MP1 and MP2:
figure(2);
clf;
hold on;
sgtitle('PRN 10')
subplot(2,1,1)
plot(MP110, '-', linewidth=2)
ylabel('MP1 Error (m)')
xlabel('Time (in seconds past midnight on 10/19/19)')
subplot(2,1,2)
plot(MP210, '-', linewidth = 2)
ylabel('MP2 Error (m)')
xlabel('Time (in seconds past midnight on 10/19/19)')

% From the results of PRN10, I can see that the satelite gives good data
% from approximately 500 to 1200 seconds past midnight on October 19, 2019.
% However, after that time the receiver loses visability of the satelite in
% the sky from 1200 seconds to 2400 seconds. I know that the receiver at
% this time is still online from PRN3 data. However, the data after 1800
% seconds gets noisier for PRN 3 and from 2400 from PRN 10(with a small gap in data around 2500 seconds on the PRN 10 data). 

% I think that the noise level could be due to the fact that tropospheric
% delay was NOT included in my calculations.