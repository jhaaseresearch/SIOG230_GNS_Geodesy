clearvars
close all
clc

% eph = read_rinexn('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/brdc2910.19n');
[eph,alpha,beta,dow] = read_rinexn('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/epgga9.292');

 %% 1 a-f 
t = 6 * 24 * 3600 + 24 * 3600 - 15 * 60;
sv =  31
rho = get_satpos(t, sv, eph);
fprintf('Satellite position in orbital coordinates (X, Y, Z): %.3f, %.3f, %.3f\n', rho(1), rho(2), rho(3));


%% 2) 
%comparison to sp3 file
t = 6*3600*24
sv = 31;
rho = get_satpos(t, sv, eph)
[sp3, sv, excl_sat]= read_sp3('/Users/hlowesbicay/Classes/SIOG239_2024/SIOG230_GNS_Geodesy/products/igs20756.sp3');

prn1_data = sp3.prn1;
rho = rho'; 

residual = sqrt((rho(:,1)-prn1_data(:,2)).^2 +(rho(:,2)-prn1_data(:,3)).^2 + (rho(:,3)-prn1_data(:,4)).^2);

plot(residual, 'LineWidth', 2);

title("Residuals of ECEF coordinates of broadcast orbit and sp3 for SV31");
xlabel("time of day (s)");
ylabel("Residual position (m)");
abs_diff_X = abs(rho(:,1)-prn1_data(:,2));
abs_diff_Y = abs(rho(:,2)-prn1_data(:,3));
abs_diff_Z = abs(rho(:,3)-prn1_data(:,4));

% figure
% plot(abs_diff_X)
% figure
% plot(abs_diff_Y)
% figure
% plot(abs_diff_Z)
%% 3 
% I'm not sure why my error is so huge, there is definitely something wrong
% in my get_satpos function that's calculating the orbit at the wrong
% position but I'm having trouble finding where. Especially at 50 seconds. 


% Also it looks like my residuals should look like oscillations instead of
% a curve. 


% It looks like the greatest difference is at 50 seconds. The error
% increases and has a miximum at 50 s then decreases afterwards. 

% My error is way too big. GOing to submit this for now and keep tryring to
% debug my get_satpos function. 