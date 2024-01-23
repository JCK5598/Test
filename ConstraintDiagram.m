clear all; close all; clc; 
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultAxesFontSize',20);
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLegendFontSize',18,'DefaultLegendFontSizeMode','manual');
set(groot,'defaultfigureposition',[300 150 800 500])
 
rho = 0.002377; % density at sea level
g = 32.17; % gravity
 
 
v_stall_guess = 35; % ft/s
CL_max = 1.5;
dTO = 10;
n_max = 1.75;
AR = 2.781;
e0 = 0.748;
k = 1 / (pi * AR * e0);
v_cruise = 55; % ft/s
q = 0.5*rho*v_cruise^2;
CD0 = 0.02;
 
W_S = 0:0.05:5.0;
 
 
WS_max_1 = 0.5*rho*v_stall_guess^2*CL_max;
WS_max_2 = 0.5*rho*(v_stall_guess*1.1)^2*(CL_max);
WS_max_3 = 0.5*rho*(v_stall_guess*0.9)^2*(CL_max);
 
T_W_TO_1 = 1.21*W_S/(rho*g*CL_max*dTO);
T_W_TO_2 = 1.21*W_S/(rho*g*CL_max*1.1*dTO);
T_W_TO_3 = 1.21*W_S/(rho*g*CL_max*0.9*dTO);
 
T_W_Turn_1 = n_max^2*k*W_S./q + q*CD0./W_S;
 
figure(2)
hold on
plot([WS_max_1 WS_max_1], [0 3])
% plot([WS_max_2 WS_max_2], [0 3])
% plot([WS_max_3 WS_max_3], [0 3])
plot(W_S,T_W_TO_1)
% plot(W_S,T_W_TO_2)
% plot(W_S,T_W_TO_3)
plot(W_S,T_W_Turn_1)
xlabel('$W/S (lb/ft^2)$')
ylabel('$T/W$')
% xlim([0 2.0])
ylim([0 3])
legend('Stall Speed','Takeoff','Turn','location','southeast')
 

B = 20;
W_empty = 0.1:0.1:5;
W_payload = 4.2;
W = W_empty + W_payload;
S = W ./ 1.9;

L_Dmax = 1 / (2 * sqrt(k * CD0));

CL_max_calc = sqrt(3 * CD0 / k);

% FS = (3 .* W_payload .* 11 ./ ((W_empty - 1).^4 + 8.9)) + B - S.^1.5;
% 
% figure
% plot(W, FS);

W_TO_S = 0:0.01:5;

V_stall = sqrt((2 .* W_TO_S ./ (0.970075 * 0.002377 * 1.5)));

V_TO = 1.2 * V_stall;

V_min = sqrt(2 .* W_TO_S ./ (0.002377 * CL_max_calc));

T_W = (1 / 0.970075) * ((q * CD0 ./ W_TO_S) + (n_max^2 .* k .* W_TO_S ./ q));

% takeoff in 10 feet
P_W_takeoff = (0.7 .* V_stall).^3 ./ (2 * 550 * 0.970075 * 0.75^2 * 10 * 32.17);

% landing in under 200 ft
P_W_landing = (V_stall.^3) ./ (550 * 0.970075 * 0.7^2 * 200 * 32.17);

% rate of climb at 5 ft/s
P_W_rotc = (1 / (550 * 0.970075 * 0.7^2)) .* (5 + (V_min ./ (0.866 * L_Dmax * sqrt(0.970075))));

% handle 4g turn
P_W_turn = v_cruise .* T_W ./ (550 * 0.7^2);

figure(1)
plot(W_TO_S, P_W_takeoff)
hold on
plot([WS_max_1 WS_max_1], [0 1])
plot(W_TO_S, P_W_landing)
plot(W_TO_S, P_W_rotc)
plot(W_TO_S, P_W_turn)
plot(7/3.455, 0.0862, 'o', 'LineWidth', 1, 'Color', [0, 0, 0])

% Some power estimation divided by MTOW
% plot([0 5], [0.6/7  0.6/7], '--')
% plot([0 5], [0.6/8  0.6/8], '-o')
% plot([0 5], [0.6/9  0.6/9], ':')

ylim([0 0.25])
xlim([0 3])
xlabel('$W/S (lb/ft^2)$')
ylabel('$P/W (hp/lb)$')
legend('Takeoff', 'Stall Speed', 'Landing', 'Rate of Climb', 'Turn', 'Design Point', 'location', 'northeast')
title('Constraint Diagram')




% payload = (1:4:13) * 4.448;
% laps = 0:1:15;
% R = laps*900;
% g = 9.8;
% LD = 3.5;
% eta = 0.65;
% e = 300000;
% WeWTO = 0.25;
% Max_M2 = laps(end)*payload(end);
% % Max_M2 = 200*4.448; % 10 laps with 20 lbs
% figure
% hold on
% for p = 1:length(payload)
% %     for L = 1:length(laps)
%         WbWTO = g*R/(e*LD*eta)+0.02;
%         WTO = payload(p)./(1 - WeWTO - WbWTO);
%         M2 = 1 + payload(p).*laps./Max_M2;
% %     end
%         plot(laps, M2)
% end
% xlabel('Laps')
% ylabel('Mission 2 Score')
% % Assuming a maximum weight to keep things reasonable
% for p = 1:length(payload)
%     for L = 1:length(laps)
%         WbWTO = g*R/(e*LD*eta);
%         payload_10=10*4.48*(1 - WeWTO - WbWTO);
%         M2_10 = 1 + payload_10.*laps./Max_M2;
%      end
%         plot(laps, M2_10, '--')
% 
%         payload_20=20*4.48*(1 - WeWTO - WbWTO);
%         M2_20 = 1 + payload_20.*laps./Max_M2;
% %      end
%         plot(laps, M2_20, '-.')
% 
%         payload_12=12*4.48*(1 - WeWTO - WbWTO);
%         M2_12 = 1 + payload_12.*laps./Max_M2;
% %      end
%         plot(laps, M2_12, ':')
% end


