clc
clear
close all

%% User parameters
testNum = 3;
anomaly = true;
calib = true;
range = 'all'; % all, beg, mid, rest

%% Hard coded values
% x limits
if testNum == 2
    if strcmp(range, 'all')
        xlim_2 = [0 35];
    elseif strcmp(range, 'beg')
        xlim_2 = [0 15.5];
    elseif strcmp(range, 'mid')
        xlim_2 = [16 28.5];
    elseif strcmp(range, 'rest')
        xlim_2 = [12 35];
    end
else
    % do nothing
end

% calibration
if calib
    off1 = -0.0765;
    off2 = 0.149;
else
    off1 = 0;
    off2 = 0;
end
%% Import
% Copter
%copData = csvread('Coptersonde1_Data_2017-07-18_16h30m13s.csv',1,1);
if testNum == 2
    copData = csvread('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/Coptersonde1_Data_2017-09-21_14h09m53s.csv',1,1);
elseif testNum == 3
    copData = csvread('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/Coptersonde1_Data_2017-09-21_14h43m56s.csv',1,1);
end
tCop = copData(:,1)/(24*60*60) + datenum(1970,1,1,0,0,0);
T1 = copData(:,24) + off1;
T2 = copData(:,25) + off2;
RH1 = copData(:,14);
RH2 = copData(:,15);
RHT1 = copData(:,18);
RHT2 = copData(:,19);

% Raw copter - iMet, Current
if testNum == 2
    iMet = load('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/00000050.BIN-634998.mat', 'IMET');
    CURR = load('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/00000050.BIN-634998.mat', 'CURR');
elseif testNum == 3
    iMet = load('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/00000051.BIN-360329.mat', 'IMET');
    CURR = load('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/00000051.BIN-360329.mat', 'CURR');
end
T1raw = iMet.IMET(:, 7) + off1;
T2raw = iMet.IMET(:, 9) + off2;
curr = CURR.CURR(:, 4);
tImet = iMet.IMET(:, 2);

% Motor position
if testNum == 2
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/motor_Thu_Sep_21_14_20_33_2017.csv');
elseif testNum == 3
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/motor_Thu_Sep_21_14_46_10_2017.csv');
end
motData = textscan(fid, '%f %s %s','headerlines',1,'Delimiter',',');
fclose(fid);
motPos = motData{1,1}(:); % in
timeArv = datenum(motData{1,2},'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);
timeDep = datenum(motData{1,3},'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);

tMot = nan(2*length(motPos), 1);
xMot = nan(length(tMot), 1);
count = 1;
for i = 1:length(timeDep)
    tMot(count) = timeArv(i);
    xMot(count) = motPos(i);
    tMot(count + 1) = timeDep(i);
    xMot(count + 1) = motPos(i);
    count = count + 2;
end

% Wind speed
if testNum == 2
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/wind_Thu_Sep_21_14_20_33_2017.csv');
elseif testNum == 3
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/wind_Thu_Sep_21_14_46_10_2017.csv');
end
winds = textscan(fid, '%s %f','headerlines',1,'Delimiter',',');
fclose(fid);
tWind = datenum(winds{1,1}(1:2:end),'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);
speed = winds{1,2}(1:2:end);

% import both speeds for testNum == 3
if testNum == 3
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/wind_Thu_Sep_21_14_20_33_2017.csv');
    winds2 = textscan(fid, '%s %f','headerlines',1,'Delimiter',',');
    fclose(fid);
    speed2 = winds2{1,2}(1:2:end);
end

% MM
if testNum == 1
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/SensorPlacementTest.csv');
    mmText = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %f %f %f %f %f %f','headerlines',4,'Delimiter',',');
    fclose(fid);
    mmtime = datenum(mmText{1,1},'mm/dd/yyyy HH:MM');
    count = 0;
    for i = 1:length(mmtime)
        if count == 60
            count = 0;
        end
        mmtime(i) = mmtime(i) + datenum(0,0,0,0,0,count);
        count = count + 1;
    end

    mmTuav = mmText{1,8}(:);
    mmTenv = mmText{1,7}(:);
else
    fid = fopen('/Users/briangreene/Nextcloud/Projects/RILChamber/V2/NSSL data 9_21_2017_185745_obs.txt');
    mmText = textscan(fid, '%f %s %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %f %f %f %f %f', 'Delimiter', ',');
    fclose(fid);
    mmtime = datenum(mmText{1, 2}, 'yyyy-mm-dd HH:MM:SS.fff');
    mmTuav = mmText{1, 8}(:);
    mmTenv = mmText{1, 7}(:);
end

%% Filter bad data
% T1(T1<10) = nan;
% T2(T2<10) = nan;
% RH1(RH1<10) = nan;
% RH2(RH2<10) = nan;

% T1 = T1(1:4214);
% T2 = T2(1:4214);
% tCop = tCop(1:4214);

%% Temp running average
T1run = nan(length(T1),1);
T2run = nan(length(T2),1);
for i = 100:length(T1)
    T1run(i) = nanmean(T1(i-99:i));
    T2run(i) = nanmean(T2(i-99:i));
end

T1rawrun = nan(length(T1raw),1);
T2rawrun = nan(length(T2raw),1);
for i = 100:length(T1raw)
    T1rawrun(i) = nanmean(T1raw(i-99:i));
    T2rawrun(i) = nanmean(T2raw(i-99:i));
end

%% Align raw data with transmitted with corrcoef
% n = length(T2);
% m = length(T2raw);
% R = nan(m-n, 1);
% 
% for i = 1:length(R)
%     tempCorr = corrcoef(T2, T2raw(i:i+n-1));
%     R(i) = tempCorr(1,2);
% end
% [val, offset] = nanmax(R);
% 
% % create new time array for raw data
% baseT = tCop(1);
% deltaT = datenum(0,0,0,0,0,0.1);
% tCopRaw = nan(m, 1);
% 
% for i = 1:m
%     tCopRaw(i) = baseT - (offset - i) * deltaT;
% end

%% Align raw data by eye test
n = length(T2);
m = length(T2raw);
if testNum == 2
    argmax = 1846;
    argmaxraw = 2645;
elseif testNum == 3
    [~, argmax] = nanmax(T2);
    [~, argmaxraw] = nanmax(T2raw);
end 

baseT = tCop(argmax);
deltaTrel = (tImet - tImet(argmaxraw)) / (1000 * 86400); % fraction of day
tCopRaw = nan(m, 1);

for i = 1:m
    tCopRaw(i) = baseT + deltaTrel(i);
end
tCurr = linspace(tCopRaw(1), tCopRaw(end), length(curr));
imm = find(tCopRaw(1) < mmtime & mmtime < tCopRaw(end));
tMMraw = mmtime(imm);

%% Convert everything to relative times in min - mm booted up first
offset = tMMraw(1);
conversion = 1440; % 86400 s/day * 1 day/24 hr * 1 hr/60 min
tCopRel = (tCopRaw - offset) * conversion;
tMMrel = (tMMraw - offset) * conversion;
tCurrRel = (tCurr - offset) * conversion;
tWindRel = (tWind - offset) * conversion;

if testNum == 2
    tCurrOn = tCurrRel(7250);
    %tCurrOn = tCurrRel(7285);
    %tCurrOn1 = tCurrRel(4640);
    tCurrOn1 = tCurrRel(4750);
    tCurrOff = tCurrRel(18200);
    tCurrOff1 = tCurrRel(5783);
    if anomaly
        ydash = [-0.2, 1];
    else
        ydash = [22.5, 24.5];
    end
elseif testNum == 3
    tCurrOn = tCurrRel(1513);
    tCurrOff = tCurrRel(10226);
    if anomaly
        ydash = [-.2, 0.8];
    else
        ydash = [23, 24.5];
    end
end

%% Correct for linearly increasing background temp
% [~, iMMon] = min(abs(tMMrel - tCurrOn));
% [~, iMMoff] = min(abs(tMMrel - tCurrOff));
% 
% dTdt = (mmTenv(iMMoff) - mmTenv(iMMon)) / (tCurrOff - tCurrOn);
% mmTenvCorr = mmTenv(imm) - dTdt * tMMrel;
% mmTuavCorr = mmTuav(imm) - dTdt * tMMrel;
% T1rawrunCorr = T1rawrun - dTdt * tCopRel;
% T2rawrunCorr = T2rawrun - dTdt * tCopRel;

%% Calculate departures from background temp to plot as anomalies
% loop through all iMet temps, find nearest MM timestep, find difference
T1anom = nan(length(tCopRel), 1);
T2anom = nan(length(tCopRel), 1);
for i = 1:length(tCopRel)
    [~, iMMClose] = min(abs(tMMrel - tCopRel(i)));
    T1anom(i) = T1rawrun(i) - mmTenv(imm(iMMClose));
    T2anom(i) = T2rawrun(i) - mmTenv(imm(iMMClose));
end

mmTanom = mmTuav - mmTenv;

%% Average across bins to plot versus distance
% T1mean = nan(length(motPos),1);
% T2mean = nan(length(motPos),1);
% speedmean = nan(length(motPos),1);
% % mmTuavmean = nan(length(motPos),1);
% % mmTenvmean = nan(length(motPos),1);
% for j = 1:length(T1mean)-1
%     % indices for each source
%     iT = find(timeArv(j) < tCop & tCop < timeDep(j));
%     % take speed from between arrival points
%     iS = find(timeArv(j) < tWind & tWind < timeArv(j+1));
% %     iM = find(timeArv(j) < mmtime & mmtime < timeDep(j));
%     
%     % average for each position
%     T1mean(j) = nanmean(T1run(iT));
%     T2mean(j) = nanmean(T2run(iT));
%     speedmean(j) = nanmean(speed(iS));
% %     mmTuavmean(j) = nanmean(mmTuav(iM));
% %     mmTenvmean(j) = nanmean(mmTenv(iM));
% end

%% Checkpoints under props - distance
% offSen = -2.3; % all in cm
% %offSen = 0;
% 
% tip1start = (offSen + 11.2) .* ones(2,1);
% mount1start = (tip1start + 9.5);
% mount1end = (mount1start + 4.2);
% tip1end = (mount1end + 9.5);
% 
% tip2start = (tip1end + 3.6);
% mount2start = (tip2start + 9.5);
% mount2end = (mount2start + 4.2);
% tip2end = (mount2end + 9.5);
% 
% ydash = [22.0 27.0];

%% Checkpoints under props - time
% t_tip1start = timeArv((tip1start(1) - 0.1) < motPos & motPos < (tip1start(1) + 0.1)) .* ones(2,1);
% t_mount1start = timeArv((mount1start(1) - 0.12) < motPos & motPos < (mount1start(1) + 0.12)) .* ones(2,1);
% t_mount1end = timeArv((mount1end(1) - 0.12) < motPos & motPos < (mount1end(1) + 0.12)) .* ones(2,1);
% t_tip1end = timeArv((tip1end(1) - 0.1) < motPos & motPos < (tip1end(1) + 0.1)) .* ones(2,1);
% 
% t_tip2start = timeArv((tip2start(1) - 0.1) < motPos & motPos < (tip2start(1) + 0.1)) .* ones(2,1);
% t_mount2start = timeArv((mount2start(1) - 0.1) < motPos & motPos < (mount2start(1) + 0.1)) .* ones(2,1);
% t_mount2end = timeArv((mount2end(1) - 0.1) < motPos & motPos < (mount2end(1) + 0.1)) .* ones(2,1);
% t_tip2end = timeArv((tip2end(1) - 0.1) < motPos & motPos < (tip2end(1) + 0.1)) .* ones(2,1);

%% Create x-coordinate for each source
% x_tCop = linspace(0,63.5,length(tCop(tWind(1) < tCop & tCop < tWind(end))));
% x_tWind = linspace(0,63.5,length(tWind));
% x_tmm = linspace(0,63.5,length(mmtime(tWind(1) < mmtime & mmtime < tWind(end))));

%% Plot
ax = figure(1);
ax.Position = [500 500 1440 810];

left_color = [0 0 0];
%right_color = [239/255 118/255 4/255];
right_color = [73/255, 0/255, 146/255];
set(ax,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left

hold on
if anomaly
    h1 = plot(tCopRel, T1anom, 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'iMet 1 - Below rwUAS', 'LineWidth', 1.5);
    h2 = plot(tCopRel, T2anom, 'Color', [0, 0, 230/255], 'LineStyle', '-', 'DisplayName', 'iMet 2 - In Tube', 'LineWidth', 2.5);
    h3 = plot(tMMrel, mmTanom(imm), 'Color', [230/255, 0, 0], 'LineStyle', '-', 'DisplayName', 'NSSL 1 - In Tube', 'LineWidth', 2.5);
    lzero = line([0, ceil(tCopRel(end))], [0, 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2.5);
    bx = gca;
    bx.XAxis.MinorTickValues = min(bx.XAxis.TickValues):0.5:max(bx.XAxis.TickValues);
    bx.YAxis(1).MinorTickValues = min(ydash):0.05:max(ydash);
else
    h1 = plot(tCopRel, T1rawrun, 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'iMet 1 - Below rwUAS', 'LineWidth', 1.5);
    h2 = plot(tCopRel, T2rawrun, 'Color', [0, 0, 230/255], 'LineStyle', '-', 'DisplayName', 'iMet 2 - In Tube', 'LineWidth', 2.5);
    h3 = plot(tMMrel, mmTuav(imm), 'Color', [230/255, 0, 0], 'LineStyle', '-', 'DisplayName', 'NSSL 1 - In Tube', 'LineWidth', 2.5);
    h4 = plot(tMMrel, mmTenv(imm), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 3, 'DisplayName', 'NSSL 2 - Background');
    bx = gca;
    bx.XAxis.MinorTickValues = min(bx.XAxis.TickValues):0.5:max(bx.XAxis.TickValues);
    bx.YAxis(1).MinorTickValues = min(ydash):0.1:max(ydash);
end

l1 = line([tCurrOn, tCurrOn], ydash, 'Color', [24/225 145/225 36/225], 'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
l2 = line([tCurrOff, tCurrOff], ydash, 'Color', [206/225 26/225 26/225], 'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

if testNum == 2
    l3 = line([tCurrOn1, tCurrOn1], ydash, 'Color', [24/225 145/225 36/225], 'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
    l4 = line([tCurrOff1, tCurrOff1], ydash, 'Color', [206/225 26/225 26/225], 'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');
end

% h1 = plot(tCopRaw,T1rawrun(tWind(1) < tCopRaw & tCopRaw < tWind(end)), 'DisplayName', 'iMet Sensor 1');
% h2 = plot(tCopRaw,T2rawrun(tWind(1) < tCopRaw & tCopRaw < tWind(end)),'Color','g','LineStyle','-', 'DisplayName', 'iMet Sensor 2');
% h3 = plot(x_tmm,mmTuav(tWind(1) < mmtime & mmtime < tWind(end)),'Color','r','LineStyle','-', 'DisplayName', 'NSSL Tfast');
% h4 = plot(x_tmm,mmTenv(tWind(1) < mmtime & mmtime < tWind(end)),'Color','k','LineStyle','-','LineWidth',2, 'DisplayName', 'NSSL Background Temp');
% plot([tip1start mount1start mount1end tip1end tip2start mount2start mount2end tip2end],ydash,'k--')
hold off
if anomaly
    ylabel('Temperature (^oC) Relative to NSSL Background Temperature','FontSize',22)
else
    ylabel('Temperature (^oC)','FontSize',22)
end
grid on
grid minor
bx.GridAlpha = 0.75;
bx.MinorGridAlpha = 0.5;
set(bx, 'FontSize', 22);

if testNum == 2
    yyaxis right
    hold on
    h5 = plot(tWindRel, speed, 'Color', [73/255, 0/255, 146/255], 'DisplayName', 'Wind Speed', 'LineWidth', 2.5);
    hold off
    h5.Color(4) = 0.5;
    ylim([0 10])
    xlim(xlim_2)
    ylabel('Wind Speed (m s^{-1})', 'FontSize', 22)
elseif testNum == 3
    yyaxis right
    h6 = plot(tWindRel, speed2(1:length(tWindRel)), 'DisplayName', 'Winds from Test 1', 'Color', [73/255, 0/255, 146/255], 'LineWidth', 2, 'LineStyle', ':');
    h6.Color(4) = 0.5;
    ylim([0 10])
    ylabel('Wind Speed (m s^{-1})', 'FontSize', 22)
end

if testNum == 2 || testNum == 3
    ylim([0, 20])
end
xlabel('Time Elapsed (minutes)','FontSize',22)

titleStr = {'', 'Wind Probe in Tube', 'No Wind Probe'};
%title(sprintf('Coptersonde Sensor Placement Experiment, %s', titleStr{testNum}), 'FontSize', 32)

% legend
if anomaly
    if testNum == 2
        leg = legend([h1 h2 h3 h5 l1 l2]);
        leg.Location = 'northwest';
        leg.FontSize = 20;
    elseif testNum == 3
        leg = legend([h1 h2 h3 h6 l1 l2]);
        leg.Location = 'northwest';
        leg.FontSize = 20;
    end
else
    if testNum == 2
        leg = legend([h1 h2 h3 h4 h5 l1 l2]);
        leg.Location = 'northwest';
        leg.FontSize = 20;
    elseif testNum == 3
        leg = legend([h1 h2 h3 h4 h6 l1 l2]);
        leg.Location = 'northwest';
        leg.FontSize = 20;
    end
end