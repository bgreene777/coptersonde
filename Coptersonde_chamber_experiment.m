clc
clear
close all

%% User parameters
testNum = 2;

%% Import
% Copter
%copData = csvread('Coptersonde1_Data_2017-07-18_16h30m13s.csv',1,1);
if testNum == 2
    copData = csvread('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/Coptersonde1_Data_2017-09-21_14h09m53s.csv',1,1);
elseif testNum == 3
    copData = csvread('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/Coptersonde1_Data_2017-09-21_14h43m56s.csv',1,1);
end
tCop = copData(:,1)/(24*60*60) + datenum(1970,1,1,0,0,0);
T1 = copData(:,24);
T2 = copData(:,25);
RH1 = copData(:,14);
RH2 = copData(:,15);
RHT1 = copData(:,18);
RHT2 = copData(:,19);

% Raw copter - iMet, Current
if testNum == 2
    iMet = load('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/00000050.BIN-634998.mat', 'IMET');
    CURR = load('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/00000050.BIN-634998.mat', 'CURR');
elseif testNum == 3
    iMet = load('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/00000051.BIN-360329.mat', 'IMET');
    CURR = load('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/00000051.BIN-360329.mat', 'CURR');
end
T1raw = iMet.IMET(:, 7);
T2raw = iMet.IMET(:, 9);
curr = CURR.CURR(:, 4);

% Motor position
if testNum == 2
    fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/motor_Thu_Sep_21_14_20_33_2017.csv');
elseif testNum == 3
    fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/motor_Thu_Sep_21_14_46_10_2017.csv');
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
    fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/wind_Thu_Sep_21_14_20_33_2017.csv');
elseif testNum == 3
    fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/wind_Thu_Sep_21_14_46_10_2017.csv');
end
winds = textscan(fid, '%s %f','headerlines',1,'Delimiter',',');
fclose(fid);
tWind = datenum(winds{1,1}(1:2:end),'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);
speed = winds{1,2}(1:2:end);

% MM
if testNum == 1
    fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/SensorPlacementTest.csv');
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
    fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/NSSL data 9_21_2017_185745_obs.txt');
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
    [min, argmax] = nanmax(T2);
    [min2, argmaxraw] = nanmax(T2raw);
end
    

baseT = tCop(argmax);
deltaT = datenum(0,0,0,0,0,0.15);
tCopRaw = nan(m, 1);

for i = 1:m
    tCopRaw(i) = baseT - (argmaxraw - i) * deltaT;
end
tCurr = linspace(tCopRaw(1), tCopRaw(end), length(curr));
imm = find(tCopRaw(1) < mmtime & mmtime < tCopRaw(end));

if testNum == 2
    tCurrOn = tCurr(7250);
    tCurrOn1 = tCurr(4640);
    tCurrOff = tCurr(18200);
    tCurrOff1 = tCurr(5783);
    ydash = [22.5, 24.5];
elseif testNum == 3
    tCurrOn = tCurr(1526);
    tCurrOff = tCurr(10226);
    ydash = [23, 24.5];
end

% hold on
% plot(tCopRaw, T2raw)
% plot(tCop, T2)
% hold off
% 
% figure
% hold on
% plot(tCopRaw, T2rawrun)
% plot(tCop, T2run)
% hold off

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
ax.Position = [500 500 1400 700];
yyaxis left
hold on
h1 = plot(tCopRaw, T1rawrun, 'DisplayName', 'iMet Sensor 1');
h2 = plot(tCopRaw, T2rawrun, 'Color', 'g', 'LineStyle', '-', 'DisplayName', 'iMet Sensor 2');
h3 = plot(mmtime(imm), mmTuav(imm), 'Color', 'r', 'LineStyle', '-', 'DisplayName', 'NSSL Tfast');
h4 = plot(mmtime(imm), mmTenv(imm), 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'NSSL Background Temperature');
line([tCurrOn, tCurrOn], ydash, 'Color', 'k', 'LineStyle', '--');
line([tCurrOff, tCurrOff], ydash, 'Color', 'k', 'LineStyle', '--');
if testNum == 2
    line([tCurrOn1, tCurrOn1], ydash, 'Color', 'k', 'LineStyle', '--');
    line([tCurrOff1, tCurrOff1], ydash, 'Color', 'k', 'LineStyle', '--');
end

% h1 = plot(tCopRaw,T1rawrun(tWind(1) < tCopRaw & tCopRaw < tWind(end)), 'DisplayName', 'iMet Sensor 1');
% h2 = plot(tCopRaw,T2rawrun(tWind(1) < tCopRaw & tCopRaw < tWind(end)),'Color','g','LineStyle','-', 'DisplayName', 'iMet Sensor 2');
% h3 = plot(x_tmm,mmTuav(tWind(1) < mmtime & mmtime < tWind(end)),'Color','r','LineStyle','-', 'DisplayName', 'NSSL Tfast');
% h4 = plot(x_tmm,mmTenv(tWind(1) < mmtime & mmtime < tWind(end)),'Color','k','LineStyle','-','LineWidth',2, 'DisplayName', 'NSSL Background Temp');
% plot([tip1start mount1start mount1end tip1end tip2start mount2start mount2end tip2end],ydash,'k--')
hold off
ylabel('Temperature(^oC)','FontSize',14)
% annotate
% annotation('textarrow',[0.195714285714286 0.237142857142857],...
%     [0.749 0.704285714285714],'String',{'Begin Prop 1'});
% annotation('textarrow',[0.337857142857143 0.351428571428571],...
%     [0.161428571428571 0.181428571428571],'String',{'Begin Motor Mount 1'});
% annotation('textarrow',[0.417142857142857 0.407142857142857],...
%     [0.162857142857143 0.184285714285714],'String',{'End Motor Mount 1'});
% annotation('textarrow',[0.48 0.520714285714286],...
%     [0.873285714285714 0.844285714285714],'String',{'End Tip 1'});
% annotation('textarrow',[0.606428571428571 0.567857142857143],...
%     [0.873285714285714 0.844285714285714],'String',{'Begin Tip 2'});
% annotation('textarrow',[0.663571428571429 0.681428571428571],...
%     [0.165714285714286 0.194285714285714],'String',{'Begin Motor Mount 2'});
% annotation('textarrow',[0.755 0.735],...
%     [0.164285714285714 0.194285714285714],'String',{'End Motor Mount 2'});
% annotation('textarrow',[0.82 0.847857142857143],...
%     [0.866142857142857 0.845714285714286],'String',{'End Tip 2'});

if testNum == 2
    yyaxis right
    hold on
    h5 = plot(tWind, speed, 'DisplayName', 'Wind Speed');
    h6 = plot(tMot, xMot, 'DisplayName', 'Position', 'Color', 'k', 'LineStyle', '-');
    hold off
    ylim([0 10])
    ylabel('Wind Speed (m/s) & Position (in)', 'FontSize', 14)
elseif testNum == 3
    yyaxis right
    hold on
    h6 = plot(tMot, xMot, 'DisplayName', 'Position', 'Color', 'k', 'LineStyle', '-');
    hold off
    ylim([0 10])
    ylabel('Position (in)', 'FontSize', 14)
end
    
ylim([0, 28])    
h = get(gca,'Children');
%xlim([min(h.XData), max(h.XData)])
xlabel('Time UTC','FontSize',14)
datetick('x','HH:MM:SS')


title(sprintf('Coptersonde Sensor Placement Experiment %d', testNum), 'FontSize',16)
grid on

if testNum == 2
    legend([h1 h2 h3 h4 h5 h6],'Location', 'northwest')
elseif testNum == 3
    legend([h1 h2 h3 h4 h6], 'Location', 'northwest')
end