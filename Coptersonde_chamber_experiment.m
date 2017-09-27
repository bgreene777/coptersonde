clc
clear
close all

%% Import
% Copter
%copData = csvread('Coptersonde1_Data_2017-07-18_16h30m13s.csv',1,1);
copData = csvread('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/Coptersonde1_Data_2017-09-21_14h09m53s.csv',1,1);
tCop = copData(:,1)/(24*60*60) + datenum(1970,1,1,0,0,0);
T1 = copData(:,24);
T2 = copData(:,25);
RH1 = copData(:,14);
RH2 = copData(:,15);
RHT1 = copData(:,18);
RHT2 = copData(:,19);

% Raw copter - iMet
iMet = load('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/00000050.BIN-634998.mat', 'IMET');
T1raw = iMet.IMET(:, 7);
T2raw = iMet.IMET(:, 9);

% Motor position
fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/motor_Thu_Sep_21_14_20_33_2017.csv');
motData = textscan(fid, '%f %s %s','headerlines',1,'Delimiter',',');
fclose(fid);
motPos = motData{1,1}(:) * 2.54; % cm
timeArv = datenum(motData{1,2},'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);
timeDep = datenum(motData{1,3},'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);

% Wind speed
fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/wind_Thu_Sep_21_14_20_33_2017.csv');
winds = textscan(fid, '%s %f','headerlines',1,'Delimiter',',');
fclose(fid);
tWind = datenum(winds{1,1}(1:2:end),'ddd mmm dd HH:MM:SS yyyy') + datenum(0,0,0,5,0,0);
speed = winds{1,2}(:);

% % MM
% fid = fopen('/Users/briangreene/Nextcloud/thermo/data/RILChamber/SensorPlacementTest.csv');
% mmText = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %f %f %f %f %f %f','headerlines',4,'Delimiter',',');
% fclose(fid);
% mmtime = datenum(mmText{1,1},'mm/dd/yyyy HH:MM');
% count = 0;
% for i = 1:length(mmtime)
%     if count == 60
%         count = 0;
%     end
%     mmtime(i) = mmtime(i) + datenum(0,0,0,0,0,count);
%     count = count + 1;
% end
% 
% mmTuav = mmText{1,8}(:);
% mmTenv = mmText{1,7}(:);

%% Filter bad data
T1(T1<10) = nan;
T2(T2<10) = nan;
RH1(RH1<10) = nan;
RH2(RH2<10) = nan;

%% Temp running average
T1run = nan(length(T1),1);
T2run = nan(length(T2),1);
for i = 100:length(T1)
    T1run(i) = nanmean(T1(i-99:i));
    T2run(i) = nanmean(T2(i-99:i));
end

T1rawrun = nan(length(T1raw),1);
T2rawrun = nan(length(T2raw),1);
%% Average across bins to plot versus distance
T1mean = nan(length(motPos),1);
T2mean = nan(length(motPos),1);
speedmean = nan(length(motPos),1);
mmTuavmean = nan(length(motPos),1);
mmTenvmean = nan(length(motPos),1);
for j = 1:length(T1mean)-1
    % indices for each source
    iT = find(timeArv(j) < tCop & tCop < timeDep(j));
    % take speed from between arrival points
    iS = find(timeArv(j) < tWind & tWind < timeArv(j+1));
    iM = find(timeArv(j) < mmtime & mmtime < timeDep(j));
    
    % average for each position
    T1mean(j) = nanmean(T1run(iT));
    T2mean(j) = nanmean(T2run(iT));
    speedmean(j) = nanmean(speed(iS));
    mmTuavmean(j) = nanmean(mmTuav(iM));
    mmTenvmean(j) = nanmean(mmTenv(iM));
end

%% Checkpoints under props - distance
offSen = -2.3; % all in cm
%offSen = 0;

tip1start = (offSen + 11.2) .* ones(2,1);
mount1start = (tip1start + 9.5);
mount1end = (mount1start + 4.2);
tip1end = (mount1end + 9.5);

tip2start = (tip1end + 3.6);
mount2start = (tip2start + 9.5);
mount2end = (mount2start + 4.2);
tip2end = (mount2end + 9.5);

ydash = [22.0 27.0];

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
x_tCop = linspace(0,63.5,length(tCop(tWind(1) < tCop & tCop < tWind(end))));
x_tWind = linspace(0,63.5,length(tWind));
x_tmm = linspace(0,63.5,length(mmtime(tWind(1) < mmtime & mmtime < tWind(end))));

%% Plot
ax = figure(1);
ax.Position = [500 500 1400 700];
yyaxis left
hold on
h1 = plot(x_tCop,T1run(tWind(1) < tCop & tCop < tWind(end)), 'DisplayName', 'iMet Sensor 1');
h2 = plot(x_tCop,T2run(tWind(1) < tCop & tCop < tWind(end)),'Color','g','LineStyle','-', 'DisplayName', 'iMet Sensor 2');
h3 = plot(x_tmm,mmTuav(tWind(1) < mmtime & mmtime < tWind(end)),'Color','r','LineStyle','-', 'DisplayName', 'NSSL Tfast');
h4 = plot(x_tmm,mmTenv(tWind(1) < mmtime & mmtime < tWind(end)),'Color','k','LineStyle','-','LineWidth',2, 'DisplayName', 'NSSL Background Temp');
plot([tip1start mount1start mount1end tip1end tip2start mount2start mount2end tip2end],ydash,'k--')
hold off
ylabel('Temperature(^oC)','FontSize',14)
% annotate
annotation('textarrow',[0.195714285714286 0.237142857142857],...
    [0.749 0.704285714285714],'String',{'Begin Prop 1'});
annotation('textarrow',[0.337857142857143 0.351428571428571],...
    [0.161428571428571 0.181428571428571],'String',{'Begin Motor Mount 1'});
annotation('textarrow',[0.417142857142857 0.407142857142857],...
    [0.162857142857143 0.184285714285714],'String',{'End Motor Mount 1'});
annotation('textarrow',[0.48 0.520714285714286],...
    [0.873285714285714 0.844285714285714],'String',{'End Tip 1'});
annotation('textarrow',[0.606428571428571 0.567857142857143],...
    [0.873285714285714 0.844285714285714],'String',{'Begin Tip 2'});
annotation('textarrow',[0.663571428571429 0.681428571428571],...
    [0.165714285714286 0.194285714285714],'String',{'Begin Motor Mount 2'});
annotation('textarrow',[0.755 0.735],...
    [0.164285714285714 0.194285714285714],'String',{'End Motor Mount 2'});
annotation('textarrow',[0.82 0.847857142857143],...
    [0.866142857142857 0.845714285714286],'String',{'End Tip 2'});

yyaxis right
h5 = plot(x_tWind,speed, 'DisplayName', 'Wind Speed');
ylim([0 10])
ylabel('Wind Speed (m/s)','FontSize',14)
%datetick('x','HH:MM:SS')

h = get(gca,'Children');
xlim([min(h.XData), max(h.XData)])
ylim([0, 22])
xlabel('Distance moved along linear actuator arm (cm)','FontSize',14)


title('Coptersonde Sensor Placement Experiment','FontSize',16)
grid on
legend([h1 h2 h3 h4 h5],'Location', 'northwest')
