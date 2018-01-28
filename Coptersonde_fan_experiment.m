clc
clear
close all

%% Import
fpathSens = '/Users/briangreene/Desktop/Coptersonde1_Data_2018-01-26_18h44m29s.csv';
sensData = csvread(fpathSens, 1, 1);

load('/Users/briangreene/Desktop/00000017.BIN-595661.mat');

%% Assign
T1 = sensData(:, 22) - 273.15;
T3 = sensData(:, 24) - 273.15;
T4 = sensData(:, 25) - 273.15;
tRel = (sensData(:, 1) - sensData(1, 1)) / 60;

%% Offsets
isCalib = 0;
if isCalib
    T1 = T1 + 0.145;
    T3 = T3 + 0.141;
    T4 = T4 + 0.0288;
end

%% Switch on
ion1 = 3059;
ioff1 = 4790;
ion2 = 6277;
ioff2 = 7381;
ion3 = 9078;
ioff3 = 10231;

tOn1 = tRel(ion1);
tOff1 = tRel(ioff1);
tOn2 = tRel(ion2);
tOff2 = tRel(ioff2);
tOn3 = tRel(ion3);
tOff3 = tRel(ioff3);

%% Vertical lines
ydash = [20.5 24.5];

%% Plot
ax = figure(1);
ax.Position = [500 500 1440 810];

h1 = plot(tRel, T1, 'Color', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'iMet Temperature 1');
hold on
h2 = plot(tRel, T3, 'Color', 'r', 'LineWidth', 1.5, ...
    'DisplayName', 'iMet Temperature 2');
h3 = plot(tRel, T4, 'Color', 'b', 'LineWidth', 1.5, ...
    'DisplayName', 'iMet Temperature 3');

lOn1 = line([tOn1, tOn1], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff1 = line([tOff1, tOff1], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

lOn2 = line([tOn2, tOn2], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff2 = line([tOff2, tOff2], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

lOn3 = line([tOn3, tOn3], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff3 = line([tOff3, tOff3], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

hold off

leg = legend([h1 h2 h3 lOn1 lOff1]);
leg.Location = 'southwest';
leg.FontSize = 20;

xlabel('Time Elapsed (minutes)','FontSize',22)
ylabel('Temperature (^oC)', 'FontSize', 22)
title('Coptersonde Fan Aspiration Experiment', 'FontSize', 32)
grid on
grid minor
bx = gca;
bx.GridAlpha = 0.75;
bx.MinorGridAlpha = 0.5;
set(bx, 'FontSize', 22);
%xlim([0 8])