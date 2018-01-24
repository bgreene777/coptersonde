clc
clear
close all

%% Import
myNextcloud = '/Users/briangreene/Nextcloud/thermo/data/RILChamber/V2/';

fnameCop = 'Coptersonde1_Data_2018-01-23_18h27m47s.csv';
fpathCop = [myNextcloud fnameCop];
copData = csvread(fpathCop, 1, 1);

fnameMeso = '20180124.NWCM.1min.csv';
fpathMeso = [myNextcloud fnameMeso];
mesoData = csvread(fpathMeso, 1);

%% Assign
T = copData(:, 25);
t = copData(:, 1);
az = copData(:, 13);
ay = copData(:, 12);

T2m = mesoData(:, 3);
T9m = mesoData(:, 10);
tmeso = mesoData(:, 1);

%% Convert t
tBase = t(1);
conversion = 1/60;
tRel = (t - tBase) * conversion;

imeso = [28:58];
tmesoRel = tmeso(imeso) - tmeso(28);

%% Only use valid data
iend = 14942;

%% Manually determine rotor throttle up and down
is1 = 1190;
ie1 = 1736;
is2 = 2258;
ie2 = 2795;
is3 = 3210;
ie3 = 3732;
is4 = 4232;
ie4 = 4780;
is5 = 5820;
ie5 = 6371;
is6 = 7387;
ie6 = 7943;
is7 = 8937;
ie7 = 11472;
is8 = 12500;
ie8 = 14580;

tOn1 = tRel(is1);
tOff1 = tRel(ie1);
tOn2 = tRel(is2);
tOff2 = tRel(ie2);
tOn3 = tRel(is3);
tOff3 = tRel(ie3);
tOn4 = tRel(is4);
tOff4 = tRel(ie4);
tOn5 = tRel(is5);
tOff5 = tRel(ie5);
tOn6 = tRel(is6);
tOff6 = tRel(ie6);
tOn7 = tRel(is7);
tOff7 = tRel(ie7);
tOn8 = tRel(is8);
tOff8 = tRel(ie8);

ydash = [5.5 9];

%% Plot
ax = figure(1);
ax.Position = [500 500 1440 810];

h1 = plot(tRel(1:iend), T(1:iend), 'Color', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'iMet Temperature');
hold on
h2 = plot(tmesoRel, T2m(imeso), 'b*', 'LineWidth', 2, 'DisplayName',... 
    'NWC Mesonet 2m T');
h3 = plot(tmesoRel, T9m(imeso), 'r*', 'LineWidth', 2, 'DisplayName',... 
    'NWC Mesonet 9m T');

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

lOn4 = line([tOn4, tOn4], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff4 = line([tOff4, tOff4], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

lOn5 = line([tOn5, tOn5], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff5 = line([tOff5, tOff5], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

lOn6 = line([tOn6, tOn6], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff6 = line([tOff6, tOff6], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

lOn7 = line([tOn7, tOn7], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff7 = line([tOff7, tOff7], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');

lOn8 = line([tOn8, tOn8], ydash, 'Color', [24/225 145/225 36/225], ...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors On');
lOff8 = line([tOff8, tOff8], ydash, 'Color', [206/225 26/225 26/225],...
    'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Rotors Off');
hold off

leg = legend([h1 h2 h3 lOn1 lOff1]);
leg.Location = 'northeast';
leg.FontSize = 20;

xlabel('Time Elapsed (minutes)','FontSize',22)
ylabel('Temperature (^oC)', 'FontSize', 22)
title('Coptersonde Outdoor Aspiration Experiment', 'FontSize', 32)
grid on
grid minor
bx = gca;
bx.GridAlpha = 0.75;
bx.MinorGridAlpha = 0.5;
set(bx, 'FontSize', 22);