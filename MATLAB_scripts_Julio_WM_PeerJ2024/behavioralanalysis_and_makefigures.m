
% load human EEG data
load('erp_data_2to2_andstats.mat')


%% figure 1 - BEHAVIOR PER BLOCK

%for the behavior, the matrix i already sent you (behav_psychometric, size 137 x8) is split into blocks like this:
%TM1: rows 1:52
%TM2: rows 53:101
%TM3: rows 102:137


% Define logistic function
logistic_function = @(b, c, x) 1 ./ (1 + exp(-b * (x - c)));
% c is the sigmoid midpoint 
% b is the steepness of the curve

%sigmoid midpoint is not the same as P0.5 midpoint. need to plot to
%estimate. 

initial_guess = [1, 1];

%for stats 
for b = 1:3 %loop blocks
    for s = 1:24 %loop humans
        x = test_durations(b,:);
        y = plongb(s,:,b);
        %estimate psych curve and plot
        [fitted_model gof] = fit(x', y', logistic_function,'StartPoint', [1 1]);
        x_fit = linspace(min(x), max(x), 100);
        y_fit = feval(fitted_model, x_fit);
%estimate PSE
[val loc] = min(abs(y_fit-0.5));

%save point of subjective quality and slope
PSE_humans(s,b) = x_fit(loc);
slope_humans(s,b) = fitted_model.b;
    end

end

 for s = 1:52 %loop monkeys b1
 x = test_durations(1,:);
y = behav_dat(s,:);
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
%estimate PSE
[val loc] = min(abs(y_fit-0.5));
PSE_monkeys(s,1) =x_fit(loc);
slope_monkeys(s,1) = fitted_model.b;
 end

 for s = 53:101 %loop monkeys b2
x = test_durations(2,:);
y = behav_dat(s,:);
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
%estimate PSE
[val loc] = min(abs(y_fit-0.5));
PSE_monkeys(s,1) =x_fit(loc);
slope_monkeys(s,1) = fitted_model.b;
 end

 for s = 102:137 %loop monkeys b3
x = test_durations(3,:);
y = behav_dat(s,:);
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
%estimate PSE
[val loc] = min(abs(y_fit-0.5));
PSE_monkeys(s,1) =x_fit(loc);
slope_monkeys(s,1) = fitted_model.b;
 end

%exclude bad performers
index = nanmean(correct,2)<0.6;

%do stats in humans
%slope
[p,tbl_anova_humans_slope,stats] = anova1(slope_humans);
[h p ci stats] = ttest(slope_humans(:,1), slope_humans(:,2))
[h p ci stats] = ttest(slope_humans(:,2), slope_humans(:,3))
[h p ci stats] = ttest(slope_humans(:,1), slope_humans(:,3))
%CE
CE_humans = [PSE_humans(:,1)-0.35, PSE_humans(:,2)-0.685,PSE_humans(:,3)-1.195] ;
[p,tbl_anova_humans_CE,stats] = anova1(CE_humans);
[h p ci stats] = ttest(CE_humans(:,1), CE_humans(:,2))
[h p ci stats] = ttest(CE_humans(:,2), CE_humans(:,3))
[h p ci stats] = ttest(CE_humans(:,1), CE_humans(:,3))


%do stats in monkey
%slope
temp = nan(52,3);
temp(1:52,1) = slope_monkeys(1:52);
temp(1:49,2) = slope_monkeys(53:101);
temp(1:36,3) = slope_monkeys(102:end);
[p,tbl_anova_monkey_slope,stats] = anova1(temp)
[h p ci stats] = ttest2(temp(:,1), temp(:,2))
[h p ci stats] = ttest2(temp(:,2), temp(:,3))
[h p ci stats] = ttest2(temp(:,1), temp(:,3))

temp = nan(52,3);
temp(1:52,1) = PSE_monkeys(1:52);
temp(1:49,2) = PSE_monkeys(53:101);
temp(1:36,3) = PSE_monkeys(102:end);
%CE
CE_monkeys = [temp(:,1)-0.35, temp(:,2)-0.685,temp(:,3)-1.195] ;
[p,tbl_anova_humans_CE,stats] = anova1(CE_monkeys);
[h p ci stats] = ttest(CE_humans(:,1), CE_humans(:,2))
[h p ci stats] = ttest(CE_humans(:,2), CE_humans(:,3))
[h p ci stats] = ttest(CE_humans(:,1), CE_humans(:,3))

%compare monkeys and humans in PSE
[h p ci stats] = ttest2(PSE_humans(:,1),PSE_monkeys(1:52)) 
[h p ci stats] = ttest2(PSE_humans(:,2),PSE_monkeys(53:101)) 
[h p ci stats] = ttest2(PSE_humans(:,3),PSE_monkeys(102:end)) 


%compare monkeys and humans in slope
[h p ci stats] = ttest2(slope_humans(:,1),temp(1:52)) 
[h p ci stats] = ttest2(slope_humans(:,2),temp(53:101)) 
[h p ci stats] = ttest2(slope_humans(:,3),temp(102:end)) 


%for plotting
 for b = 1:3 %loop blocks

%human
x = test_durations(b,:);
y = mean(plongb(:,:,b),1);
err = std(plongb(:,:,b),1) / sqrt(24);
%plot data 
e = errorbar(x,y,err,'o', 'LineWidth', 2, 'Color', 'Black');
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = 'black';
e.CapSize = 10;
hold on 
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
plot(x_fit,y_fit, 'LineWidth', 2, 'Color', 'Black')


% monkey
if b ==1
y = mean(behav_dat(1:52,:),1);
err = std(behav_dat(1:52,:),1) / sqrt(52);
e = errorbar(x,y,err,'o', 'LineWidth', 2, 'Color', 'Green');
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = 'green';
e.CapSize = 10;
hold on
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
plot(x_fit,y_fit, 'LineWidth', 2, 'Color', 'green')


end

if b ==2
y = mean(behav_dat(53:101,:),1);
err = std(behav_dat(53:101,:),1) / sqrt(52);
e = errorbar(x,y,err,'o', 'LineWidth', 2, 'Color', 'Green');
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = 'green';
e.CapSize = 10;
hold on 
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
plot(x_fit,y_fit, 'LineWidth', 2, 'Color', 'green')
end

if b ==3
y = mean(behav_dat(102:137,:),1);
err = std(behav_dat(102:137,:),1) / sqrt(52);
e = errorbar(x,y,err,'o', 'LineWidth', 2, 'Color', 'Green');
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = 'green';
e.CapSize = 10;
hold on 
%estimate psych curve and plot
fitted_model = fit(x', y', logistic_function,'StartPoint', initial_guess);
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitted_model, x_fit);
plot(x_fit,y_fit, 'LineWidth', 2, 'Color', 'green')
end


 end


xlabel ('Stimulus duration');
ylabel ('P (long)');
legend('', 'Humans','', 'Monkey')
ax = gca;
ax.FontSize = 20;
box off


%% Figure 2 - HUMANS 2 PCAS
load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\erp_and_stats_pca2.mat')

bluecolor = [0.07,0.62,1.00];
orangecolor = [1.00,0.41,0.16];

%plot pca2 (already loaded)
subplot(2,4,5) 
topoplot(coeff(:,2),OUTEEG.chanlocs)
c = colorbar('southoutside', 'FontSize', 8);
c.Label.String = 'Weights (a.u.)';


subplot(2,4,6) 
a = stdshade(co_erp_long,0.1,bluecolor,EEG_epoched.times,0);
box off
hold on
a = stdshade(co_erp_short,0.1,orangecolor,EEG_epoched.times,0);
box off
hold on
sigtimes = logical([zeros(1,500), stat_co.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short', 'Location', 'northeast')
ylabel('Amplitude (µV)')
xlabel('time (s)')
xlim([-56 1996])
ax = gca;
ax.FontSize = 10;


subplot(2,4,7) 
a = stdshade(inco_erp_long,0.1,bluecolor,EEG_epoched.times,0);
box off
hold on
a = stdshade(inco_erp_short,0.1,orangecolor,EEG_epoched.times,0);
box off
hold on
sigtimes = logical([zeros(1,500), stat_inco.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short', 'Location', 'northeast')
ylabel('Amplitude (µV)')
xlabel('time(s)')
xlim([-56 1996])
ax = gca;
ax.FontSize = 10;

subplot(2,4,8)
a = stdshade(temp11,0.1,bluecolor,EEG_epoched.times,0);
hold on
a = stdshade(temp22,0.1,orangecolor,EEG_epoched.times,0);
hold on
sigtimes = logical([zeros(1,500), stat_co_sameduration.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short')
ylabel('Amplitude (µV)')
xlabel('time(s)')
ax = gca;
ax.FontSize = 10;
xlim([-56 1996])
box off

%load pca1
load('erp_and_stats_pca1.mat')
subplot(2,4,1) 
topoplot(coeff(:,1),OUTEEG.chanlocs)
c = colorbar('southoutside', 'FontSize', 8);
c.Label.String = 'Weights (a.u.)';
title('PCA-derived spatial filter', 'FontSize',12)

subplot(2,4,2) 
a = stdshade(co_erp_long,0.1,bluecolor,EEG_epoched.times,0);
box off
hold on
a = stdshade(co_erp_short,0.1,orangecolor,EEG_epoched.times,0);
box off
hold on
sigtimes = logical([zeros(1,500), stat_co.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short', 'Location', 'northeast')
ylabel('Amplitude (µV)')
xlabel('time (s)')
xlim([-56 1996])
ax = gca;
ax.FontSize = 10;
title('Correct all', 'FontSize',12)

subplot(2,4,3) 
a = stdshade(inco_erp_long,0.1,bluecolor,EEG_epoched.times,0);
box off
hold on
a = stdshade(inco_erp_short,0.1,orangecolor,EEG_epoched.times,0);
box off
hold on
sigtimes = logical([zeros(1,500), stat_inco.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short', 'Location', 'northeast')
ylabel('Amplitude (µV)')
xlabel('time(s)')
xlim([-56 1996])
ax = gca;
ax.FontSize = 10;
title('Incorrect all', 'FontSize',12)

subplot(2,4,4)
a = stdshade(temp11,0.1,bluecolor,EEG_epoched.times,0);
hold on
a = stdshade(temp22,0.1,orangecolor,EEG_epoched.times,0);
hold on
sigtimes = logical([zeros(1,500), stat_co_sameduration.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short')
ylabel('Amplitude (µV)')
xlabel('time(s)')
ax = gca;
ax.FontSize = 10;
title('Correct same duration', 'FontSize',12)
xlim([-56 1996])
box off


%% Figure 3 - HUMANS FRONTAL CLUSTER VS MONKEYS

bluecolor = [0.07,0.62,1.00];
orangecolor = [1.00,0.41,0.16];

%load results frontal cluster
load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\erp_and_stats_frontalcluster_fzf1f2fczfc1fc2.mat')
load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\erp_and_stats_frontalcluster_sameduration.mat')
subplot(3,4,1)
topoplot(ones(1,8),EEG_epoched.chanlocs(e),'electrodes', 'on')
colorbar

subplot(3,4,2) 
a = stdshade(co_erp_long,0.1,bluecolor,EEG_epoched.times,0);
box off
hold on
a = stdshade(co_erp_short,0.1,orangecolor,EEG_epoched.times,0);
box off
hold on
sigtimes = logical([zeros(1,500), stat_co.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short', 'Location', 'northeast')
ylabel('Amplitude (µV)')
xlabel('time (s)')
xlim([-56 1996])
ax = gca;
ax.FontSize = 10;
title('Correct all', 'FontSize',12)

subplot(3,4,3) 
a = stdshade(temp11,0.1,bluecolor,EEG_epoched.times,0);
hold on
a = stdshade(temp22,0.1,orangecolor,EEG_epoched.times,0);
hold on
sigtimes = logical([zeros(1,500), stat_co_sameduration.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
ylabel('Amplitude (µV)')
xlabel('time(s)')
ax = gca;
ax.FontSize = 10;
title('Correct same duration', 'FontSize',12)
xlim([-56 1996])
box off

subplot(3,4,4)
a = stdshade(inco_erp_long,0.1,bluecolor,EEG_epoched.times,0);
box off
hold on
a = stdshade(inco_erp_short,0.1,orangecolor,EEG_epoched.times,0);
box off
hold on
sigtimes = logical([zeros(1,500), stat_inco.mask]);
scatter(EEG_epoched.times(sigtimes),repmat(0,1,size(EEG_epoched.times(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
ylabel('Amplitude (µV)')
xlabel('time(s)')
xlim([-56 1996])
ax = gca;
ax.FontSize = 10;
title('Incorrect all', 'FontSize',12)


%%%%%%%%%%%%%%%%%%%%%% Monkeys LFP %%%%%%%%%%
load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\monkey_data\ERP_m1_preSMA.mat')
time = time*1000;

filename = 'C:\Users\bjjr014\Desktop\cattask paper\Figures\monkeymri.png';
y = imread(filename, 'BackgroundColor', [1 1 1]);
subplot(3,4,[5 9]) 
imshow(y);


subplot(3,4,[6]) 
a = stdshade(long_correct,0.1,bluecolor,time,0);
box off
hold on
a = stdshade(short_correct,0.1,orangecolor,time,0);
box off
hold on
sigtimes = logical([zeros(1,300), stat_co_monkey.mask]);
scatter(time(sigtimes),repmat(0,1,size(time(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
legend('', 'long', '', 'short', 'Location', 'northeast')
ylabel('Amplitude (µV)')
xlabel('time (s)')
xlim([-56 500])
ax = gca;
ax.FontSize = 10;


subplot(3,4,8) 
a = stdshade(long_incorrect,0.1,bluecolor,time,0);
box off
hold on
a = stdshade(short_incorrect,0.1,orangecolor,time,0);
box off
hold on
sigtimes = logical([zeros(1,300), stat_inco_monkey.mask]);
scatter(time(sigtimes),repmat(0,1,size(time(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
ylabel('Amplitude (µV)')
xlabel('time (s)')
xlim([-56 500])
ax = gca;
ax.FontSize = 10;


subplot(3,4,7) 
time_samed = (-1:1/500:1)*1000;
a = stdshade(longsamedur,0.1,bluecolor,time_samed,0);
hold on
a = stdshade(shortsamedur,0.1,orangecolor,time_samed,0);
hold on
sigtimes = logical([zeros(1,500), stat_co_monkey_samedur.mask]);
scatter(time_samed(sigtimes),repmat(0,1,size(time_samed(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')
ylabel('Amplitude (µV)')
xlabel('time(s)')
title('correct')
ax = gca;
ax.FontSize = 10;
xlim([-56 500])
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%% Monkeys firing rate

load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\monkey_data\firingrate_m1_preSMA.mat')
time = time*1000;
time2save = find(time==0):find(time==500);

firrat_long_correct = movmean(zscore(firrat_long_correct,0,2),meansamples,2);
firrat_long_incorrect = movmean(zscore(firrat_long_incorrect,0,2),meansamples,2);
firrat_short_correct = movmean(zscore(firrat_short_correct,0,2),meansamples,2);
firrat_short_incorrect = movmean(zscore(firrat_short_incorrect,0,2),meansamples,2);



subplot(3,4,10) 
a = stdshade(firrat_long_correct,0.1,bluecolor,time,0);
box off
hold on
a = stdshade(firrat_short_correct,0.1,orangecolor,time,0);
box off
hold on
sigtimes = logical([zeros(1,find(time==0)), stat_firrat_co_monkey.mask]);
scatter(time(sigtimes),repmat(0,1,size(time(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')

ylabel('Mean firing rate (zscore)')
xlabel('time(s)')
xlim([-56 500])
ax = gca;
ax.FontSize = 10;

subplot(3,4,12) 
a = stdshade(firrat_long_incorrect,0.1,bluecolor,time,0);
hold on
a = stdshade(firrat_short_incorrect,0.1,orangecolor,time,0);
hold on
sigtimes = logical([zeros(1,find(time==0)), stat_firrat_inco_monkey.mask]);
scatter(time(sigtimes),repmat(0,1,size(time(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')

ylabel('Mean firing rate (zscore)')
xlabel('time(s)')
xlim([-56 500])
box off
ax = gca;
ax.FontSize = 10;

subplot(3,4,11) 
load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\monkey_data\ERP_m1_preSMA_overlaps.mat')
time_samed = (-0.5:1/500:0.5)*1000;
a = stdshade(longsamedur_firrat,0.1,bluecolor,time_samed,0);
hold on
a = stdshade(shortsamedur_firrat,0.1,orangecolor,time_samed,0);
hold on
sigtimes = logical([zeros(1,251), stat_co_monkey_samedur_firrat.mask]);
scatter(time_samed(sigtimes),repmat(0,1,size(time_samed(sigtimes),2)),10, 'Marker', '*', 'MarkerEdgeColor', 'black')

ylabel('Mean firing rate (zscore)')
xlabel('time(s)')
ax = gca;
ax.FontSize = 10;

xlim([-56 500])
box off


