%% load burst variables

load('detectburstfinal_2to40.mat')

%% average beta frequency range and estimate mean frequency

index = find(frequencies==15):find(frequencies==40); %index of frequencies to average

%divide mean duration by duration of 1cycle per frequency to get number of cycles 
for s = 1:size(bdur_1stdelay,1)
    for e = 1:size(bdur_1stdelay,2)
        for t = 1:size(bdur_1stdelay,3)
            bdur_1stdelay(s,e,t,:) =   squeeze(bdur_1stdelay(s,e,t,:)) ./  (1./frequencies');
            bdur_2nddelay(s,e,t,:) =   squeeze(bdur_2nddelay(s,e,t,:)) ./  (1./frequencies');
            bdur_baseline(s,e,t,:) =   squeeze(bdur_baseline(s,e,t,:)) ./  (1./frequencies');
        end
    end
end


% averaging frequency dimension
mean_bamp_delay = nanmean(bamp_1stdelay(:,:,:,index),4);
mean_bamp_delay2 = nanmean(bamp_2nddelay(:,:,:,index),4);
mean_bamp_fix = nanmean(bamp_baseline(:,:,:,index),4);

mean_bcov_delay = nanmean(bcov_1stdelay(:,:,:,index),4);
mean_bcov_delay2 = nanmean(bcov_2nddelay(:,:,:,index),4);
mean_bcov_fix = nanmean(bcov_baseline(:,:,:,index),4);

mean_bdur_delay = nanmean(bdur_1stdelay(:,:,:,index),4);
mean_bdur_delay2 = nanmean(bdur_2nddelay(:,:,:,index),4);
mean_bdur_fix = nanmean(bdur_baseline(:,:,:,index),4);

mean_brate_delay = nanmean(bnumb_1stdelay(:,:,:,index),4);
mean_brate_delay2 = nanmean(bnumb_2nddelay(:,:,:,index),4);
mean_brate_fix = nanmean(bnumb_baseline(:,:,:,index),4);

mean_bwidth_delay = nanmean(bwidth_1stdelay(:,:,:,index),4);
mean_bwidth_delay2 = nanmean(bwidth_2nddelay(:,:,:,index),4);
mean_bwidth_fix = nanmean(bwidth_baseline(:,:,:,index),4);


% getting mean frequency of bursts
%multiply number of bursts by their frequency and divide by number of
%bursts
for s = 1:size(bnumb_1stdelay,1)
    for e = 1:size(bnumb_1stdelay,2)
        for t = 1:size(bnumb_1stdelay,3)
            
            %get number of burst per frequency
            temp1  = squeeze(bnumb_2nddelay (s,e,t,index))';
            %estimate mean frequency
            mean_bfreq_delay2 (s,e,t) = sum(temp1.*frequencies(index)) / sum(temp1);
            

            %get number of burst per frequency
            temp1  = squeeze(bnumb_1stdelay (s,e,t,index))';
            %estimate mean frequency
            mean_bfreq_delay (s,e,t) = sum(temp1.*frequencies(index)) / sum(temp1);
            
            
             %get number of burst per frequency
            temp1  = squeeze(bnumb_baseline (s,e,t,index))';
            %estimate mean frequency
            mean_bfreq_fix (s,e,t) = sum(temp1.*frequencies(index)) / sum(temp1);
            
            
        end
    end
    
end


%% load trial information per subject
cd '/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/preprocessed_EEG/';
files = dir('*.mat');
clearvars cload accuracy instruction reactiontime
for f = 1:length(files)
     load(files(f).name, 'OUTEEG')
cload(f,:) = [OUTEEG.load(OUTEEG.load==1)';OUTEEG.load(OUTEEG.load==3)'];
accuracy(f,:) = [OUTEEG.accuracy(OUTEEG.load==1);OUTEEG.accuracy(OUTEEG.load==3)] / 100;
instruction(f,:) = [OUTEEG.condition(OUTEEG.load==1)';OUTEEG.condition(OUTEEG.load==3)'];
rt(f,:) = [OUTEEG.rt(OUTEEG.load==1);OUTEEG.rt(OUTEEG.load==3)];

end


%% plot descriptives Figure 2B
% restoredefaultpath
 addpath ('/Users/juliorodriguezlarios/Documents/MATLAB/stdshade.m')
 addpath ('/Users/juliorodriguezlarios/Documents/MATLAB/eeglab2022.0')
 eeglab;close

bamp = nanmean(cat(4,mean_bamp_delay,mean_bamp_delay2,mean_bamp_fix),4);
bfreq = nanmean(cat(4,mean_bfreq_delay,mean_bfreq_delay2,mean_bfreq_fix),4);
brate =  nanmean(cat(4,mean_brate_delay,mean_brate_delay2,mean_brate_fix),4);
bdur = nanmean(cat(4,mean_bdur_delay,mean_bdur_delay2,mean_bdur_fix),4);

%exclude outliers for descriptives
outliers = isoutlier(bamp);
bamp(outliers) = nan;
outliers = isoutlier(bfreq);
bfreq(outliers) = nan;
outliers = isoutlier(brate);
brate(outliers) = nan;
outliers = isoutlier(bdur );
bdur (outliers) = nan;


%plot
subplot(1,4,1)
temp = nanmean(nanmean(bfreq,3),1);
topoplot(temp,OUTEEG.chanlocs)
colorbar 
title('Beta burst frequency')
ax = gca;
ax.FontSize = 15;
ax.CLim = [20 25];

subplot(1,4,2)
temp = nanmean(nanmean(bamp,3),1);
topoplot(temp,OUTEEG.chanlocs)
colorbar 
title('Beta burst amplitude')
ax = gca;
ax.FontSize = 15;
ax.CLim(1) = 0;

subplot(1,4,3)
temp = nanmean(nanmean(brate,3),1);
topoplot(temp,OUTEEG.chanlocs)
colorbar 
title('Beta burst rate')
ax = gca;
ax.FontSize = 15;
ax.CLim(1) = 0;

subplot(1,4,4)
temp = nanmean(nanmean(bdur,3),1);
topoplot(temp,OUTEEG.chanlocs)
colorbar 
title('Beta burst duration')
ax = gca;
ax.FontSize = 15;
ax.CLim(1) = 2;

colormap parula


%% PERMUTATION TEST FOR CONDITION COMPARISON
% restoredefaultpath
% addpath ('/Users/juliorodriguezlarios/Documents/MATLAB/eeglab2022.0')
% eeglab;close;
% addpath('/Users/juliorodriguezlarios/Documents/MATLAB/fieldtrip-master')


%put group data in fieldtrip structure

%load fieldtrip structure
load('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/variables/ft_structure_fft.mat')
load('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/variables/foof_per_trial_laplacianfixed_justspectrum_nogaussian.mat', 'data1')

%create empty fieldtrip structures
freq_delay = freq_baseline;
freq_delay.label = data1.label;
freq_baseline.label = data1.label;
freq_delay.elec = data1.elec;
freq_baseline.elec = data1.elec;
freq_delay.freq = 1;
freq_baseline.freq = 1;
freq_delay1 = freq_delay;
freq_baseline1= freq_baseline;
freq_delay2 = freq_delay;
freq_baseline2= freq_baseline;
freq_delay4 = freq_delay;
freq_baseline4= freq_baseline;
freq_delay5 = freq_delay;
freq_baseline5= freq_baseline;


%%% DEFINE CONDITIONS %%%%
%%%%%%%%%%  delay vs fixation %%%%%%%%%%

freq_delay1.powspctrm = nanmean(mean_bamp_delay,3);
freq_baseline1.powspctrm = nanmean(mean_bamp_fix,3);

freq_delay2.powspctrm = nanmean(mean_bdur_delay,3);
freq_baseline2.powspctrm = nanmean(mean_bdur_fix,3);

freq_delay4.powspctrm = nanmean(mean_bfreq_delay,3);
freq_baseline4.powspctrm = nanmean(mean_bfreq_fix,3);

freq_delay5.powspctrm = nanmean(mean_brate_delay,3);
freq_baseline5.powspctrm = nanmean(mean_brate_fix,3);



%%%%%%%%%%%%%% load %%%%%%%%%%%%%%

% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_bamp_delay(s,:,cload(s,:)==3),3);
% temp3(s,:) = nanmean(mean_bamp_delay(s,:,cload(s,:)==1),3);
% end
% freq_delay1.powspctrm =temp2;
% freq_baseline1.powspctrm = temp3;
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_bdur_delay(s,:,cload(s,:)==3),3);
% temp3(s,:) = nanmean(mean_bdur_delay(s,:,cload(s,:)==1),3);
% end
% freq_delay2.powspctrm =temp2;
% freq_baseline2.powspctrm = temp3;
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_bfreq_delay(s,:,cload(s,:)==3),3);
% temp3(s,:) = nanmean(mean_bfreq_delay(s,:,cload(s,:)==1),3);
% end
% freq_delay4.powspctrm =temp2;
% freq_baseline4.powspctrm = temp3;
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_brate_delay(s,:,cload(s,:)==3),3);
% temp3(s,:) = nanmean(mean_brate_delay(s,:,cload(s,:)==1),3);
% end
% freq_delay5.powspctrm =temp2;
% freq_baseline5.powspctrm = temp3;



%%%%%%%%%%%%%% instruction %%%%%%%%%%%%%%
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_bamp_delay2(s,:,instruction(s,:)==2),3);
% temp3(s,:) = nanmean(mean_bamp_delay2(s,:,instruction(s,:)==1),3);
% end
% freq_delay1.powspctrm =temp2;
% freq_baseline1.powspctrm = temp3;
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_bdur_delay2(s,:,instruction(s,:)==2),3);
% temp3(s,:) = nanmean(mean_bdur_delay2(s,:,instruction(s,:)==1),3);
% end
% freq_delay2.powspctrm =temp2;
% freq_baseline2.powspctrm = temp3;
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_bfreq_delay2(s,:,instruction(s,:)==2),3);
% temp3(s,:) = nanmean(mean_bfreq_delay2(s,:,instruction(s,:)==1),3);
% end
% freq_delay4.powspctrm =temp2;
% freq_baseline4.powspctrm = temp3;
% 
% clearvars temp2 temp3
% for s = 1:length(files)
% temp2(s,:) = nanmean(mean_brate_delay2(s,:,instruction(s,:)==2),3);
% temp3(s,:) = nanmean(mean_brate_delay2(s,:,instruction(s,:)==1),3);
% end
% freq_delay5.powspctrm =temp2;
% freq_baseline5.powspctrm = temp3;



%permutation test
%for some reason it only works with this fieldtrip version
%restoredefaultpath
%addpath('/Users/juliorodriguezlarios/Documents/MATLAB/fieldtrip-master')
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 0;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freq_delay);
subj = size(freq_delay1.powspctrm,1);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;


%%%% amplitude
[stat1] = ft_freqstatistics(cfg, freq_delay1, freq_baseline1);

%%% duration
[stat2] = ft_freqstatistics(cfg, freq_delay2, freq_baseline2);

%%%frequency
[stat4] = ft_freqstatistics(cfg, freq_delay4, freq_baseline4);

%%%rate
[stat5] = ft_freqstatistics(cfg, freq_delay5, freq_baseline5);


%% PLOT CONDITION EFFECT 
f = figure;
f.Position(3) = 1000;
f.Position(4) = 500;
set(gcf,'color','w');

%define labels
condition1 = 'Delay';
condition2 = 'Fixation';
cstring = 't-value';


% get spectrum to plot

%delay vs fix
spectrum1 = squeeze(nanmean(zscore(allbursts_spectrum_1stdelay,0,4),3));
spectrum2 = squeeze(nanmean(zscore(allbursts_spectrum_baseline,0,4),3));

%load effect
% temp = zscore(allbursts_spectrum_1stdelay,0,4);
% for s = 1:size(allbursts_spectrum_1stdelay,1)
% spectrum1(s,:,:) = squeeze(nanmean(temp(s,:,cload(s,:)==3,:),3));
% spectrum2(s,:,:) = squeeze(nanmean(temp(s,:,cload(s,:)==1,:),3));
% end

% %instruction effect
% temp = zscore(allbursts_spectrum_2nddelay,0,4);
% for s = 1:size(allbursts_spectrum_2nddelay,1)
% spectrum1(s,:,:) = squeeze(nanmean(temp(s,:,instruction(s,:)==2,:),3));
% spectrum2(s,:,:) = squeeze(nanmean(temp(s,:,instruction(s,:)==1,:),3));
% end


%%%% amplitude
stat = stat1;
sig_m = find(stat.mask);
subplot(2,6,1)
topoplot(stat.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-10 10]);
c = colorbar('SouthOutside');
c.Label.String = cstring;
c.Position(3) = 0.05;
c.FontSize = 13;
if isempty(find(stat.mask))
stat.mask = ones(1,96);
end
temp1 = nanmean(freq_delay1.powspctrm(:,stat.mask),2);
temp2 = nanmean(freq_baseline1.powspctrm(:,stat.mask),2);
subplot(2,6,2)
coordLineStyle = 'k.';
boxplot([temp1, temp2],'Labels',{condition1,condition2}, 'Symbol', coordLineStyle); 
axes = get(gca,'Children');
ylabel('Amplitude (a.u.)')
for s = 3:length(axes.Children)
axes.Children(s).LineWidth = 2;
axes.Children(s).LineStyle = '-';
end
hold on;
parallelcoords([temp1, temp2], 'Color', [0.7 0.7 0.7], 'LineStyle', '-','LineWidth',0.2,...
  'Marker', '.', 'MarkerSize', 10);
 title('Beta burst amplitude')
 axes = gca;
axes.FontSize = 13;


%%% duration
stat = stat2;
sig_m = find(stat.mask);
subplot(2,6,3)
topoplot(stat.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
c = colorbar('SouthOutside');
c.Label.String = cstring;
c.FontSize = 13;
c.Position(3) = 0.05;
if isempty(find(stat.mask))
stat.mask = ones(1,96);
end
temp1 = nanmean(freq_delay2.powspctrm(:,stat.mask),2);
temp2 = nanmean(freq_baseline2.powspctrm(:,stat.mask),2);
subplot(2,6,4)
coordLineStyle = 'k.';
boxplot([temp1, temp2],'Labels',{condition1,condition2}, 'Symbol', coordLineStyle); 
ylabel('# cycles')
axes = get(gca,'Children')
for s = 3:length(axes.Children)
axes.Children(s).LineWidth = 2;
axes.Children(s).LineStyle = '-';
end
hold on;
parallelcoords([temp1, temp2], 'Color', [0.7 0.7 0.7], 'LineStyle', '-','LineWidth',0.2,...
  'Marker', '.', 'MarkerSize', 10);
 title('Beta burst duration')
 axes = gca;
axes.FontSize = 13;
 
%%%frequency
stat = stat4;
sig_m = find(stat.mask);
subplot(2,6,7)
topoplot(stat.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
c = colorbar('SouthOutside');
c.Label.String = cstring;
c.FontSize = 13;
c.Position(3) = 0.05;
if isempty(find(stat.mask))
stat.mask = ones(1,96);
end
temp1 = nanmean(freq_delay4.powspctrm(:,stat.mask),2);
temp2 = nanmean(freq_baseline4.powspctrm(:,stat.mask),2);
subplot(2,6,8)
coordLineStyle = 'k.';
boxplot([temp1, temp2],'Labels',{condition1,condition2}, 'Symbol', coordLineStyle); 
ylabel('Peak frequency (Hz)')
axes = get(gca,'Children')
for s = 3:length(axes.Children)
axes.Children(s).LineWidth = 2;
axes.Children(s).LineStyle = '-';
end
hold on;
parallelcoords([temp1, temp2], 'Color', [0.7 0.7 0.7], 'LineStyle', '-','LineWidth',0.2,...
  'Marker', '.', 'MarkerSize', 10);
 title('Beta burst frequency')
 axes = gca;
axes.FontSize = 13;

%%%rate
stat = stat5;
sig_m = find(stat.mask);
subplot(2,6,9)
topoplot(stat.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
c = colorbar('SouthOutside');
c.Label.String = cstring;
c.FontSize = 13;
c.Position(3) = 0.05;
if isempty(find(stat.mask))
stat.mask = ones(1,96);
end
temp1 = nanmean(freq_delay5.powspctrm(:,stat.mask),2);
temp2 = nanmean(freq_baseline5.powspctrm(:,stat.mask),2);
subplot(2,6,10)
coordLineStyle = 'k.';
boxplot([temp1, temp2],'Labels',{condition1,condition2}, 'Symbol', coordLineStyle); 
ylabel('# bursts')
axes = get(gca,'Children')
for s = 3:length(axes.Children)
axes.Children(s).LineWidth = 2;
axes.Children(s).LineStyle = '-';
end
hold on;
parallelcoords([temp1, temp2], 'Color', [0.7 0.7 0.7], 'LineStyle', '-','LineWidth',0.2,...
  'Marker', '.', 'MarkerSize', 10);
 title('Beta burst rate')
axes = gca;
axes.FontSize = 13;

%Spectrum
subplot(2,6,[5 6 11 12])
e = 1:96;%find(stat1.mask &   stat4.mask);
addpath ('/Users/juliorodriguezlarios/Documents/MATLAB/stdshade.m')
ga_spectrum = squeeze(nanmean(spectrum1(:,e,:),2));
a = stdshade(ga_spectrum,0.1,'green',frequencies,0);
a.LineWidth = 3;
hold on
ga_spectrum = squeeze(nanmean(spectrum2(:,e,:),2));
a = stdshade(ga_spectrum,0.1,[0.8500 0.3250 0.0980],frequencies,0);
a.LineWidth = 3;
legend('',condition1,'',condition2, 'location', 'northeast')
xlabel('Frequency (Hz)')
ylabel('Relative amplitude (zscore)')
axes = gca;
axes.FontSize = 13;
title('Mean burst spectrum')

 colormap parula



%% compare interindividual variability behavior with bursts

%define variables of interest (behavior and condition)

%if you want to test reaction time
clearvars behavior bamp bdur bfreq brate 
index = accuracy>-1;
for s = 1:size(files,1)
behavior(s,1) = mean(rt(s,index(s,:)),2);
bamp(s,:) = nanmean(mean_bamp_delay(s,:,index(s,:)),3);
bdur (s,:) = nanmean(mean_bdur_delay(s,:,index(s,:)),3);
bfreq (s,:) = nanmean(mean_bfreq_delay(s,:,index(s,:)),3);
brate (s,:) = nanmean(mean_brate_delay(s,:,index(s,:)),3);
end

% if you want to test accuracy
% behavior = mean(accuracy,2);
% bamp = nanmean(mean_bamp_delay,3);
% bdur = nanmean(mean_bdur_delay,3);
% bfreq = nanmean(mean_bfreq_delay,3);
% brate = nanmean(mean_brate_delay,3);

 restoredefaultpath
 addpath('/Users/juliorodriguezlarios/Documents/MATLAB/fieldtrip-master')

%put group data in fieldtrip structure

load('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/variables/ft_structure_fft.mat')
load('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/variables/foof_per_trial_laplacianfixed_justspectrum_nogaussian.mat', 'data1')

freq_baseline.label = data1.label;
freq_baseline.elec = data1.elec;
freq_baseline.freq = 1;

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_correlationT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg.computestat = 'yes';
% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freq_baseline);
 cfg.design   = behavior;
 cfg.type = 'Pearson';

%define var
freq_baseline.powspctrm = bamp;
%compare 
[stat1] = ft_freqstatistics(cfg, freq_baseline);
%plot 
sig_m = find(stat1.mask);
subplot(2,2,1)
topoplot(stat1.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
title('Beta burst amplitude')

% %define var
freq_baseline.powspctrm = bdur;
%compare 
[stat2] = ft_freqstatistics(cfg, freq_baseline);
%plot 
sig_m = find(stat2.mask);
subplot(2,2,2)
topoplot(stat2.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
title('Beta burst duration')

% %define var
freq_baseline.powspctrm = bfreq;
%compare 
[stat4] = ft_freqstatistics(cfg, freq_baseline);
%plot 
sig_m = find(stat4.mask);
subplot(2,2,3)
topoplot(stat4.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
title('Beta burst frequency')

% %define var
freq_baseline.powspctrm = brate;
compare 
[stat5] = ft_freqstatistics(cfg, freq_baseline);
subplot(2,2,5)
topoplot(stat5.stat,OUTEEG.chanlocs,'emarker2',{sig_m,'*','black',5,1})
caxis ([-6 6]);
title('Beta burst rate')



