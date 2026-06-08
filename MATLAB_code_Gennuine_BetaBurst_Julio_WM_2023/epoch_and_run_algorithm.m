clc
clearvars
%add necessary paths
addpath ('/Users/juliorodriguezlarios/Documents/MATLAB/eeglab2022.0')
eeglab;close
addpath('/Users/juliorodriguezlarios/Documents/MATLAB/FieldTrip')
addpath('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/Paper/scripts')

%go to folder with preprocessed data 
cd '/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/preprocessed_EEG/';
files = dir('*.mat');

%define parameters for detect burst function
frequencies = 2:40;%frequencies to analyse and resolution
nmbcycles = 1; %minimum number of cycles
noiselimit = [-10 inf]; %beta power has to be above the mean and below x std after 1/f correction
lowerfreq = 15; %lowest frequency to save parameters


for subject = 1:size(files,1) %loop subjects

%load subject file
 load(files(subject).name, 'OUTEEG')
  
 %get number of rejected components
 rejected_c(subject,1) = 95 - size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
 
%get behavioral info after reordering
cload = [OUTEEG.load(OUTEEG.load==1)';OUTEEG.load(OUTEEG.load==3)'];
accuracy = [OUTEEG.accuracy(OUTEEG.load==1);OUTEEG.accuracy(OUTEEG.load==3)];
instruction = [OUTEEG.condition(OUTEEG.load==1)';OUTEEG.condition(OUTEEG.load==3)'];
reactiontime = [OUTEEG.rt(OUTEEG.load==1);OUTEEG.rt(OUTEEG.load==3)];

%resample to make computations faster
[OUTEEG] = pop_resample( OUTEEG, 250);

%epoching
%epoch length depends of condition
%trials with load 1 last 11 seconds:
%fixation(3), 1stcue(1), 1stdelay(3), 2ndcue(1),2nddelay(3)
%trials with load 3 last 2 seconds more (they got 2 more 1stcues)
clearvars epochs_l1 epochs_l3

%create index of events marking fixation for load 1 and 3
index_fix = ismember(extractfield(OUTEEG.event, 'type'), 'S 11');
index_load = zeros(1, length(index_fix));
index_load(index_fix==1) = OUTEEG.load;
%epoch load 1 and 3 trials separatedly
[epochs_l1, indices] = pop_epoch(OUTEEG, [], [0 11],'eventindices',find(index_load==1'));
[epochs_l3, indices] = pop_epoch(OUTEEG, [], [0 13],'eventindices',find(index_load==3'));

%%% combine l3 and l1 in fieldtrip; load 1 goes first%%
data1 = eeglab2fieldtrip( epochs_l1, 'raw', 'none' );
data3 = eeglab2fieldtrip( epochs_l3, 'raw', 'none' );
data_all = data1;
data_all.trial = [data1.trial, data3.trial];
data_all.time = [data1.time, data3.time];
data_all.trialinfo = [data1.trialinfo; data3.trialinfo];

%update trialinfo in fieldtrip structure
data_all.trialinfo.load(1:192) = cload;
data_all.trialinfo.accuracy(1:192) = accuracy;
data_all.trialinfo.instruction(1:192) = instruction;
data_all.trialinfo.rt(1:192) = reactiontime;

%%%%% epoch data %%%%%
cfg = [];
cfg.toilim = [0 3];
data_baseline = ft_redefinetrial(cfg, data_all);

cfg = [];
%different time points for delay1 in load 1 and 3
cfg.toilim = [repmat([4 7],size(data1.trial,2),1); repmat([6 9],size(data3.trial,2),1)];
data_1stdelay = ft_redefinetrial(cfg, data_all);

cfg = [];
%different time points for delay2 in load 1 and 3
cfg.toilim = [repmat([8 11],size(data1.trial,2),1); repmat([10 13],size(data3.trial,2),1)];
data_2nddelay = ft_redefinetrial(cfg, data_all);

%run detect burst function in each time period

[btime_baseline(subject,:,:,:),bamp_baseline(subject,:,:,:), bnumb_baseline(subject,:,:,:), bdur_baseline(subject,:,:,:), allbursts_spectrum_baseline(subject,:,:,:),allbursts_shape_baseline(subject,:,:,:),btime_lowfreq_baseline(subject,:,:,:),...
   bamp_lowfreq_baseline(subject,:,:,:), bnumb_lowfreq_baseline(subject,:,:,:), bdur_lowfreq_baseline(subject,:,:,:), allbursts_spectrum_lowfreq_baseline(subject,:,:,:),allbursts_shape_lowfreq_baseline(subject,:,:,:)]...
    = detect_bursts_maxpeak(data_baseline,frequencies,nmbcycles,noiselimit,lowerfreq);

[btime_1stdelay(subject,:,:,:),bamp_1stdelay(subject,:,:,:), bnumb_1stdelay(subject,:,:,:), bdur_1stdelay(subject,:,:,:), allbursts_spectrum_1stdelay(subject,:,:,:),allbursts_shape_1stdelay(subject,:,:,:),btime_lowfreq_1stdelay(subject,:,:,:),...
   bamp_lowfreq_1stdelay(subject,:,:,:), bnumb_lowfreq_1stdelay(subject,:,:,:), bdur_lowfreq_1stdelay(subject,:,:,:), allbursts_spectrum_lowfreq_1stdelay(subject,:,:,:),allbursts_shape_lowfreq_1stdelay(subject,:,:,:)]...
    = detect_bursts_maxpeak(data_1stdelay,frequencies,nmbcycles,noiselimit,lowerfreq);

[btime_2nddelay(subject,:,:,:),bamp_2nddelay(subject,:,:,:), bnumb_2nddelay(subject,:,:,:), bdur_2nddelay(subject,:,:,:), allbursts_spectrum_2nddelay(subject,:,:,:),allbursts_shape_2nddelay(subject,:,:,:),btime_lowfreq_2nddelay(subject,:,:,:),...
   bamp_lowfreq_2nddelay(subject,:,:,:), bnumb_lowfreq_2nddelay(subject,:,:,:), bdur_lowfreq_2nddelay(subject,:,:,:), allbursts_spectrum_lowfreq_2nddelay(subject,:,:,:),allbursts_shape_lowfreq_2nddelay(subject,:,:,:)]...
    = detect_bursts_maxpeak(data_2nddelay,frequencies,nmbcycles,noiselimit,lowerfreq);


end


save('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/variables/detectburstfinal_2to40','-v7.3')







