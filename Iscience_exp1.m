%% timelock and trial select

% load ~epx1_instfreq_decdelay.mat'

for si = 1:numel(decdelay_instfreq)
    
   trl = decdelay_instfreq{si}.trialinfo;
   
%   F1_trls = find(trl(:,5)==1 & trl(:,1)==1 & trl(:,4)==1); % correct short
 %  F2_trls = find(trl(:,5)==1 & trl(:,1)==1 & trl(:,4)==2); % correct long
 
    F1_trls = find( (trl(:,3)==1 & trl(:,4)==1) | (trl(:,3)==2 & trl(:,4)==0) ); % subjective short
    F2_trls = find( (trl(:,3)==2 & trl(:,4)==1) | (trl(:,3)==1 & trl(:,4)==0) ); % correct long
 
% F1_trls = find(trl(:,3)==1 & trl(:,4)==1); % correct short
% F2_trls = find(trl(:,3)==2 & trl(:,4)==1); % correct long

   min_trls = min(numel(F1_trls), numel(F2_trls));
   
   F1 = randsample(F1_trls, min_trls);
   F2 = randsample(F2_trls, min_trls);
   
%    cfg = [];
%    cfg.demean = 'yes';
%    cfg.baselinewindow = [-.3 0];
%    decdelay_instfreqsi = ft_preprocessing(cfg, decdelay_instfreq{si})

   
   cfg = []; 
   cfg.latency = [-.5 2];
   cfg.trials = F1;
   F1_if{si} = ft_timelockanalysis(cfg, decdelay_instfreq{si});
   cfg.trials = F2;
   F2_if{si} = ft_timelockanalysis(cfg, decdelay_instfreq{si});
   clear decdelay_instfreqsi
%    cfg = [];
%    cfg.baseline = [-.1 0];
%    F1_if{si} = ft_timelockbaseline(cfg, F1_if{si})
%    F2_if{si} = ft_timelockbaseline(cfg, F2_if{si})
   
   si
end

%% GA

cfg = [];
%cfg.keepindividual = 'yes';
GAF1 = ft_timelockgrandaverage(cfg, F1_if{:});
GAF2 = ft_timelockgrandaverage(cfg, F2_if{:});

%% plot

figure;
plot(GAF1.time, GAF1.avg(2,:))
hold on; 
plot(GAF2.time, GAF2.avg(2,:))
xlim([-.3 2.2])

%% stat

label = F1_if{1}.label;

% neighbours(1).label = label{1};
% neighbours(1).neighblabel = {[label{2}, '; ' label{3}, '; ' label{4}]};
% neighbours(2).label = label{2};
% neighbours(2).neighblabel = {[label{1}, '; ' label{3}, '; ' label{4}]};
% neighbours(3).label = label{3};
% neighbours(3).neighblabel = {[label{1}, '; ' label{2}, '; ' label{4}]};
% neighbours(4).label = label{4};
% neighbours(4).neighblabel = {[label{1}, '; ' label{2}, '; ' label{3}]};

neighbours(1).label = label{1};
neighbours(1).neighblabel = label(3);
neighbours(2).label = label{2};
neighbours(2).neighblabel = label(4);
neighbours(3).label = label{3};
neighbours(3).neighblabel = label(1);
neighbours(4).label = label{4};
neighbours(4).neighblabel = label(2);

%
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;

%cfg.frequency = [4 30];
cfg.latency = [0 2];

cfg.clusterstatistic = 'wcm';          
cfg.tail = 0;               % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.05;               % alpha level of the permutation test
cfg.numrandomization = 5000;

% design
ll = numel(F1_if);
design = [];
design = repmat(1:ll, 1, 2);
design(2,:) = [repmat(1, 1, ll) repmat(2, 1, ll)];
cfg.design = design; 
%cfg.minnbchan = 1;

cfg.ivar  = 2;
cfg.uvar = 1; 
cfg.neighbours = neighbours;

cfg.channel = [1:4];

stat = ft_timelockstatistics(cfg, F1_if{:}, F2_if{:}) % negclu
%stat = ft_timelockstatistics(cfg, GAF1, GAF2)
% subj decision:
% p=4e-4
% 
%% get cluster values

clustertimes = stat.time(stat.negclusterslabelmat(2,:)==1);
start_idx = find(stat.time==clustertimes(1));
end_idx = find(stat.time==clustertimes(end));

for si=1:numel(F1_if)
    
ind_diff(si) = mean(F1_if{si}.avg(2, start_idx:end_idx))-mean(F2_if{si}.avg(2, start_idx:end_idx));

end
%% plot
timevec = [-2:1/250:1.996];
figure;
plot(timevec, squeeze(nanmean(long_correct_if(:,4,:),1)))
hold on;
plot(timevec, squeeze(nanmean(short_correct_if(:,4,:),1)))
xlim([-.5 1.8])


%% burst detection
addpath /project/3035003.01/JURIQUILLA/toolbox/SpectralEvents-master/SpectralEvents-master/
tic
clear
src_data_folder = '/project/3015079.02/categorization EEG/elie/svs_data_500Hz/';
src_data_files = dir(fullfile(src_data_folder, '*mat'));

d=1
for si = [1:numel(src_data_files)] % 7 is nan
    
    load([src_data_folder src_data_files(si).name])

% cfg = [];
% cfg.derivative = 'yes';
% svs_data = ft_preprocessing(cfg, svs_data);

trl = []; trl = svs_data.trialinfo;

   F1_trls = find(trl(:,3)==1 & trl(:,4)==1); % correct short
   F2_trls = find(trl(:,3)==2 & trl(:,4)==1); % correct long

% min_trls = min(numel(F1_trls), numel(F2_trls));
%    
% F1 = randsample(F1_trls, min_trls);
% F2 = randsample(F2_trls, min_trls);

cfg = [];
cfg.latency = [0 2];
cfg.channel = 2;
%cfg.trials = [F1 F2];
%cfg.avgoverchan = 'yes';
data = ft_selectdata(cfg, svs_data);

% update trl
%trl = data.trialinfo;

match = zeros(size(trl,1),1);
match(trl(:,3)==1 & trl(:,4)==1)=1;
match(trl(:,3)==2 & trl(:,4)==1)=2;

% put data in time x trials matrix
for ti = 1:numel(data.trial)
    
    x{d}(ti, :) = cell2mat(data.trial(ti));
    
end

x{d} = x{d}';
% long vs short 
classLabels{d} = match; clear match
%classLabels{1} = 1;

d = d+1;
end

eventBand = [13,35]; % freq range of bursts
fVec = 12:.5:36; % freqs for TFR
Fs = 500; % sampling rate
findMethod = 1;
vis = false; % visualize

[specEvents, TFRs, timeseries] = spectralevents(eventBand, fVec, Fs, findMethod, vis, x, classLabels);

toc
%% extract params

for fi = 1:numel(specEvents)
    
    maxfreq_short(fi) = mean(specEvents(fi).Events.Events.maximafreq(specEvents(fi).Events.Events.classLabels==1));
    maxfreq_long(fi) = mean(specEvents(fi).Events.Events.maximafreq(specEvents(fi).Events.Events.classLabels==2));
    
    burstrate_short(fi) = sum(specEvents(fi).Events.Events.classLabels==1) / sum(specEvents(fi).TrialSummary.TrialSummary.classLabels==1);
    burstrate_long(fi) = sum(specEvents(fi).Events.Events.classLabels==2) / sum(specEvents(fi).TrialSummary.TrialSummary.classLabels==2);
    
    Fspan_short(fi) = mean(specEvents(fi).Events.Events.Fspan(specEvents(fi).Events.Events.classLabels==1));
    Fspan_long(fi) = mean(specEvents(fi).Events.Events.Fspan(specEvents(fi).Events.Events.classLabels==2));
    
    maxtiming_short(fi) = mean(specEvents(fi).Events.Events.maximatiming(specEvents(fi).Events.Events.classLabels==1));
    maxtiming_long(fi) = mean(specEvents(fi).Events.Events.maximatiming(specEvents(fi).Events.Events.classLabels==2));
    
    duration_short(fi) = mean(specEvents(fi).Events.Events.duration(specEvents(fi).Events.Events.classLabels==1));
    duration_long(fi) = mean(specEvents(fi).Events.Events.duration(specEvents(fi).Events.Events.classLabels==2));
    
    maxpow_short(fi) = mean(specEvents(fi).Events.Events.maximapower(specEvents(fi).Events.Events.classLabels==1));
    maxpow_long(fi) = mean(specEvents(fi).Events.Events.maximapower(specEvents(fi).Events.Events.classLabels==2));
    
%     numevents_short(fi) = mean(specEvents(fi).TrialSummary.TrialSummary.eventnumber(specEvents(fi).TrialSummary.TrialSummary.classLabels==1));
%     numevents_long(fi) = mean(specEvents(fi).TrialSummary.TrialSummary.eventnumber(specEvents(fi).TrialSummary.TrialSummary.classLabels==2));
    % this was same as burst rate
end

% maxfreq is significant with findMethod=1
% even moreso with findMethod = 2 (T=4)
% less so but stsill sig with method 3 (T=2.3, p=.03)

% to report for revision
% mean burst rate short 3.00 +/- .300 (T=48.9, p=0)
% mean burst rate long 2.87 +/- .384 (T=36.7, p=0)
%% tfr

clear
src_data_folder = '/project/3015079.02/categorization EEG/elie/svs_data/';
src_data_files = dir(fullfile(src_data_folder, '*mat'));

for fi = 1:numel(src_data_files)
    
    load([src_data_folder src_data_files(fi).name])
    
    trl = svs_data.trialinfo;
   
trl = []; trl = svs_data.trialinfo;

F1_trls = find(trl(:,3)==1 & trl(:,4)==1); % correct match
F2_trls = find(trl(:,3)==2 & trl(:,4)==1); % correct mismatch
   
   min_trls = min(numel(F1_trls), numel(F2_trls));
   
   F1 = randsample(F1_trls, min_trls);
   F2 = randsample(F2_trls, min_trls);
      
    cfg = [];
    cfg.latency = [-.1 2];
    data = ft_selectdata(cfg, svs_data)

    cfg = [];
    cfg.method = 'wavelet';
%    cfg.output = 'fractal';
 %   cfg.taper = 'hanning';
    cfg.foi = [4:36];
    cfg.pad = 3;
    cfg.toi = 0:.1:2;
 %   cfg.t_ftimwin = ones(1, length(cfg.foi))*.4;
    cfg.channel = [1 2 3 4];
    cfg.keeptrials = 'yes';
    
    cfg.trials = F1; % select F1 mot trials
    F1_fft{fi} = ft_freqanalysis(cfg, data)
    F1_fft{fi}.powspctrm = log10(F1_fft{fi}.powspctrm);
    F1_fft{fi} = ft_freqdescriptives([], F1_fft{fi});
    
    cfg.trials = F2; % select F1 aud trials
    F2_fft{fi} = ft_freqanalysis(cfg, data)
    F2_fft{fi}.powspctrm = log10(F2_fft{fi}.powspctrm);
    F2_fft{fi} = ft_freqdescriptives([], F2_fft{fi});
    
end

% GA
GAF1 = ft_freqgrandaverage([], F1_fft{:})
GAF2 = ft_freqgrandaverage([], F2_fft{:})

%% plt

figure;
plot(GAF1.freq, GAF1.powspctrm(4,:))
hold on; plot(GAF2.freq, GAF2.powspctrm(4,:))

%% stat
% this way of making neighbours might not be working
neighbours = [];
label = GAF1.label;

% neighbours(1).label = label{1};
% neighbours(1).neighblabel = {[label{2}, '; ' label{3}, '; ' label{4}]};
% neighbours(2).label = label{2};
% neighbours(2).neighblabel = {[label{1}, '; ' label{3}, '; ' label{4}]};
% neighbours(3).label = label{3};
% neighbours(3).neighblabel = {[label{1}, '; ' label{2}, '; ' label{4}]};
% neighbours(4).label = label{4};
% neighbours(4).neighblabel = {[label{1}, '; ' label{2}, '; ' label{3}]};

neighbours(1).label = label{1};
neighbours(1).neighblabel = label(3);
neighbours(2).label = label{2};
neighbours(2).neighblabel = label(4);
neighbours(3).label = label{3};
neighbours(3).neighblabel = label(1);
neighbours(4).label = label{4};
neighbours(4).neighblabel = label(2);
%
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';

cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;

cfg.frequency = [13 35];
cfg.latency = [0 2];

cfg.clusterstatistic = 'maxsum';          
cfg.tail = 0;               % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.05;               % alpha level of the permutation test
cfg.numrandomization = 10000;

% design
ll = numel(F1_fft);
design = [];
design = repmat(1:ll, 1, 2);
design(2,:) = [repmat(1, 1, ll) repmat(2, 1, ll)];
cfg.design = design; 

cfg.ivar  = 2;
cfg.uvar = 1; 
cfg.neighbours = neighbours;
%cfg.channel = 5;

stat = ft_freqstatistics(cfg, F1_fft{:}, F2_fft{:}) 

%% plot stat
poscluster = stat.posclusterslabelmat==1;
negcluster = stat.negclusterslabelmat==1;

figure; imagesc(stat.time, stat.freq, squeeze(poscluster(2,:,:))); axis xy 
figure; imagesc(stat.time, stat.freq, squeeze(negcluster(4,:,:))); axis xy 

%% decoding

clear
addpath /project/3035003.01/MVPA-Light-master/startup
startup_MVPA_Light

load('/project/3015079.02/categorization EEG/elie/inst_freq_sens/instfreq_decdelay.mat')

for fi = 1:numel(decdelay_instfreq)
    

% load dataset from 1 subj to try
trl = []; trl = decdelay_instfreq{fi}.trialinfo;

F1_trls = find(trl(:,3)==1 & trl(:,4)==1); % correct short
F2_trls = find(trl(:,3)==2 & trl(:,4)==1); % correct long

   cfg = [];
   cfg.latency = [-.4 1.9];
cfg.trials = F1_trls;
data_F1 = ft_selectdata(cfg, decdelay_instfreq{fi})
cfg.trials = F2_trls;
data_F2 = ft_selectdata(cfg, decdelay_instfreq{fi})

cfg = [];
cfg.method           = 'mvpa';
cfg.features         = 'chan';
%cfg.features         = [];
cfg.mvpa.classifier  = 'lda'; % or lda
cfg.mvpa.metric      = 'auc';
cfg.mvpa.k           = 8;
cfg.mvpa.repeat      = 2;
%cfg.neighbours        = neighbours;

cfg.design = [ones(numel(F1_trls),1); 2*ones(numel(F2_trls),1)];
cfg.mvpa.preprocess = 'zscore';

statx{fi} = ft_timelockstatistics(cfg, data_F1, data_F2)

fi
end

%% average the stat?

% below chance ! -> do on sens level (it works there)

for fi = 1:24
    
    auc(fi, :) = statx{fi}.auc;
    
end

figure; plot(statx{1}.time, smooth(mean(auc), 7))

%% do we see freqshift on the fft spectra?
clear
folder = '/project/3015079.02/categorization EEG/elie/svs_data_dec_2.5/';
files = dir(fullfile(folder, '*mat'))

for fi = 1:numel(files)
    
    load([folder files(fi).name])
    
    trl = []; trl = svs_data.trialinfo;

    F1_trls = find(trl(:,3)==1 & trl(:,4)==1); % correct short
    F2_trls = find(trl(:,3)==2 & trl(:,4)==1);
    
    min_trls = min(numel(F1_trls), numel(F2_trls));
   
   F1 = randsample(F1_trls, min_trls);
   F2 = randsample(F2_trls, min_trls);
   
   cfg = [];
   cfg.derivative = 'yes';
   svs_data = ft_preprocessing(cfg, svs_data)
          
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
%    cfg.output = 'fooof_peaks';
%     cfg.taper = 'dpss';
%     cfg.tapsmofrq = 2;
    cfg.pad = 4;
    cfg.foilim = [8 35];
    
    cfg.trials = F1;
    fft_short{fi} = ft_freqanalysis(cfg, svs_data)
    
    cfg.trials = F2;
    fft_long{fi} = ft_freqanalysis(cfg, svs_data)
     
end

%% GA

GA_short = ft_freqgrandaverage([], fft_short{:})
GA_long = ft_freqgrandaverage([], fft_long{:})

%% plot

figure; plot(GA_short.freq, GA_short.powspctrm(2,:))
hold on; plot(GA_long.freq, GA_long.powspctrm(2,:))
xlim([8 35])