%% load epoched clean data
%load ERP data of humans
load ('/Users/juliorodriguezlarios/Desktop/Projects/categorization EEG/variables/erp_data_2to2.mat')

%% add paths
addpath('C:\Users\bjjr014\Documents\MATLAB\eeglab_current\eeglab2023.1')
 eeglab;close
addpath ('C:\Users\bjjr014\Documents\MATLAB')
%% do group pca for humans

%estimate erp per subject and electrode excluding bad trials
for f = 1:24
erp(f,:,:) = squeeze(nanmean(EEGdata(f,:,:, bad_trials(f,:)==0),4));
end
%zscore
erp = zscore(erp,0,3);
%run pca
 temp = permute(erp,[2,3,1]);
 temp2 = reshape(temp,96,1000*24);
%restoredefaultpath 
[coeff,score,latent,tsquared,explained,mu] = pca(temp2');


%% get erp stats in humans

%%%%%%%%%% define spatial filter based on group pca
weights_central = repmat(coeff(:,2),1,24)';

%estimate time series from weights
%clearvars components
for s = 1:24
ntrials = length(find(~isnan(squeeze(EEGdata(s,1,1,:)))));
data = reshape(EEGdata(s,:,:,1:ntrials),1,96,1000*ntrials);
temp = squeeze(data);
components(s,:,1:ntrials) = reshape(weights_central(s,:) * temp, 1000, ntrials);  
end

% in case you want to select electrodes instead 
%  e = ismember(extractfield(OUTEEG.chanlocs, 'labels'), {'F1', 'F2', 'FC1', 'FC2','FCz', 'Fz', 'FFC1h', 'FFC2h'});
%components = squeeze(nanmean(EEGdata(:,e,:,:),2));

time2save = find(EEG_epoched.times==0):find(EEG_epoched.times==1996);

co_erp  = nan(24,1000);
inco_erp  = nan(24,1000);
co_erp_long  = nan(24,1000);
co_erp_short  = nan(24,1000);
inco_erp_long  = nan(24,1000);
inco_erp_short  = nan(24,1000);

clearvars co_erp inco_erp co_erp_long inco_erp_long co_erp_short inco_erp_short data1 data2
for f = 1:24
co_erp_long(f,:) = squeeze(nanmean(components(f,:,correct(f,:)==1 & bad_trials(f,:)==0 & long(f,:)==1 ),3));
inco_erp_long(f,:) =squeeze(nanmean(components(f,:,correct(f,:)==0 & bad_trials(f,:)==0 & long(f,:)==1),3));

co_erp_short(f,:) = squeeze(nanmean(components(f,:,correct(f,:)==1 & bad_trials(f,:)==0 & long(f,:)==0 ),3));
inco_erp_short(f,:) =squeeze(nanmean(components(f,:,correct(f,:)==0 & bad_trials(f,:)==0 & long(f,:)==0),3));
end

% Exclude subjects with bad performance
 acc = nanmean(correct,2);
index_exc = acc<0.6;
co_erp_long(index_exc,:,:) = [];
co_erp_short(index_exc,:,:) = [];
inco_erp_long(index_exc,:,:) = [];
inco_erp_short(index_exc,:,:) = [];

%permutation test
%load fieldtrip structure as data1
data = eeglab2fieldtrip(OUTEEG, 'raw');
data1.label = data.label(1);
%load fieldtrip structure as data1
data1.trialinfo = [];
data1.trial = [];
data1.time = [];
data1.fsample = EEG_epoched.srate;
data2 = data1;

%fill ft structure
for s = 1:size(inco_erp_long,1)
  data1.trial{s} = squeeze(co_erp_long(s,time2save));
  data2.trial{s} = squeeze(co_erp_short(s,time2save));

  data1.time{s} = EEG_epoched.times(time2save); 
  data2.time{s} = EEG_epoched.times(time2save); 
end

%permutation test
addpath('C:\Users\bjjr014\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118')
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
if length(data1.label)==96
cfg_neighb.method    = 'triangulation';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, data1);
end
subj = size(data1.trial,2);
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


[stat_co] = ft_timelockstatistics(cfg, data1, data2);


for s = 1:size(inco_erp_long,1)
  data1.trial{s} = squeeze(inco_erp_long(s,time2save));
  data2.trial{s} = squeeze(inco_erp_short(s,time2save));
  data1.time{s} = EEG_epoched.times(time2save); 
  data2.time{s} = EEG_epoched.times(time2save); 
end

[stat_inco] = ft_timelockstatistics(cfg, data1, data2);



%% get erp same block duration

temp1 = test_durations(1,:);
temp2 = test_durations(2,:);
temp3 = test_durations(3,:);

erp_b1 = nan(24,8,1000);
erp_b2= nan(24,8,1000);
erp_b3= nan(24,8,1000);
erp_b1_inco= nan(24,8,1000);
erp_b2_inco= nan(24,8,1000);
erp_b3_inco= nan(24,8,1000);

for s = 1:24
    for d = 1:length(temp1)
        erp_b1(s,d,:) = squeeze(nanmean(components(s,:,duration(s,:) == temp1(d)& bad_trials(s,:)==0 & correct(s,:)==1 & block_type(s,:)==1),3));
        erp_b2(s,d,:) = squeeze(nanmean(components(s,:,duration(s,:) == temp2(d)& bad_trials(s,:)==0 & correct(s,:)==1 & block_type(s,:)==2),3));
        erp_b3(s,d,:) = squeeze(nanmean(components(s,:,duration(s,:) == temp3(d)& bad_trials(s,:)==0 & correct(s,:)==1 & block_type(s,:)==3),3));

        erp_b1_inco(s,d,:) = squeeze(nanmean(components(s,:,duration(s,:) == temp1(d)& bad_trials(s,:)==0 & correct(s,:)==0 & block_type(s,:)==1),3));
        erp_b2_inco(s,d,:) = squeeze(nanmean(components(s,:,duration(s,:) == temp2(d)& bad_trials(s,:)==0 & correct(s,:)==0 & block_type(s,:)==2),3));
        erp_b3_inco(s,d,:) = squeeze(nanmean(components(s,:,duration(s,:) == temp3(d)& bad_trials(s,:)==0 & correct(s,:)==0 & block_type(s,:)==3),3));
   
    end
  
end


% permutation test ERP same duration
temp11 = squeeze(nanmean(nanmean(cat(4, erp_b1(:,7:8,:), erp_b2(:,7:8,:)) ,2),4));
temp22 = squeeze(nanmean(nanmean(cat(4, erp_b2(:,1:2,:), erp_b3(:,1:2,:)) ,2),4));
temp111 = squeeze(nanmean(nanmean(cat(4, erp_b1_inco(:,7:8,:), erp_b2_inco(:,7:8,:)) ,2),4));
temp222 = squeeze(nanmean(nanmean(cat(4, erp_b2_inco(:,1:2,:), erp_b3_inco(:,1:2,:)) ,2),4));


%fill ft structure
for s = 1:size(inco_erp_long,1)
  data1.trial{s} = temp11(s,time2save);
  data2.trial{s} =temp22(s,time2save);

  data1.time{s} = EEG_epoched.times(time2save); 
  data2.time{s} = EEG_epoched.times(time2save); 
end


%permutation test
%restoredefaultpath
%addpath('/Users/juliorodriguezlarios/Documents/MATLAB/fieldtrip-master')
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
%cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
if length(data1.label)==96
cfg_neighb.method    = 'triangulation';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, data1);
end
subj = size(data1.trial,2);
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


[stat_co_sameduration] = ft_timelockstatistics(cfg, data1, data2);

for s = 1:size(inco_erp_long,1)
  data1.trial{s} = temp111(s,time2save);
  data2.trial{s} =temp222(s,time2save);

end

[stat_inco_sameduration] = ft_timelockstatistics(cfg, data1, data2);


%% monkeys all blocks LFP

load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\monkey_data\ERP_m1_preSMA.mat')
time = time*1000;
% permutation test ERP
time2save = find(time==0):find(time==500);

%fill ft structure
for s = 1:size(long_correct,1)
  data1.trial{s} = squeeze(long_correct(s,time2save));
  data2.trial{s} = squeeze(short_correct(s,time2save));

  data1.time{s} =time(time2save); 
  data2.time{s} = time(time2save); 
end

%permutation test
addpath('C:\Users\bjjr014\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118')
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
subj = size(data1.trial,2);
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


[stat_co_monkey] = ft_timelockstatistics(cfg, data1, data2);


for s = 1:size(long_incorrect,1)
  data1.trial{s} = squeeze(long_incorrect(s,time2save));
  data2.trial{s} = squeeze(short_incorrect(s,time2save));

  data1.time{s} =time(time2save); 
  data2.time{s} = time(time2save); 
end

[stat_inco_monkey] = ft_timelockstatistics(cfg, data1, data2);


%% monkeys same duration 
load('C:\Users\bjjr014\Desktop\cattask paper\matlab variables\monkey_data\ERP_m1_preSMA_overlaps.mat')
time_samed = (-1:1/500:1)*1000;

time2save = find(time_samed ==0):find(time_samed ==500);

longsamedur= [ERP_TM1long; ERP_TM2long];
shortsamedur = [ERP_TM3short;ERP_TM2short];
data1.trial = [];
data2.trial = [];
data1.time = [];
data2.time = [];
%fill ft structure
for s = 1:size(longsamedur,1)
  data1.trial{s} = squeeze(longsamedur(s,time2save));
  data1.time{s} =time_samed (time2save); 
end
for s = 1:size(shortsamedur,1)
  data2.trial{s} = squeeze(shortsamedur(s,time2save));
    data2.time{s} = time_samed (time2save); 
end

%permutation test
addpath('C:\Users\bjjr014\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118')
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.design=[ones(1,size(longsamedur,1)), ones(1,size(shortsamedur,1)).*2];
cfg.ivar  = 1; 
cfg.uvar  = []; 


[stat_co_monkey_samedur] = ft_timelockstatistics(cfg, data1, data2);



%% monkeys all blocks firing rate

load('C:\Users\bjjr014\Desktop\cattask paper\monkey_data\firingrate_m1_preSMA.mat')
time = time*1000;
% permutation test ERP
time2save = find(time==0):find(time==500);

meansamples = 30; %no samples for moving mean

firrat_long_correct = movmean(zscore(firrat_long_correct,0,2),meansamples,2);
firrat_long_incorrect = movmean(zscore(firrat_long_incorrect,0,2),meansamples,2);
firrat_short_correct = movmean(zscore(firrat_short_correct,0,2),meansamples,2);
firrat_short_incorrect = movmean(zscore(firrat_short_incorrect,0,2),meansamples,2);


%fill ft structure
for s = 1:size(long_correct,1)
  data1.trial{s} = squeeze(firrat_long_correct(s,time2save));
  data2.trial{s} = squeeze(firrat_short_correct(s,time2save));

  data1.time{s} =time(time2save); 
  data2.time{s} = time(time2save); 
end


%permutation test
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

subj = size(data1.trial,2);
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


[stat_firrat_co_monkey] = ft_timelockstatistics(cfg, data1, data2);


data1.trial = [];
data2.trial = [];
data1.time = [];
data2.time = [];
for s = 1:size(long_incorrect,1)
  data1.trial{s} = squeeze(firrat_long_incorrect(s,time2save));
  data2.trial{s} = squeeze(firrat_short_incorrect(s,time2save));

  data1.time{s} =time(time2save); 
  data2.time{s} = time(time2save); 
end

subj = size(data1.trial,2);
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

[stat_firrat_inco_monkey] = ft_timelockstatistics(cfg, data1, data2);


%% firing ratemonkeys same duration
load('C:\Users\bjjr014\Desktop\cattask paper\monkey_data\firingrate_m1_preSMA_overlaps.mat')
time_samed = (-0.5:1/500:0.5)*1000;

time2save = find(time_samed ==0):find(time_samed ==500);

longsamedur_firrat= movmean(zscore([firrat_TM1_450_500_long; firrat_TM2_870_920_long],0,2),30,2);
shortsamedur_firrat = movmean(zscore([firrat_TM2_450_500_short;firrat_TM3_870_920_short],0,2),30,2);
data1.trial = [];
data2.trial = [];
data1.time = [];
data2.time = [];
%fill ft structure
for s = 1:size(longsamedur_firrat,1)
  data1.trial{s} = squeeze(longsamedur_firrat(s,time2save));
  data1.time{s} =time_samed (time2save); 
end
for s = 1:size(shortsamedur_firrat,1)
  data2.trial{s} = squeeze(shortsamedur_firrat(s,time2save));
    data2.time{s} = time_samed (time2save); 
end

%permutation test
addpath('C:\Users\bjjr014\Documents\MATLAB\fieldtrip-20230118\fieldtrip-20230118')
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.design=[ones(1,size(longsamedur_firrat,1)), ones(1,size(shortsamedur_firrat,1)).*2];
cfg.ivar  = 1; 
cfg.uvar  = []; 


[stat_co_monkey_samedur_firrat] = ft_timelockstatistics(cfg, data1, data2);


%% save everything
save('C:\Users\bjjr014\Desktop\cattask paper\alldata_andstats','-v7.3')


