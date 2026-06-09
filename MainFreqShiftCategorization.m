% Main code for the temporal categorization task while recording human EEG.
% This dataset and task were used in a published paper Eli, iScience, 2025 and Julio -PeerJ). The aim is to
% take preprocessed data, perform epoching, and then apply further analyses.
% For example, my approach is to examine inter-trial variability in beta frequency.
% In addition, other time windows can be considered, such as those related to
% working memory processes. NOTE: This analysis is based on preprocessed data.
% Author: Ehsan Eqlimi, columbia university
%Date: 8 June 2026
clc;clear;close all
%Add eeglab toolbox
addpath '/Users/ehsaneqlimi/eeglab2026.0.0'
eeglab;close

% Add fieldtrip toolbox
addpath('/Users/ehsaneqlimi/fieldtrip-20260518');
ft_defaults; % -> really important!

%Add my freq. analysis codes
% !ssh -T git@github.com
% system('git clone https://github.com//EhsanEqlimi/ComplexPCA-PD.git');
addpath(genpath(fullfile(pwd,'ComplexPCA-PD')));
% preprocessed data path
PreProcseedDataPath='/Users/ehsaneqlimi/FreqShiftDynData/preprocessed_EEG/';
% get file list
Dir=dir(fullfile(PreProcseedDataPath,'*.mat'));

% epoch relative to the start of the delay (stimulus offset)
delaystart = {'S101','S102','S103', 'S104', 'S105', 'S106','S107','S108',...
    'S201','S202','S203', 'S204', 'S205', 'S206','S207','S208',...
    'S 45','S 46','S 47', 'S 48', 'S 49', 'S 50','S 51','S 52'};
% These are the intervals presented as stimuli in seconds in three blocks ...
% (T1, T2, T3,). Each block contains 8 numbers: the first four correspond...
% to the short category, and the last four correspond to the long category.
test_durations = [0.2 0.250 0.319 0.331 0.369 0.381 0.450 0.5;...
    0.45 0.5 0.619 0.669 0.706 0.756 0.870 0.920;...ß
    0.870 0.920 0.981 1.169 1.231 1.419 1.470 1.520];
% EEGdat--> channel x time x trial
EEGdata = nan(length(Dir), 96,  550,  336);
bad_trials = nan(length(Dir),336);
block_type= nan(length(Dir),336);
duration= nan(length(Dir),336);
long= nan(length(Dir),336);
correct= nan(length(Dir),336);
rt= nan(length(Dir),336);


for f=2:length(Dir)

    %add load EEG data and loop for subject
    load([PreProcseedDataPath Dir(f).name]);
    %in case you want to delete more components
    % [component,type] = find(OUTEEG.etc.ic_classification.ICLabel.classifications(:,[2 3 4 6])>iclabel_threshold);
    % if ~isempty(unique(component))
    %  OUTEEG= pop_subcomp(OUTEEG, unique(component),0);
    % end
    %  %in case you want to interpolate more noisy channels
    %   OUTEEG.etc = [];
    %  clearvars OUTEEG2
    % OUTEEG2 = clean_channels(OUTEEG, cleanchannels_threshold);
    % OUTEEG = pop_interp(OUTEEG2,OUTEEG.chanlocs, 'spherical');


    [EEG_epoched, indices] = pop_epoch( OUTEEG, delaystart, [-2 5]);

    for e = 1:size(EEG_epoched.data,3) %loop epochs to get info

        %get events of each epoch
        %             temp_index = find(extractfield(EEG_epoched.event, 'epoch') ==e);
        %             temp_events = extractfield(EEG_epoched.event(temp_index), 'type');

        temp_index = find([EEG_epoched.event.epoch] == e);
        temp_events = {EEG_epoched.event(temp_index).type};
        %extract the name of delaystart (first number reflects block and third the
        %duration)
        temp = delaystart(ismember(delaystart,temp_events));

        %get block type and duration
        ind = ismember(delaystart,temp);
        ind = [ind(1:8); ind(9:16);ind(17:24)];

        [b,dur] = find(ind);

        block_type(f,e) = b;
        duration(f,e) = test_durations(block_type(f,e),dur);

        %get condition
        long(f,e) = logical(sum(duration(f,e) > median(test_durations(block_type(f,e),:))));
        %get accuracy
        if logical(find(ismember(temp_events,'S 81')))
            correct(f,e)= 1;
        end
        if logical(find(ismember(temp_events,'S 80')))
            correct(f,e)= 0;
        end

        %get rt
        %response mappings onset
        map_index = or(ismember(temp_events,'S 71'),ismember(temp_events,'S 72'));
        %motor response
        response_index = or(ismember(temp_events,'S190'),ismember(temp_events,'S188'));
        %get latencies and subtract
        %         temp_latencies = extractfield(EEG_epoched.event(temp_index), 'latency');
        temp_latencies = [EEG_epoched.event(temp_index).latency];
        try
            rt(f,e) = (temp_latencies(response_index) - temp_latencies(map_index)) / 500;
        catch
        end

    end


    %% epoching for investigating the delay window
    [EEG_epoched, indices] = pop_epoch( OUTEEG, delaystart, [-0.2 2]);

    %mark bad epochs
    [EEG_epoched,index] = pop_eegthresh(EEG_epoched, 1, 1:96, -100, ...
        100, 0, 2, 0, 0);

    bad_trials(f,1:size(EEG_epoched.data,3)) = EEG_epoched.reject.rejthresh;

    co_erp_long(f,:,:) = mean(EEG_epoched.data(:,:,correct(f,:)==1 & bad_trials(f,:)==0 & long(f,:)==1),3);
    co_erp_short(f,:,:) = mean(EEG_epoched.data(:,:,correct(f,:)==1 & bad_trials(f,:)==0 & long(f,:)==0),3);
    inco_erp_long(f,:,:) = mean(EEG_epoched.data(:,:,correct(f,:)==0 & bad_trials(f,:)==0 & long(f,:)==1),3);
    inco_erp_short(f,:,:) = mean(EEG_epoched.data(:,:,correct(f,:)==0 & bad_trials(f,:)==0 & long(f,:)==0),3);


    EEGdata(f,:,:,1:size(EEG_epoched.data,3)) =  EEG_epoched.data;
    %% eeglab spectopo (averaged over trial/epoch)
    EpochLength=size(EEG_epoched.data,2);
    Fs=EEG_epoched.srate;
    WinSize=EpochLength;% Fs*2; %samples
    NFFT=2.0 ^(ceil(log2(WinSize)));
    FOI=[1 40]; %beta
    figure,
    [Spec,Freq]=spectopo(EEG_epoched.data, 0, Fs,'nfft',NFFT, 'winsize',WinSize,...
        'overlap',0,'freqrange',FOI,'wintype','hamming');
    %Note: spectopo is 10log10​(power) (aka db)

    %plot  PSD
    figure;
    plot(Freq,Spec');
    xlabel('Frequency (Hz)');
    ylabel('Power (dB, corrected)');
    title(' Spectrum');
    xlim([2 40]);
    %% 1/f effect (simple log-log)
    %Note:% Computed PSD for each epoch, average across epochs (via spectopo), then fit/remove 1/f
    % (preferred over fitting 1/f separately for each epoch)
    FitIdx=Freq >= 2 & Freq <= 40;
    FreqSel=Freq(FitIdx);
    logF=log10(FreqSel);
    SpecSel=Spec(:,FitIdx)/10; % convert dB → log10(power)
    for Ch=1:size(SpecSel,1) %channel loop
        ChannelPower=SpecSel(Ch,:)'; %no need to log10 (spectopo is db)
        Powfitted=polyfit(logF,ChannelPower,1);
        AperiodicFit = polyval(Powfitted,logF);
        SpecFlat(:,Ch) = ChannelPower-AperiodicFit;
        Alpha(Ch)=-Powfitted(1);

    end
    %plot falttend PSD
    figure;
    plot(FreqSel,mean(SpecFlat',1), 'k', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB, corrected)');
    title('Flattened Spectrum');
    xlim([2 40]);


    %plot falttend PSD
    figure;
    plot(FreqSel,SpecFlat');
    xlabel('Frequency (Hz)');
    ylabel('Power (dB, corrected)');
    title('Flattened Spectrum');
    xlim([2 40]);
    %% my code
    Params=struct('Fs',Fs,'tapers', [2, 3],'Fpass',[1, 40],'itc',0,'pad',1);
    % create data: channel x time x trial 3-way
    EEG=EEG_epoched.data;
    %% MTPLV
    % No 1/f correction for PLV/ITC (phase-based measures); applied only to power.
    %PLV
    Params.itc=0;
    [PLV,SelFreqs]=FnMultiTaperFreqPLV(EEG,Params);
    figure;
    plot(SelFreqs,mean(PLV,1), 'k', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('PLV');
    xlim([1 40]);
    figure;
    plot(SelFreqs,PLV', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('PLV');
    xlim([1 40]);
    %ITC
    Params.itc=1;
    [PLV,SelFreqs]=FnMultiTaperFreqPLV(EEG,Params);
    figure;
    plot(SelFreqs,mean(PLV,1), 'k', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('ITC');
    xlim([1 40]);
    figure;
    plot(SelFreqs,PLV', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('ITC');
    xlim([1 40]);
    %% MTSpec
    [SpecTapFIND,SelFreqs]=FnMultiTaperFreqSpec(EEG,Params);
    figure;
    plot(SelFreqs,mean(SpecTapFIND,1), 'k', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Spec');
    title('PLV');
    xlim([1 40]);
    figure;
    plot(SelFreqs,10*log10(SpecTapFIND'), 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Spec');
    xlim([1 40]);
    %% Chronux
    S=[];
    params = struct( ...
        'Fs', Fs, ...
        'tapers', [2 3], ...
        'fpass', [1 40], ...
        'pad', 2, ...
        'trialave', 1 ...
        );
    for Ch=1:size(EEG,1)
        Data=squeeze(EEG(Ch,:,:));   % time*trial
        [S(Ch,:),freqChron]=mtspectrumc(Data, params);
    end
    figure;
    plot(freqChron,10*log10(S'), 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('ITC');
    xlim([1 40]);
    %%
    win=2.2 ;   % full epoch length = 550 sample
    EEGPerm=permute(EEG,[2,1,3]);
    [Sc, Cmat, Ctot, Cvec, Cent, fr] = CrossSpecMatc(EEGPerm, win, params);
    %  Sc	cross-spectral matrix (channels × channels × freq)
    % Cmat	coherence matrix (freq-resolved)
    % Ctot	total coherence
    % Cvec	vectorized coherence (for stats)
    % Cent	centered coherence
    % f	frequency vector
    % Assumes: f, Ctot, Cvec, Cent already exist

figure('Color','w','Position',[100 100 900 600]);

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Optional: colorblind-friendly paletteß
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.4660 0.6740 0.1880];

% ---------- Ctot ----------
nexttile
plot(fr, Ctot, 'LineWidth', 2.5, 'Color', c1);
grid on
box off
ylabel('C_{tot}')
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% % ---------- Cvec ----------
% nexttile
% plot(nanmean(abs(Cvec),1), 'LineWidth', 2.5, 'Color', c2);
% grid on
% box off
% ylabel('C_{vec}')
% set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% ---------- Cent ----------
nexttile
plot(fr, Cent, 'LineWidth', 2.5, 'Color', c3);
grid on
box off
ylabel('C_{ent}')
xlabel('Frequency (Hz)')
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% Global title
sgtitle('Cross-Spectral Coherence Measures','FontSize',14,'FontWeight','bold')


% long(f,:)==1
% [Sc, Cmat, Ctot, Cvec, Cent, f] = CrossSpecMatc(EEGPerm, win, params);

[~,idx] = min(abs(f-33));
figure,topoplot(abs(Cvec(idx,:)),OUTEEG.chanlocs);

EEGPermLong=permute(EEG(:,:,find(long(f,:)==1)),[2,1,3]);
[Sc, Cmat, Ctot, Cvec, Cent, fr]=CrossSpecMatc(EEGPermLong,win,params);
    %  Sc	cross-spectral matrix (channels × channels × freq)

figure('Color','w','Position',[100 100 900 600]);

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Optional: colorblind-friendly paletteß
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.4660 0.6740 0.1880];

% ---------- Ctot ----------
nexttile
plot(fr, Ctot, 'LineWidth', 2.5, 'Color', c1);
grid on
box off
ylabel('C_{tot}')
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% % ---------- Cvec ----------
% nexttile
% plot(nanmean(abs(Cvec),1), 'LineWidth', 2.5, 'Color', c2);
% grid on
% box off
% ylabel('C_{vec}')
% set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% ---------- Cent ----------
nexttile
plot(fr, Cent, 'LineWidth', 2.5, 'Color', c3);
grid on
box off
ylabel('C_{ent}')
xlabel('Frequency (Hz)')
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% Global title
sgtitle('Cross-Spectral Coherence Measures','FontSize',14,'FontWeight','bold')




EEGPermShort=permute(EEG(:,:,find(long(f,:)==0)),[2,1,3]);
[Sc, Cmat, Ctot, Cvec, Cent, fr]=CrossSpecMatc(EEGPermShort,win,params);
    %  Sc	cross-spectral matrix (channels × channels × freq)


    figure('Color','w','Position',[100 100 900 600]);

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Optional: colorblind-friendly paletteß
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.4660 0.6740 0.1880];

% ---------- Ctot ----------
nexttile
plot(fr, Ctot, 'LineWidth', 2.5, 'Color', c1);
grid on
box off
ylabel('C_{tot}')
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% % ---------- Cvec ----------
% nexttile
% plot(nanmean(abs(Cvec),1), 'LineWidth', 2.5, 'Color', c2);
% grid on
% box off
% ylabel('C_{vec}')
% set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% ---------- Cent ----------
nexttile
plot(fr, Cent, 'LineWidth', 2.5, 'Color', c3);
grid on
box off
ylabel('C_{ent}')
xlabel('Frequency (Hz)')
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% Global title
sgtitle('Cross-Spectral Coherence Measures','FontSize',14,'FontWeight','bold')


end