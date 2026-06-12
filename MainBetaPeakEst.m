clc;clear;
%close
% This is the main file to estimate the beta peaks using the different
% methods in temporal categorization task within the decision delay window
% based on human multichannel EEG signals
%% Initialization
% Add eeglab toolbox
addpath '/Users/ehsaneqlimi/eeglab2026.0.0'
eeglab;close
% Add fieldtrip toolbox
addpath('/Users/ehsaneqlimi/fieldtrip-20260518');
ft_defaults; % -> really important!
% Add my freq. analysis codes
% !ssh -T git@github.com
% system('git clone https://github.com//EhsanEqlimi/ComplexPCA-PD.git');
addpath(genpath(fullfile(pwd,'ComplexPCA-PD')));

% Prreprocessed data path (human EEG with temporal categorization task)
PreProcseedDataPath='/Users/ehsaneqlimi/FreqShiftDynData/preprocessed_EEG/';
% Get file list
Dir=dir(fullfile(PreProcseedDataPath,'S*.mat'));
%% The following lines (whole section) are taken from Julio's code on task parameters (changed to PascalCase for consistency)
% + some cosmetic changes
% epoch relative to the start of the delay (stimulus offset)
DelayStart={'S101','S102','S103', 'S104', 'S105', 'S106','S107','S108',...
    'S201','S202','S203', 'S204', 'S205', 'S206','S207','S208',...
    'S 45','S 46','S 47', 'S 48', 'S 49', 'S 50','S 51','S 52'};
% These are the intervals presented as stimuli in seconds in three blocks ...
% (T1, T2, T3,). Each block contains 8 numbers: the first four correspond...
% to the short category, and the last four correspond to the long category.
TestDurations=[0.2 0.250 0.319 0.331 0.369 0.381 0.450 0.5;...
    0.45 0.5 0.619 0.669 0.706 0.756 0.870 0.920;...
    0.870 0.920 0.981 1.169 1.231 1.419 1.470 1.520];
% EEGData--> channel x time x trial
EEGData=nan(length(Dir),96,550,336);
BadTrials = nan(length(Dir),336);
BlockType= nan(length(Dir),336);
Duration= nan(length(Dir),336);
Long= nan(length(Dir),336);
Correct= nan(length(Dir),336);
RT= nan(length(Dir),336);
%%
for i=2:length(Dir)
    % Add load EEG data and loop for subject
    load([PreProcseedDataPath Dir(i).name]);
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

    [EEGEpoched,Indices]=pop_epoch(OUTEEG,DelayStart, [-2 5]);

    for e=1:size(EEGEpoched.data,3) %loop epochs to get info

        % Get events of each epoch
        %             temp_index = find(extractfield(EEGEpochedd.event, 'epoch') ==e);
        %             temp_events = extractfield(EEGEpochedd.event(temp_index), 'type');

        TempIndex=find([EEGEpoched.event.epoch] == e);
        TempEvents={EEGEpoched.event(TempIndex).type};
        %extract the name of delaystart (first number reflects block and third the
        %duration)
        Temp=DelayStart(ismember(DelayStart,TempEvents));

        %get block type and duration
        Inds=ismember(DelayStart,Temp);
        Inds=[Inds(1:8);Inds(9:16);Inds(17:24)];

        [B,Dur]=find(Inds);

        BlockType(i,e)=B;
        Duration(i,e)=TestDurations(BlockType(i,e),Dur);

        % Get condition
        Long(i,e)=logical(sum(Duration(i,e)>median(TestDurations(BlockType(i,e),:))));
        % get accuracy
        if logical(find(ismember(TempEvents,'S 81')))
            Correct(i,e)= 1;
        end
        if logical(find(ismember(TempEvents,'S 80')))
            Correct(i,e)= 0;
        end

        % Get rt
        % Response mappings onset
        MapIndex=or(ismember(TempEvents,'S 71'),ismember(TempEvents,'S 72'));
        % Motor response
        ResponseIndex=or(ismember(TempEvents,'S190'),ismember(TempEvents,'S188'));
        % Get latencies and subtract
        %         temp_latencies = extractfield(EEGEpochedd.event(temp_index), 'latency');
        TemoLatencies=[EEGEpoched.event(TempIndex).latency];
        try
            RT(i,e)=(TemoLatencies(ResponseIndex)-TempLatencies(MapIndex))/500; %original smapling freq=500;
        catch
        end

    end


    %% epoching for investigating the delay window
    [EEGEpoched,Indices]=pop_epoch(OUTEEG,DelayStart,[-0.2 2]);

    % Mark bad epochs
    [EEGEpoched,Index]=pop_eegthresh(EEGEpoched,1,1:96,-100, ...
        100,0,2,0,0);

    BadTrials(i,1:size(EEGEpoched.data,3))=EEGEpoched.reject.rejthresh;

    CoERPLong(i,:,:)=mean(EEGEpoched.data(:,:,Correct(i,:)==1 & BadTrials(i,:)==0 & Long(i,:)==1),3);
    CoERPShort(i,:,:)=mean(EEGEpoched.data(:,:,Correct(i,:)==1 & BadTrials(i,:)==0 & Long(i,:)==0),3);
    InCoERPLong(i,:,:)=mean(EEGEpoched.data(:,:,Correct(i,:)==0 & BadTrials(i,:)==0 & Long(i,:)==1),3);
    IncoERPShort(i,:,:)=mean(EEGEpoched.data(:,:,Correct(i,:)==0 & BadTrials(i,:)==0 & Long(i,:)==0),3);

    %%EE
    LongCondBin(i,:)= BadTrials(i,:)==0 & Long(i,:)==1;
    ShortCondBin(i,:)=BadTrials(i,:)==0 & Long(i,:)==0;

    LongCondCorrectBin(i,:)=BadTrials(i,:)==0 & Long(i,:)==1 & Correct(i,:)==1;
    LongCondInCorrectBin(i,:)=BadTrials(i,:)==0 & Long(i,:)==1 & Correct(i,:)==0;


    ShortCondCorrectBin(i,:)=BadTrials(i,:)==0 & Long(i,:)==0 & Correct(i,:)==1;
    ShortCondInCorrectBin(i,:)=BadTrials(i,:)==0 & Long(i,:)==0 & Correct(i,:)==0;


    SubjectiveLongBin(i,:)=LongCondCorrectBin(i,:) | ShortCondInCorrectBin(i,:);

    SubjectiveShortBin(i,:)= ShortCondCorrectBin(i,:) | LongCondInCorrectBin(i,:);






    %% EEGlab spectopo based on welch (averaged over trial/epoch)
    [ChannNum,EpochLength,TrialNum]=size(EEGEpoched.data); % channel by samples/time by trail/epoch
    ChanDim=1;TimeDim=1;TrailDim=1;
    Fs=EEGEpoched.srate;
    WinSize=EpochLength;% Fs*2; %samples
    NFFT=2.0 ^(ceil(log2(WinSize)));
    FOI=[1 40]; %beta
    figure,
    [Spec,Freq]=spectopo(EEGEpoched.data, 0, Fs,'nfft',NFFT, 'winsize',WinSize,...
        'overlap',0,'freqrange',FOI,'wintype','hamming');
    % Note: spectopo is 10log10​(power) (aka db)
    %% Remove 1/f effect (based on a simple log-log approach)
    % Note: Compute PSD for each epoch, then average across epochs (via spectopo), and fit/remove the 1/f component.
    % This is preferred over fitting 1/f separately for each epoch.)
    FitIdx=Freq >= 2 & Freq <= 40;
    FreqSel=Freq(FitIdx);
    LogFreq=log10(FreqSel);
    SpecSel=Spec(:,FitIdx)/10; % convert dB → log10(power)
    for Ch=1:size(SpecSel,1) % channel loop
        ChannelPower=SpecSel(Ch,:)'; %no need to log10 (spectopo is db)
        FitParam=polyfit(LogFreq,ChannelPower,1);
        AperiodicFit=polyval(FitParam,LogFreq);
        FlatPower(:,Ch)=ChannelPower-AperiodicFit;% Note:FlatPower is already a logarithmic quantity
        Alpha(Ch)=-FitParam(1); %kinda scaling exponent
    end
    % Optional: you can use this function, same as above lines, just in case
    % of calling it in other main files :)
    % %%%%%%%%FitParam,AperiodicFit,FlatPower,Alpha]=FnOneOverFLogLog(PowerMat,FreqVec,FOI,IsLogPower)%%%%%%%%%%%%

    % Plot flat power
    FlatPower_dB = 10*FlatPower;
    figure('Color','w','Position',[100 100 700 500]);
    hold on
    plot(FreqSel,FlatPower_dB,...
        'Color',[0.8 0.8 0.8],...
        'LineWidth',0.5);
    plot(FreqSel,mean(FlatPower_dB,2),...
        'k','LineWidth',3);

    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    title('1/f-removed Spectrum')

    xlim([2 40])
    set(gca,...
        'FontSize',12,...
        'LineWidth',1.2,...
        'Box','off',...
        'TickDir','out')
    % Topoplot- spectral exponent
    figure,topoplot(Alpha,EEGEpoched.chanlocs);
    colormap(flipud(ft_colormap('RdBu')))
    title('Spectral exponent (\alpha)')
%     %% Frequency analysis based on my 2018's codes (multi-taper freq analysis)
%     %Add my freq. analysis codes
%     % !ssh -T git@github.com
%     % system('git clone https://github.com//EhsanEqlimi/ComplexPCA-PD.git');
%     addpath(genpath(fullfile(pwd,'ComplexPCA-PD')));
%     Params=struct('Fs',Fs,'tapers', [2, 3],'Fpass',[1, 40],'itc',0,'pad',1);
%     % create data: channel x time x trial 3-way
      EEG=EEGEpoched.data;
%     %% MTPLV (multi-taper)
%     % No 1/f correction for PLV/ITC (phase-based measures); applied only to power.
%     %--->PLV
%     Params.itc=0;
%     [PLV,SelFreqs]=FnMultiTaperFreqPLV(EEG,Params);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     plot(SelFreqs,PLV,...
%         'Color',[0.8 0.8 0.8],...
%         'LineWidth',0.5);
%     plot(SelFreqs,mean(PLV,1),...
%         'k','LineWidth',3);
% 
%     xlabel('Frequency (Hz)')
%     ylabel('PLV')
% 
%     xlim([2 40])
% 
%     %---->ITC
%     Params.itc=1;
%     [ITC,SelFreqs]=FnMultiTaperFreqPLV(EEG,Params);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     plot(SelFreqs,ITC,...
%         'Color',[0.8 0.8 0.8],...
%         'LineWidth',0.5);
%     plot(SelFreqs,mean(ITC,1),...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('ITC')
% 
%     xlim([2 40])
%     % --> MTSpec
%     [SpecTapFIND,SelFreqs]=FnMultiTaperFreqSpec(EEG,Params);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     plot(SelFreqs,SpecTapFIND,...
%         'Color',[0.8 0.8 0.8],...
%         'LineWidth',0.5);
%     plot(SelFreqs,mean(SpecTapFIND,1),...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('MT-Spec')
%     xlim([2 40])
%     %--> CPCA-Time Domain
%     Params.itc=0;
%     [YCPC,CWTS,CWTSfIND,SelFreqs,CPCFreqfIND]=FnMTCPCATIMEDomain(EEG,Params);
%     figure,plot(YCPC);
%     title(['YCPC']);
%     xlabel('Sample No.');
%     ylabel('YCPC');
% 
% 
%     figure,plot(SelFreqs,abs(CPCFreqfIND));
%     title(['CPCFreq-' Name(1:end-4)]);
%     xlabel('Freq (Hz)');
%     ylabel('CPCFreq');
%     %% CPCA- PLV
%     Params.itc=0;
%     [CPCAPLV,~]=FnMultiTaperFreqCPCA(EEG,Params);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     %     plot(SelFreqs,CPCAPLV,...
%     %         'Color',[0.8 0.8 0.8],...
%     %         'LineWidth',0.5);
%     plot(SelFreqs,mean(CPCAPLV,1),...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('Complex PCA-PLV')
%     xlim([2 40])
%     %% CPCA- ITC
%     Params.itc=1;
%     [CPCAITC,~]=FnMultiTaperFreqCPCA(EEG,Params);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     %     plot(SelFreqs,CPCAPLV,...
%     %         'Color',[0.8 0.8 0.8],...
%     %         'LineWidth',0.5);
%     plot(SelFreqs,mean(CPCAITC,1),...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('Complex PCA-PLV')
%     xlim([2 40]);
%     %% Time domian PCA
%     [YPC,EVec]=FnTimeDomainPCA(EEG);
%     figure,plot(YPC);
%     title(['YPC']);
%     xlabel('Sample No.');
%     ylabel('YPC');
%     Params.itc=0;
%     %% MT-Phase
%     [PhaseTapFind,SelFreqs]=FnMultiTaperFreqPhase(EEG,Params);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     plot(SelFreqs,PhaseTapFind,...
%         'Color',[0.8 0.8 0.8],...
%         'LineWidth',0.5);
%     plot(SelFreqs,mean(PhaseTapFind,1),...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('MT-Phasew')
%     xlim([2 40]);
% 
%     %% MT-Power
%     [PowerTapFIND,SelFreqs]=FnMultiTaperPower(EEG,Params);
%     % Convert to PSD-like scaling
%      PowerTapFIND = PowerTapFIND% / (Fs * 1);
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     plot(SelFreqs,10*log10(PowerTapFIND),...
%         'Color',[0.8 0.8 0.8],...
%         'LineWidth',0.5);
%     plot(SelFreqs,mean(10*log10(PowerTapFIND),1),...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('MT-Power (dB)')
%     xlim([1 40]);
%     %% Chronux
%     S=[];
    ParamsCh = struct( ...
        'Fs', Fs, ...
        'tapers', [2 3], ...
        'fpass', [1 40], ...
        'pad', 2, ...
        'trialave', 1 ...
        );
%     for Ch=1:size(EEG,1)
%         Data=squeeze(EEG(Ch,:,:));   % time*trial
%         [S(Ch,:),freqChron]=mtspectrumc(Data,ParamsCh);
%     end
%     figure;
%     plot(freqChron,10*log10(S'), 'LineWidth', 2);
%     xlabel('Frequency (Hz)');
%     ylabel('Power');
%     xlim([1 40]);
% 
% 
% 
%     figure('Color','w','Position',[100 100 700 500]);
%     hold on
%     plot(freqChron,10*log10(S'),...
%         'Color',[0.8 0.8 0.8],...
%         'LineWidth',0.5);
%     plot(freqChron,mean(10*log10(S),1)',...
%         'k','LineWidth',3);
%     xlabel('Frequency (Hz)')
%     ylabel('Chronux Power')
%     xlim([1 40])
%     %%
%     win=2.2 ;   % full epoch length = 550 sample
%     EEGPerm=permute(EEG,[2,1,3]);
%     [Sc, Cmat, Ctot, Cvec, Cent, fr] = CrossSpecMatc(EEGPerm, win,ParamsCh);
%     %  Sc	cross-spectral matrix (channels × channels × freq)
%     % Cmat	coherence matrix (freq-resolved)
%     % Ctot	total coherence
%     % Cvec	vectorized coherence (for stats)
%     % Cent	centered coherence
%     % f	frequency vector
%     % Assumes: f, Ctot, Cvec, Cent already exist
% 
%     figure('Color','w','Position',[100 100 900 600]);
% 
%     tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% 
%     % Optional: colorblind-friendly paletten
%     c1=[0 0.4470 0.7410];
%     c2=[0.8500 0.3250 0.0980];
%     c3=[0.4660 0.6740 0.1880];
% 
%     % ---------- Ctot ----------
%     nexttile
%     plot(fr, Ctot, 'LineWidth', 2.5, 'Color', c1);
%     grid on
%     box off
%     ylabel('C_{tot}')
%     set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % % ---------- Cvec ----------
%     % nexttile
%     % plot(nanmean(abs(Cvec),1), 'LineWidth', 2.5, 'Color', c2);
%     % grid on
%     % box off
%     % ylabel('C_{vec}')
%     % set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % ---------- Cent ----------
%     nexttile
%     plot(fr, Cent, 'LineWidth', 2.5, 'Color', c3);
%     grid on
%     box off
%     ylabel('C_{ent}')
%     xlabel('Frequency (Hz)')
%     set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % Global title
%     sgtitle('Cross-Spectral Coherence Measures','FontSize',14,'FontWeight','bold')
% 
% 
%     % long(f,:)==1
%     % [Sc, Cmat, Ctot, Cvec, Cent, f] = CrossSpecMatc(EEGPerm, win, params);
% 
%     [~,idx] = min(abs(fr-33));
%     figure,topoplot(abs(Cvec(idx,:)),OUTEEG.chanlocs);
% 
%     EEGPermLong=permute(EEG(:,:,find(Long(i,:)==1)),[2,1,3]);
%     [Sc, Cmat, Ctot, Cvec, Cent, fr]=CrossSpecMatc(EEGPermLong,win,ParamsCh);
%     %  Sc	cross-spectral matrix (channels × channels × freq)
% 
%     figure('Color','w','Position',[100 100 900 600]);
% 
%     tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% 
%     % Optional: colorblind-friendly paletteß
%     c1 = [0 0.4470 0.7410];
%     c2 = [0.8500 0.3250 0.0980];
%     c3 = [0.4660 0.6740 0.1880];
% 
%     % ---------- Ctot ----------
%     nexttile
%     plot(fr, Ctot, 'LineWidth', 2.5, 'Color', c1);
%     grid on
%     box off
%     ylabel('C_{tot}')
%     set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % % ---------- Cvec ----------
%     % nexttile
%     % plot(nanmean(abs(Cvec),1), 'LineWidth', 2.5, 'Color', c2);
%     % grid on
%     % box off
%     % ylabel('C_{vec}')
%     % set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % ---------- Cent ----------
%     nexttile
%     plot(fr, Cent, 'LineWidth', 2.5, 'Color', c3);
%     grid on
%     box off
%     ylabel('C_{ent}')
%     xlabel('Frequency (Hz)')
%     set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % Global title
%     sgtitle('Cross-Spectral Coherence Measures','FontSize',14,'FontWeight','bold')
% 
% 
% 
% 
%     EEGPermShort=permute(EEG(:,:,find(Long(i,:)==0)),[2,1,3]);
%     [Sc, Cmat, Ctot, Cvec, Cent, fr]=CrossSpecMatc(EEGPermShort,win,ParamsCh);
%     %  Sc	cross-spectral matrix (channels × channels × freq)
% 
% 
%     figure('Color','w','Position',[100 100 900 600]);
% 
%     tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% 
%     % Optional: colorblind-friendly paletteß
%     c1 = [0 0.4470 0.7410];
%     c2 = [0.8500 0.3250 0.0980];
%     c3 = [0.4660 0.6740 0.1880];
% 
%     % ---------- Ctot ----------
%     nexttile
%     plot(fr, Ctot, 'LineWidth', 2.5, 'Color', c1);
%     grid on
%     box off
%     ylabel('C_{tot}')
%     set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % % ---------- Cvec ----------
%     % nexttile
%     % plot(nanmean(abs(Cvec),1), 'LineWidth', 2.5, 'Color', c2);
%     % grid on
%     % box off
%     % ylabel('C_{vec}')
%     % set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % ---------- Cent ----------
%     nexttile
%     plot(fr, Cent, 'LineWidth', 2.5, 'Color', c3);
%     grid on
%     box off
%     ylabel('C_{ent}')
%     xlabel('Frequency (Hz)')
%     set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')
% 
%     % Global title
%     sgtitle('Cross-Spectral Coherence Measures','FontSize',14,'FontWeight','bold')




win=2.2;

% ---------- All trials ----------
EEGPerm=permute(EEG,[2,1,3]);
[~,~,CtotAll,CvecAll,CentAll,fr]=CrossSpecMatc(EEGPerm,win,ParamsCh);

% ---------- Long trials ----------
EEGPermLong=permute(EEG(:,:,SubjectiveLongBin(i,:)==1),[2,1,3]);
[~,~,CtotLong,CvecLong,CentLong,~]=CrossSpecMatc(EEGPermLong,win,ParamsCh);

% ---------- Short trials ----------
EEGPermShort=permute(EEG(:,:,SubjectiveShortBin(i,:)==1),[2,1,3]);
[~,~,CtotShort,CvecShort,CentShort,~]=CrossSpecMatc(EEGPermShort,win,ParamsCh);

% ---------- Plot ----------
figure('Color','w','Position',[100 100 900 700]);

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

cAll=[0 0 0];
cLong=[0 0.4470 0.7410];
cShort=[0.8500 0.3250 0.0980];

% ---------- Total coherence ----------
nexttile
plot(fr,CtotAll,'k','LineWidth',2); hold on
plot(fr,CtotLong,'Color',cLong,'LineWidth',2);
plot(fr,CtotShort,'Color',cShort,'LineWidth',2);

ylabel('C_{tot}')
title('Total Coherence')
legend({'All','Long','Short'})
grid on
box off
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

% ---------- Entropy coherence ----------
nexttile
plot(fr,CentAll,'k','LineWidth',2); hold on
plot(fr,CentLong,'Color',cLong,'LineWidth',2);
plot(fr,CentShort,'Color',cShort,'LineWidth',2);

ylabel('C_{ent}')
xlabel('Frequency (Hz)')
title('Coherence Entropy')
legend({'All','Subj. Long','Subj. Short'})
grid on
box off
set(gca,'FontSize',12,'LineWidth',1.2,'TickDir','out')

sgtitle('Cross-Spectral Coherence: All vs Long vs Short','FontSize',14,'FontWeight','bold')










end