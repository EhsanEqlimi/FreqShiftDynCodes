clc;clear;close all
%
% !ssh -T git@github.com
% system('git clone https://github.com//EhsanEqlimi/ComplexPCA-PD.git');
% addpath(genpath(fullfile(pwd,'ComplexPCA-PD')));
addpath '/Users/ehsaneqlimi/eeglab2026.0.0'
eeglab;close

% Add fieldtrip toolbox
addpath('/Users/ehsaneqlimi/fieldtrip-20260518');
ft_defaults; % -> really important!


PreProcseedDataPath='/Users/ehsaneqlimi/FreqShiftDynData/preprocessed_EEG/';

Dir=dir(fullfile(PreProcseedDataPath,'*.mat'));
%epoch relative to the start of the delay (stimulus offset)
delaystart = {'S101','S102','S103', 'S104', 'S105', 'S106','S107','S108',...
    'S201','S202','S203', 'S204', 'S205', 'S206','S207','S208',...
    'S 45','S 46','S 47', 'S 48', 'S 49', 'S 50','S 51','S 52'};

test_durations = [0.2 0.250 0.319 0.331 0.369 0.381 0.450 0.5;...
    0.45 0.5 0.619 0.669 0.706 0.756 0.870 0.920;...
    0.870 0.920 0.981 1.169 1.231 1.419 1.470 1.520];

EEGdata = nan(length(Dir), 96,  550,  336);
bad_trials = nan(length(Dir),336);
block_type= nan(length(Dir),336);
duration= nan(length(Dir),336);
long= nan(length(Dir),336);
correct= nan(length(Dir),336);
rt= nan(length(Dir),336);


% Configure parameters for time-frequency analysis
Params=struct(...
    'Fs', 250, ...              % Sampling frequency (Hz)
    'tapers', [2, 3], ...       % Time-frequency tapers (TW, K)
    'Fpass', [12, 36], ...       % Frequency range of interest (Hz)
    'itc', 0, ...               % Flag for inter-trial coherence calculation
    'pad', 1, ...               % Padding factor for FFT
    'SegmentSec', 2);           % Length of each time segment (seconds)
ParamsChronux=struct(...
    'Fs', 250, ...              % Sampling frequency (Hz)
    'tapers', [2, 3], ...       % Time-frequency tapers (TW, K)
    'fpass', [12, 36],'SegmentSec', 2);      % Frequency range of interest (Hz)

for f=1:length(Dir)

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




    % create a channel x trial x time tensor
%     EEGDat3=permute(EEG_epoched.data,[ 1 3 2]);
%     [YCPC,CWTS,CWTSfIND,SelFreqs,CPCFreqfIND,MeanOscilatoryPowerfInd]=FnMutiTaperGlobalCoherence(EEGDat3,Params);

   % chronux
    %create samples x channels x trials data
    EEGDat3=permute(EEG_epoched.data,[2 1 3]);

    %create samples x trial x channels data
    EEGDat3_trial=permute(EEG_epoched.data,[2 3 1]);
  
   [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(EEGDat3_trial,ParamsChronux.SegmentSec,ParamsChronux);


      % create a channel x trial x time tensor
    EEGDat3_EFR=permute(EEG_epoched.data,[ 1 3 2]);
   [PLVCPCAFIND,SelFreqs]=FnMultiTaperFreqCPCA(EEGDat3_EFR,Params);

    Params.itc=1;
    [CPCAPLV,~]=FnMultiTaperFreqCPCA(EEGDat3_EFR,Params);
    figure,plot(SelFreqs,CPCAPLV);
    title(['CPCAPLV-' Name(1:end-4)]);
    xlabel('Freq (Hz)');
    ylabel('CPCAPLV');
   

    Params.itc=0;
        EEGDat3_EFR=EEG_epoched.data;

    [YCPC,CWTS,CWTSfIND,SelFreqs,CPCFreqfIND]=FnMTCPCATIMEDomain(EEGDat3_EFR,Params);
        Indf100=find(SelFreqs==20);

    figure,plot(YCPC);
   
    xlabel('Sample No.');
    ylabel('YCPC');

    
    figure,plot(SelFreqs,abs(CPCFreqfIND));
    xlabel('Freq (Hz)');
    ylabel('CPCFreq');
    
    figure,topoplot( abs(CWTSfIND(:,2)), OUTEEG.chanlocs,'electrodes','numbers')


%%
    Fs=250;
    [Z,D,beta_variability,Smat]=FnBetaTrialCPCA(EEG_epoched.data,Fs);
    

%%
    figure;
    if size(Z,2)>=3
        scatter3(real(Z(:,1)),real(Z(:,2)),real(Z(:,3)),60,1:size(Z,1),'filled');
        xlabel('CPCA1'); ylabel('CPCA2'); zlabel('CPCA3');
    else
        scatter(real(Z(:,1)),imag(Z(:,1)),60,1:size(Z,1),'filled');
        xlabel('Real CP1'); ylabel('Imag CP1');
    end
    title('Beta CPCA Trial State Space');
    colorbar;

    %%%

    figure;

%% =========================================================
% 1. CPCA STATE SPACE
%% =========================================================
subplot(2,2,1);

if size(Z,2)>=3
    scatter3(real(Z(:,1)),real(Z(:,2)),real(Z(:,3)),60,1:size(Z,1),'filled');
    xlabel('CPCA1'); ylabel('CPCA2'); zlabel('CPCA3');
else
    scatter(real(Z(:,1)),imag(Z(:,1)),60,1:size(Z,1),'filled');
    xlabel('Real CP1'); ylabel('Imag CP1');
end

title('Beta CPCA State Space');
colorbar;


%% =========================================================
% 2. DISTANCE MATRIX (VARIABILITY GEOMETRY)
%% =========================================================
subplot(2,2,2);

imagesc(D);
axis square;
colormap hot;
colorbar;

xlabel('Trial');
ylabel('Trial');
title('Inter-Trial Beta Distance');


%% =========================================================
% 3. SIMILARITY MATRIX
%% =========================================================
subplot(2,2,3);

imagesc(Smat);
axis square;
colormap parula;
colorbar;

xlabel('Trial');
ylabel('Trial');
title('Beta State Similarity');


%% =========================================================
% 4. GLOBAL VARIABILITY INDEX
%% =========================================================
subplot(2,2,4);

bar(beta_variability);
ylabel('Mean Distance');
title('Global Beta Variability');

sgtitle('Beta Oscillatory State Space Analysis');

end

