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
    0.45 0.5 0.619 0.669 0.706 0.756 0.870 0.920;...
    0.870 0.920 0.981 1.169 1.231 1.419 1.470 1.520];
% EEGdat--> channel x time x trial
EEGdata = nan(length(Dir), 96,  550,  336);
bad_trials = nan(length(Dir),336);
block_type= nan(length(Dir),336);
duration= nan(length(Dir),336);
long= nan(length(Dir),336);
correct= nan(length(Dir),336);
rt= nan(length(Dir),336);

%% Test for Emina
% load([PreProcseedDataPath '10_C10M2_ft.mat']);
% Data_Test=[];
% PSD_mean=[];
% for t=1:length(data_ft.trial)
%     Data_Test(:,:,t)=data_ft.trial{t};
% end
% x_mean = mean(Data_Test,3);   % channel × time
% fs=data_ft.fsample;
% 
% for ch = 1:size(x_mean,1)
%     [PSD_mean(:,ch), F] = pwelch(x_mean(ch,:), [], [], [], fs);
% end
% % Keep only 2-40 Hz
% idx = F >= 2 & F <= 40;
% 
% F = F(idx);
% PSD_mean = PSD_mean(idx,:);
% figure,plot(F,PSD_mean)
% figure,plot(SpecSel)
% 
% % 1/f effect (simple log-log)
% %Note:% Computed PSD for each epoch, average across epochs (via spectopo), then fit/remove 1/f
% % (preferred over fitting 1/f separately for each epoch)
% %     FitIdx=Freq >= 2 & Freq <= 40;
% logF=log10(F);
% SpecFlat=[];
% SpecSel=log10(PSD_mean); % convert dB → log10(power)
% for Ch=1:size(PSD_mean,2) %channel loop
%     ChannelPower=SpecSel(:,Ch); %no need to log10 (spectopo is db)
%     Powfitted=polyfit(logF,ChannelPower,1);
%     AperiodicFit = polyval(Powfitted,logF);
%     SpecFlat(:,Ch) = ChannelPower-AperiodicFit;
%     Alpha(Ch)=-Powfitted(1);
% 
% end
% 
% %plot falttend PSD
% figure;
% plot(F,SpecFlat);
% xlabel('Frequency (Hz)');
% ylabel('Power (dB, corrected)');
% title('Flattened Spectrum');
% xlim([2 40]);
% 
% 
% %plot falttend PSD
% figure;
% plot(FreqSel,SpecFlat');
% xlabel('Frequency (Hz)');
% ylabel('Power (dB, corrected)');
% title('Flattened Spectrum');
% xlim([2 40]);
%% End test for Emina (mice dta and 1/f effect, 9 June 2026)


for f=4:length(Dir)

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

    %% 10 June Idea
    %% =========================
    %  INTER-TRIAL BETA VARIABILITY PIPELINE
    %  NO GED, NO HMM, NO GA-CCA
    %% =========================

  

    %% -------------------------
    % INPUT DATA
    % X: trials x channels x time
    %% -------------------------

    % Example placeholders (replace with your data)
    % X = randn(nTrials, nCh, nTime);
    X=double(permute(EEG_epoched.data,[3,1,2]));

    [nTrials, nCh, nTime] = size(X);
    Fs = 250;

    %% -------------------------
    % STEP 1: BANDPASS FILTER (BETA)
    %% -------------------------

    bp = [13 30];
    [b,a] = butter(4, bp/(Fs/2), 'bandpass');

    Xf = zeros(size(X));

    for i = 1:nTrials
        for ch = 1:nCh
            Xf(i,ch,:) = filtfilt(b,a,squeeze(X(i,ch,:)));
        end
    end

    %% -------------------------
    % STEP 2: PCA SPATIAL REDUCTION
    % (NO GED, NO TEMPLATE)
    %% -------------------------
% Assumes data is loaded as a 3D matrix: [channels x time x trials] 
% For example: [nChans, nTimes, nTrials] = size(EEG_Data);

% Step 1: Reshape data to 2D (Channels x all time points/trials)
% Concatenate trials/time to preserve the spatial characteristics
Xf=permute(X,[2,3,1]);
[nChans, nTimes, nTrials] = size(Xf);
data2D = reshape(Xf, nChans, nTimes * nTrials);

% Step 2: Zero-center the data across channels
mean_data = mean(data2D, 2);
data_centered = data2D - mean_data;

% Step 3: Compute the spatial covariance matrix
% The covariance of the channels
covMatrix = (1 / (size(data_centered, 2) - 1)) * data_centered * data_centered';

% Step 4: Perform Eigenvalue Decomposition
[V, D] = eig(covMatrix);

% Step 5: Sort eigenvalues in descending order and flip eigenvectors accordingly
[eigenvalues, sortIdx] = sort(diag(D), 'descend');
V_sorted = V(:, sortIdx);

% Step 6: Project data (Spatial Filters)
% Create spatial filters using the top eigenvectors (e.g., top 5)
num_components = 5;
spatial_filters = V_sorted(:, 1:num_components);

% Project original data into the new PCA space
pca_time_series = spatial_filters' * data2D;

% Reshape back into components x time x trials to evaluate the ERP
pca_components = reshape(pca_time_series, num_components, nTimes, nTrials);

% (Optional) Scree plot to visualize variance explained
variance_explained = (eigenvalues / sum(eigenvalues)) * 100;
figure;
plot(variance_explained(1:20), 'ko-', 'LineWidth', 2);
title('Scree Plot of Spatial PCA Components');
xlabel('Component Number');
ylabel('Percentage of Variance Explained');
disp(size(pca_components))  % should be nTrials × K × nTime
Xpca=pca_components;
K=num_components;
    %% -------------------------
    % STEP 3: SLIDING WINDOWS
    %% -------------------------

    win  = round(0.4 * Fs);
    step = round(0.05 * Fs);

    tIdx = 1:step:(nTime-win);
    nWin = length(tIdx);

    %% -------------------------
    % STEP 4: FREQUENCY ESTIMATION (MUSIC-LIKE PEAK)
    %% -------------------------

    fgrid = 13:0.1:30;

    f_est = zeros(nTrials, nWin);

    for i = 1:nTrials

        Xi = squeeze(Xpca(i,:,:)); % K x time

        for w = 1:nWin

            seg = Xi(:, tIdx(w):tIdx(w)+win-1);

            % covariance
            R = (seg * seg') / size(seg,2);

            % eigen decomposition
            [E,D] = eig(R);
            [~,idx] = sort(diag(D),'descend');
            E = E(:,idx);

            En = E(:,2:end); % noise subspace (1 dominant source assumption)

            P = zeros(length(fgrid),1);

            for k = 1:length(fgrid)

                f = fgrid(k);

                % steering vector (simplified multichannel surrogate)
                t = (0:size(seg,2)-1)/Fs;
                a = exp(-1j*2*pi*f*t);
                a = repmat(a, K, 1);
                a = a(:);

                P(k) = 1 / abs(a' * (En*En') * a);
            end

            [~,ix] = max(P);
            f_est(i,w) = fgrid(ix);

        end
    end

    %% -------------------------
    % STEP 5: BUILD DISTRIBUTION p(f|t)
    %% -------------------------

    fgrid_pdf = 13:0.1:30;
    p = zeros(length(fgrid_pdf), nWin);

    for w = 1:nWin

        data = f_est(:,w);

        p(:,w) = ksdensity(data, fgrid_pdf);

    end

    %% normalize
    for w = 1:nWin
        p(:,w) = p(:,w) / sum(p(:,w));
    end

    %% -------------------------
    % STEP 6: VARIABILITY METRICS
    %% -------------------------

    mu_f  = mean(f_est,1);
    var_f = var(f_est,0,1);

    H = zeros(1,nWin);
    KL = zeros(1,nWin);

    p0 = p(:,1) + eps;
    p0 = p0 / sum(p0);

    for w = 1:nWin

        pw = p(:,w) + eps;
        pw = pw / sum(pw);

        % entropy
        H(w) = -sum(pw .* log(pw));

        % KL divergence vs baseline
        KL(w) = sum(pw .* log(pw ./ p0));

    end

    %% -------------------------
    % STEP 7: VISUALIZATION
    %% -------------------------

    time = linspace(-2,0,nWin);

    figure;
    plot(time, mu_f, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('Mean Beta Frequency Across Trials');

    figure;
    plot(time, var_f, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Variance');
    title('Inter-Trial Frequency Variability');

    figure;
    plot(time, H, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Entropy');
    title('Distribution Entropy');

    figure;
    plot(time, KL, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('KL Divergence vs Baseline');
    title('Distribution Shift Over Time');

    figure;
    imagesc(time, fgrid_pdf, p);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('p(f | t) Across Trials');
    colorbar;


    %%
    %%
    %%=========================================================
    % BETA FILTER
    %%=========================================================
    %EEG -> [NTrials x NChannels x NTimes]
    EEGDat=permute(EEG_epoched.data,[3,1,2]);
    [NTrials,NChannels,NTimes]=size(EEGDat);
    Fs=EEG_epoched.srate;
    BetaBand=[13,36];
    NComponents=3;
    [B,A]=butter(4,BetaBand/(EEG_epoched.srate/2),'bandpass');
    NTrails=size(EEG_epoched,3);
    EEGBeta=zeros(size(EEGDat));

    for ITrial=1:NTrials

        Trial=double(squeeze(EEGDat(ITrial,:,:)));

        EEGBeta(ITrial,:,:)=filtfilt(B,A,Trial')';

    end

    %%=========================================================
    % SINGLE-TRIAL COVARIANCE
    %%=========================================================

    Covs=zeros(NChannels,NChannels,NTrials);

    for ITrial=1:NTrials

        X=squeeze(EEGBeta(ITrial,:,:));

        C=cov(X');

        % regularization
        C=C+1e-6*eye(size(C));

        Covs(:,:,ITrial)=C;

    end


    %%=========================================================
    % SUBSPACE EXTRACTION
    %%=========================================================

    Subspaces=cell(NTrials,1);

    for ITrial=1:NTrials

        C=Covs(:,:,ITrial);

        [V,D]=eig(C);

        [~,Idx]=sort(diag(D),'descend');

        V=V(:,Idx);

        U=V(:,1:NComponents);

        Subspaces{ITrial}=U;

    end


    %%=========================================================
    % GRASSMANN DISTANCE MATRIX
    %%=========================================================

    DGrass=zeros(NTrials);

    for ITrial=1:NTrials

        U1=Subspaces{ITrial};

        for JTrial=1:NTrials

            U2=Subspaces{JTrial};

            Angle=subspace(U1,U2);

            DGrass(ITrial,JTrial)=Angle;

        end
    end


    %%=========================================================
    % RIEMANNIAN DISTANCE MATRIX
    %%=========================================================

    DRiem=zeros(NTrials);

    for ITrial=1:NTrials

        C1=Covs(:,:,ITrial);

        for JTrial=1:NTrials

            C2=Covs(:,:,JTrial);

            EigVals=eig(C1\C2);

            DRiem(ITrial,JTrial)=sqrt(sum(log(EigVals).^2));

        end
    end


    %%=========================================================
    % MANIFOLD VISUALIZATION
    %%=========================================================
    % enforce symmetry
    DGrass=(DGrass+DGrass')/2;

    % enforce diagonal = 0
    for ITrial=1:NTrials
        DGrass(ITrial,ITrial)=0;
    end

    % remove NaN / Inf
    DGrass(~isfinite(DGrass))=0;

    % ensure non-negative
    DGrass(DGrass<0)=0;

    % optional: enforce double precision
    DGrass=double(DGrass);
    Y=mdscale(DGrass,2);

    figure;

    scatter(Y(:,1),Y(:,2),80,'filled');

    xlabel('Dim 1');
    ylabel('Dim 2');

    title('Beta Subspace Geometry');

    grid on;


    %%=========================================================
    % CLUSTERING
    %%=========================================================

    NClusters=3;

    ClusterIdx=kmeans(Y,NClusters);

    figure;

    gscatter(Y(:,1),Y(:,2),ClusterIdx);

    xlabel('Dim 1');
    ylabel('Dim 2');

    title('Beta State Clusters');

    grid on;


    %%=========================================================
    % INTER-TRIAL VARIABILITY
    %%=========================================================

    Variability=mean(DGrass,2);

    figure;

    plot(Variability,'LineWidth',2);

    xlabel('Trial');

    ylabel('Variability');

    title('Inter-Trial Beta Variability');

    grid on;


    %%=========================================================
    % ROOT-MUSIC PER TRIAL
    %%=========================================================

    BetaPeaks=zeros(NTrials,1);

    for ITrial=1:NTrials

        X=squeeze(EEGBeta(ITrial,3,:));

        [S,F]=pmusic(X,4,[],Fs);

        BetaIdx=find(F>=15 & F<=30);

        [~,MaxIdx]=max(S(BetaIdx));

        BetaPeaks(ITrial)=F(BetaIdx(MaxIdx));

    end


    figure;

    histogram(BetaPeaks);

    xlabel('Beta Peak Frequency');

    ylabel('Count');

    title('Trial-wise Beta Peaks');

    %%CCA

    EEGdata(f,:,:,1:size(EEG_epoched.data,3)) =  EEG_epoched.data;
    % apply bad trial indices
    NumChan=size(EEG_epoched.data,1);
    NumTime=size(EEG_epoched.data,2);
    NumGoodTrial=numel(bad_trials(f,:)==0);
    EpochedCatX=reshape(EEG_epoched.data,NumChan,NumTime*NumGoodTrial)'; %this is X matrix channel x time*trial (for CCA)
    GACatY=repmat(nanmean(EEG_epoched.data(:,:,bad_trials(f,:)==0),3),[1,NumGoodTrial])'; %this is Y matrix channel x time*trial (for CCA)
    [A,B,r,U,V]=canoncorr(GACatY,EpochedCatX); %X, Y are usually are not full-rank
    %n_cca_all = size(R_all, 2);
    NumComp=6;length(r);
    CCACat=(EpochedCatX*B(:,1:NumComp))';

    % Spatial patterns
    SpatialPattern=cov(EpochedCatX)*B;



    %% Re-reshape
    CCAComps=reshape(CCACat, NumComp,NumTime,NumGoodTrial);
    x_lim = [-200 2000];
    pt_plot = dsearchn(EEG_epoched.times', x_lim'); % for getting the scaling right in interval of interest
    y_lim = [-6 6];
    hh = figure; set(gcf, 'Position', [678,233,907,647]);

    for k = 1:NumComp
        subplot(NumComp,3,k*3-2)
        plot(EEG_epoched.times(pt_plot(1):pt_plot(2)), mean(CCAComps(k,pt_plot(1):pt_plot(2),:),3))
        xlim(x_lim)
        %ylim(y_lim)

        subplot(NumComp,3,k*3-1)
        imagesc(EEG_epoched.times(pt_plot(1):pt_plot(2)), linspace(1, size(CCAComps, 3), size(CCAComps, 3)),...
            squeeze(CCAComps(k,pt_plot(1):pt_plot(2),:))');
        colorbar; xlabel('time in ms'); ylabel('trials'); axis('xy')
        %xlim(x_lim)
        %caxis([-5 5])
        caxis([-2 2])

        subplot(NumComp,3,k*3)
        topoplot(SpatialPattern(:,k), EEG_epoched.chanlocs);
        colorbar
    end
    figure,
    ShortIdx=find(bad_trials(f,:)==0 & long(f,:)==0);%1 ;
    LongIdx=find(bad_trials(f,:)==0 & long(f,:)==1);%1 ;
    for k=1:NumComp
        subplot(NumComp,3,k)

        plot(EEG_epoched.times(pt_plot(1):pt_plot(2)),...
            mean(CCAComps(k,pt_plot(1):pt_plot(2),ShortIdx),3),...
            'LineWidth',2)

        hold on

        plot(EEG_epoched.times(pt_plot(1):pt_plot(2)),...
            mean(CCAComps(k,pt_plot(1):pt_plot(2),LongIdx),3),...
            'LineWidth',2)

        xlim(x_lim)

        if k==1
            legend('Short','Long')
        end

        title(['CCA Comp ' num2str(k)])
    end

    %% LDA calssifier
    Labels=bad_trials(f,:)==0 & long(f,:)==1;%1 long zero short
    CompIdx=1:3;

    Win1=EEG_epoched.times>=150 & EEG_epoched.times<=300;
    Win2=EEG_epoched.times>=300 & EEG_epoched.times<=600;
    Win3=EEG_epoched.times>=600 & EEG_epoched.times<=1000;

    NumTrial=size(CCAComps,3);

    Features=[];

    for Tr=1:NumTrial

        FeatureVec=[];

        for Comp=CompIdx

            FeatureVec=[FeatureVec ...
                mean(CCAComps(Comp,Win1,Tr)) ...
                mean(CCAComps(Comp,Win2,Tr)) ...
                mean(CCAComps(Comp,Win3,Tr))];

        end

        Features(Tr,:)=FeatureVec;

    end

    Mdl=fitcdiscr(Features,Labels);

    CVMdl=crossval(Mdl,'KFold',10);

    Accuracy=1-kfoldLoss(CVMdl);

    fprintf('Accuracy=%.2f%%\n',100*Accuracy);
    [~,Score]=resubPredict(Mdl);

    LDAScore=Score(:,2);
    figure

    subplot(1,2,1)

    histogram(LDAScore(Labels==0))
    hold on
    histogram(LDAScore(Labels==1))

    xlabel('LDA score')
    ylabel('Count')
    legend('Short','Long')

    subplot(1,2,2)

    boxplot(LDAScore,Labels)

    ylabel('LDA score')
    xticklabels({'Short','Long'})

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