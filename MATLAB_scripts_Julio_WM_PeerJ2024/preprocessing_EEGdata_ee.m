restoredefaultpath
% addpath '/Users/juliorodriguezlarios/Documents/MATLAB/eeglab2022.0'
addpath '/Users/ehsaneqlimi/eeglab2026.0.0'
eeglab;close
% plugin_askinstall('BVA-io', [], 1);
% make sure thath you ruin the following command in terminal by ee
% open  smb://omh-smb.svc.ny.gov/OMH_Shared/NYSPI/HaegensLab

% clone Jones' burst detection repo and addapth
system('git clone https://github.com/jonescompneurolab/SpectralEvents.git');
addpath(genpath(fullfile(pwd,'SpectralEvents')))

% EEGpath = '/Users/juliorodriguezlarios/Desktop/Projects/categorization EEG/categorization_task_EEG';
EEGpath = '/Volumes/HaegensLab/DATA/project_categorization/EEG/';

% behaviorpath = '/Users/juliorodriguezlarios/Desktop/Projects/categorization EEG/behavior_categorization';
behaviorpath='/Volumes/HaegensLab/DATA/project_categorization/EEG/'; %lab shared folder

%load locations of Cz
% load('/Users/juliorodriguezlarios/FreqeShiftDynCodes/cz_locs.mat')
TemplateEEGChanLoc=pop_chanedit([], 'load', '/Users/ehsaneqlimi/eeglab2026.0.0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');

CzIdx=find(strcmpi({TemplateEEGChanLoc.labels}, 'Cz'));
cz_locs = TemplateEEGChanLoc(CzIdx);
cz_locs = rmfield(cz_locs, 'urchan');
path2save = '/Users/ehsaneqlimi/FreqShiftDynData/preprocessed_EEG/';


%Loop subjects

test_durations = [0.2 0.250 0.319 0.331 0.369 0.381 0.450 0.5;...
    0.45 0.5 0.619 0.669 0.706 0.756 0.870 0.920;...
    0.870 0.920 0.981 1.169 1.231 1.419 1.470 1.520];


for s = 7:25
    %%%% LOAD DATA %%%%%%%%%%

    if or(or(s == 24 , s== 20), s==19)%subjects with 2 sessions
        %load EEG sessions and concatenate
        temp_name = ['S_', num2str(s), '.vhdr'];
        [EEG1, com] = pop_loadbv(EEGpath, temp_name);

        temp_name = ['S_', num2str(s),'_2', '.vhdr'];
        [EEG2, com] = pop_loadbv(EEGpath, temp_name);

        EEG = pop_mergeset(EEG1, EEG2);
    end


    if s == 23%subjects with 3 sessions
        %load EEG sessions and concatenate
        temp_name = ['S_', num2str(s), '.vhdr'];
        [EEG1, com] = pop_loadbv(EEGpath, temp_name);

        temp_name = ['S_', num2str(s),'_2', '.vhdr'];
        [EEG2, com] = pop_loadbv(EEGpath, temp_name);

        temp_name = ['S_', num2str(s),'_3', '.vhdr'];
        [EEG3, com] = pop_loadbv(EEGpath, temp_name);

        EEGtemp = pop_mergeset(EEG1, EEG2);
        EEG = pop_mergeset(EEGtemp, EEG3);

    end

    %for subjects with only 1 session
    temp = 1:25;
    temp([19 23 24 20 22]) = [];
    if logical(find(s == temp))
        temp_name = ['S_', num2str(s), '.vhdr'];
        [EEG, com] = pop_loadbv(EEGpath, temp_name);
    end


    %%%%%%%% PREPROCESSING %%%%%%%%%%


    %resample to make computations faster
    [EEG] = pop_resample( EEG, 250);

    %filter data
    EEG = pop_eegfiltnew(EEG, [], 0.5, [], true, [], 0); % Highpass filter
    EEG = pop_eegfiltnew(EEG, [], 50, [], false, [], 0); % Lowpass filter

    %remove bipolar channels for the rest of the analysis
    OUTEEG = pop_select(EEG, 'channel', 1:95);

    %add reference channel
    OUTEEG.data(96,:) = 0;
    OUTEEG.chanlocs(96) = cz_locs;
    OUTEEG.nbchan = 96;

    OUTEEG = pop_reref(OUTEEG,[]);%rereference to average

    %save original electrode locations (in case you need to interpolate)
    temp_locs = OUTEEG.chanlocs;

    %detect and delete flat channels
    OUTEEG = clean_flatlines(OUTEEG,60);

    %detect and delete noisy channels
    %correlation with its robust estimate
    OUTEEG = clean_channels(OUTEEG,0.5); %delete noisy channels

    %keep track of number of bad electrodes
    OUTEEG.int_electrodes = length(temp_locs) - size(OUTEEG.data,1);

    %use asr method for cleaning
    OUTEEG = clean_asr(OUTEEG,20);

    %interpolate flat or bad electrodes if any
    OUTEEG = pop_interp(OUTEEG,temp_locs, 'spherical');

    %run ICA adjusting for data rank
    %number of components = number of non-interpolated - 1 (-1 is because of the common reference)
    OUTEEG= pop_runica(OUTEEG, 'icatype', 'runica','options',{'pca', length(temp_locs)-OUTEEG.int_electrodes-1} );
    OUTEEG = iclabel(OUTEEG); %classify components
    %mark components that are muscle, eyes, heart or channel noise (p>80%)
    [component,type] = find(OUTEEG.etc.ic_classification.ICLabel.classifications(:,[2 3 4 6])>0.8);

    %save bipolar channels (HEOG, VEOG, ECG)
    OUTEEG.bipolar = EEG.data(96:98,:);

    %save components' timecourses
    OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);

    %mark components highly correlated to bipolar channels to delete
    [r, p] = corr(OUTEEG.icaact', OUTEEG.bipolar');
    [extra_comp, ~] = find(r>0.8);

    if ~isempty(extra_comp)
        component = [component;extra_comp];
    end

    if ~isempty(component) %if there are noisy components
        OUTEEG= pop_subcomp(OUTEEG, unique(component),0);
        OUTEEG.rejected_components = length(unique(component)); %keep track of components you reject
    end


    save(strcat(path2save, temp_name, '.mat'), 'OUTEEG')



end

%% Burst detection based on jones repo










