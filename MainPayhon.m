clc;clear;close all
%i create an clean eeg env in bash terminal before
system('ls /Users/ehsaneqlimi/miniconda3/envs/eeg/bin/')
pyenv('Version','/Users/ehsaneqlimi/miniconda3/envs/eeg/bin/python')
%uncomment if yoou do not want to use eeg env below
% pyenv('Version','/usr/bin/python3')
% pyenv
% i have already installed python3, pip3 and conda, miniconda3 in bash
% now i wan to install numpy, scipy and mne via conda in matlab
% system('~/miniconda3/bin/conda --version')
% system('~/miniconda3/bin/conda info --base')
% %Conda versions require you to accept the Anaconda Terms of Service before using the default channels
% system('~/miniconda3/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main')
% system('~/miniconda3/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r')
% 
% 
% system('~/miniconda3/bin/conda install numpy scipy -y') % quite fast
% system('~/miniconda3/bin/conda install -c conda-forge mne -y') %take a while
% system('~/miniconda3/bin/python -m pip install matplotlib')
%%
%Test Python inside MATLAB
py.importlib.import_module('sys');
% mne
system('~/miniconda3/bin/python -c "import mne; print(mne.__version__)"')
% read data
DataPath = '/Users/ehsaneqlimi/FreqShiftDynData/preprocessed_EEG/';
Files=dir(fullfile(DataPath, '*.mat'));
%% these lines did not work
% ne = py.importlib.import_module('mne');
% scipy = py.importlib.import_module('scipy.io');
% np = py.importlib.import_module('numpy');
% so i create a clen eeg env using bash: conda create -n eeg python=3.11 numpy scipy mne matplotlib -c conda-forge 
% my matlab does not support python 3.11 so again
%conda create -n eeg python=3.8 numpy scipy mne matplotlib -c conda-forge
% now test eeg env (first restart matlab)
for i = 1:length(files)

    fname = fullfile(dataPath, files(i).name);

    % ---- load MATLAB EEG struct in Python ----
    mat = scipy.loadmat(fname, pyargs('struct_as_record', false, 'squeeze_me', true));
    EEG = mat{'EEG'};

    disp(['Processing: ' files(i).name])

    % =====================================================
    % Convert EEGLAB → MNE Raw (manual minimal conversion)
    % =====================================================

    data = np.array(EEG.data);          % EEG signals
    sfreq = double(EEG.srate);          % sampling rate

    % channel names
    nCh = length(EEG.chanlocs);
    ch_names = py.list();

    for c = 1:nCh
        ch_names.append(string(EEG.chanlocs(c).labels));
    end

    info = mne.create_info(ch_names, sfreq, pyargs('ch_types','eeg'));

    raw = mne.io.RawArray(data, info);

    % =====================================================
    % BASIC MNE PROCESSING
    % =====================================================

    raw.filter(1, 40);     % bandpass filter
    raw.set_eeg_reference('average');

    % save processed version
    outName = fullfile(dataPath, ['mne_' files(i).name '.fif']);
    raw.save(outName, pyargs('overwrite', true));

end
