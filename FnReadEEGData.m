function EEGData=FnReadEEGData(EEGDir)

FName=fullfile(EEGDir); % file name (without extension)

EEGFile=strcat(FName,'.eeg'); % eeglab 'EEG' file (EEG variable structure)
% % % fdtFile = strcat(fname,'.vhdr'); % eeglab 'vhdr' file (EEG signals)

load('-mat',EEGFile); % reading eeglab 'set' file

fid = fopen(fdtFile,'r','ieee-le');
EEG.data = fread(fid,[EEG.nbchan,Inf],'float32'); % reading eeglab 'fdt' file

eeg.data = EEG.data'; % complete EEG measurements data
% % % % eeg.Fs = EEG.srate; % EEG dataset sampling frequency