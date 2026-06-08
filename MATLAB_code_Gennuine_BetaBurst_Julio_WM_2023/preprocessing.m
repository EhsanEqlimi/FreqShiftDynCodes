
%define EEGlab and data path
restoredefaultpath
addpath '/Users/juliorodriguezlarios/Documents/MATLAB/eeglab2022.0'
eeglab;close
path = '/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/Raw_EEG';

%path2save = '/Users/juliorodriguezlarios/Desktop/WM_EEG_study/preprocessed_EEG/';
path2save = '/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/preprocessed_EEG_betafiltered/';

%go to folder with all data
cd '/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/Raw_EEG'

%load locations of Cz
load('/Users/juliorodriguezlarios/Desktop/Projects/WM_EEG_study/variables/cz_locs.mat')

for s = 1:31 %loop subjects
    
    temp_name = ['S_', num2str(s), '.vhdr'];
    
    if ~(s == 13 || s == 11 || s == 6)%excluding noisy subjects

        %load data
        [EEG, com] = pop_loadbv(path, temp_name);
        
        %add reference channel
        EEG.data(96,:) = 0;
        EEG.chanlocs(96) = cz_locs;
        EEG.nbchan = 96;
        
        %rereference to average
        EEG = pop_reref(OUTEEG,[]);
        
        %filter data
        EEG = pop_eegfiltnew(EEG, [], 1, [], true, [], 0); % Highpass filter 
        EEG = pop_eegfiltnew(EEG, [], 50, [], false, [], 0); % Lowpass filter
        
        %remove bipolar channels for the rest of the analysis
        OUTEEG = pop_select(EEG, 'channel', 1:95);
        
        %save original electrode locations (in case you need to interpolate)
        temp_locs = OUTEEG.chanlocs;
        
        %detect and delete flat channels
        OUTEEG = clean_flatlines(OUTEEG,60);
        
        %detect and delete noisy channels 
        %correlation with its robust estimate 
        OUTEEG = clean_channels(OUTEEG,0.5); %delete noisy channels
     
        %in the first 3 subjects delete Fcz (it was misplaced)
        if s==1 || s==2|| s==3
            OUTEEG = pop_select(OUTEEG, 'nochannel', {'FCz'});
        end
        
        %keep track of number of interpolated electrodes
        OUTEEG.int_electrodes = 95 - size(OUTEEG.data,1); %keep track of flat electrodes if any
        
        %use asr method for cleaning
        OUTEEG = clean_asr(OUTEEG,20);
        
        %interpolate flat or bad electrodes if any
        OUTEEG = pop_interp(OUTEEG,temp_locs, 'spherical');
        
        %run ICA adjusting for data rank
        %number of components = number of non-interpolated - 1 (-1 is because of the common reference)
        OUTEEG= pop_runica(OUTEEG, 'icatype', 'runica','options',{'pca', 95-OUTEEG.int_electrodes-1} );
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

        %load behavioral data
        temp_name = [path, '/S_', num2str(s), '.mat'];
        load(temp_name, 'trial_condition', 'trial_load', 'reaction_time','correct_incorrect')

        %estimate accuracy per trial
        accuracy = nan(48,3);
        rt = nan(48,3);
        for t = 1:48
            for b = 1:4
                accuracy(t,b) = mean(correct_incorrect{t,b}) * 100;
                rt(t,b) = mean(reaction_time{t,b});
            end
        end
        
        
        %save individual behavior in eeglab structure 
        OUTEEG.accuracy = accuracy(:);
        OUTEEG.condition = cell2mat(trial_condition');
        OUTEEG.load = cell2mat(trial_load');
        OUTEEG.rt = rt(:);
        
        
        %save group behavior in metavar
        accuracy_l1_stay(s,1) = mean(OUTEEG.accuracy(and(OUTEEG.load==1, OUTEEG.condition==1)));
        accuracy_l3_stay(s,1) = mean(OUTEEG.accuracy(and(OUTEEG.load==3, OUTEEG.condition==1)));
        accuracy_l1_switch(s,1) = mean(OUTEEG.accuracy(and(OUTEEG.load==1, OUTEEG.condition==2)));
        accuracy_l3_switch(s,1) = mean(OUTEEG.accuracy(and(OUTEEG.load==3, OUTEEG.condition==2)));
        
        rt_l1_stay(s,1) = mean(OUTEEG.rt(and(OUTEEG.load==1, OUTEEG.condition==1)));
        rt_l3_stay(s,1) = mean(OUTEEG.rt(and(OUTEEG.load==3, OUTEEG.condition==1)));
        rt_l1_switch(s,1) = mean(OUTEEG.rt(and(OUTEEG.load==1, OUTEEG.condition==2)));
        rt_l3_switch(s,1) = mean(OUTEEG.rt(and(OUTEEG.load==3, OUTEEG.condition==2)));
        

        %save preprocessed data in a folder
        file_name = [path2save, 'S_', num2str(s), '_preprocessed.mat'];
        
        save(file_name, 'OUTEEG')
        
        %save in metavar number of interpolated electrodes, rejected components
        %and rejected epochs (to check if you need to have a closer look to any
        %subject)
        
        int_electrodes(s,1) = OUTEEG.int_electrodes ;
        rejected_components(s,1) = OUTEEG.rejected_components;
           
        
    end
    
 
    
end

%save metavars
 save('/Users/juliorodriguezlarios/Desktop/WM_EEG_study/variables/metavars', 'int_electrodes','rejected_components', 'accuracy_l1_stay', 'accuracy_l1_switch', 'accuracy_l3_stay', 'accuracy_l3_switch', 'rt_l1_stay', 'rt_l1_switch',...
     'rt_l3_stay', 'rt_l3_switch')
 


 
 
