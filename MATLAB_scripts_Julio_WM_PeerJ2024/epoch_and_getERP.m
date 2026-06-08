
 iclabel_threshold = 0.7;
 cleanchannels_threshold = 0.8;
%go to folder with preprocessed data
cd '/Users/juliorodriguezlarios/Desktop/Projects/categorization EEG/preprocessed_EEG'
files = dir('*.mat');
%epoch relative to the start of the delay (stimulus offset)
 delaystart = {'S101','S102','S103', 'S104', 'S105', 'S106','S107','S108',...
     'S201','S202','S203', 'S204', 'S205', 'S206','S207','S208',...
     'S 45','S 46','S 47', 'S 48', 'S 49', 'S 50','S 51','S 52'};
 
test_durations = [0.2 0.250 0.319 0.331 0.369 0.381 0.450 0.5;...
    0.45 0.5 0.619 0.669 0.706 0.756 0.870 0.920;...
    0.870 0.920 0.981 1.169 1.231 1.419 1.470 1.520];
 
EEGdata = nan(length(files), 96,  550,  336);
bad_trials = nan(length(files),336);
block_type= nan(length(files),336);
duration= nan(length(files),336);
long= nan(length(files),336);
correct= nan(length(files),336);
rt= nan(length(files),336);

for f = 1:length(files)
%add load EEG data and loop for subject
load(files(f).name)
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
temp_index = find(extractfield(EEG_epoched.event, 'epoch') ==e);
temp_events = extractfield(EEG_epoched.event(temp_index), 'type');
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
temp_latencies = extractfield(EEG_epoched.event(temp_index), 'latency');
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


end

save('/Users/juliorodriguezlarios/Desktop/Projects/categorization EEG/variables/erp_data','-v7.3')


