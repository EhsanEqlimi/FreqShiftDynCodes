

%%%%%%%%% INPUT %%%%%%%%%%%%

% data - fieldtrip structure with epoched data:
% data.trial must contain a cell per trial with EEG/MEG time series
% data.fsample must contain sampling rate in Hz

% frequencies - vector including frequencies of interest and frequency resolution

%nmbcycles - minimum number of cycles to detect burst

%noiselimit - minimum and maximum power for beta bursts in z score
%e.g. [0 5] you get beta bursts in their power is above the mean beta power after 1/f correction and below 5std

% lowerfreq - lowest frequency to detect burst

%%%%%%%% OUTPUT %%%%%%%%%

%bamp - burst amplitude (prominence of peak) per trial and frequency
%bnumb - burst rate per trial and frequency
%bdur - burst duration per trial and frequency
%allbursts_spectrum - frequency spectrum of detected burst.
%allbursts_shape - waveform shape of detected bursts

function [btime,bamp, bnumb, bdur, allbursts_spectrum,allbursts_shape,btime_lowfreq,...
    bamp_lowfreq, bnumb_lowfreq, bdur_lowfreq, allbursts_spectrum_lowfreq,allbursts_shape_lowfreq]...
    = detect_bursts_maxpeak(data,frequencies,nmbcycles,noiselimit,lowerfreq)


%define time vector for later estimation of btime
timevec = linspace(0,3,length(data.time{1}));

%concatenate trials
data_conc =  cell2mat(data.trial);

%create empty output to fill
btime = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
bamp = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
bdur = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
bnumb = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
allbursts_spectrum =  nan(size(data_conc,1),size(data.trial,2), length(frequencies));
allbursts_shape =  nan(size(data_conc,1),size(data.trial,2), length(data.time{1})*2);

btime_lowfreq = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
bamp_lowfreq = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
bdur_lowfreq = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
bnumb_lowfreq = nan(size(data_conc,1), size(data.trial,2), length(frequencies));
allbursts_spectrum_lowfreq =  nan(size(data_conc,1),size(data.trial,2), length(frequencies));
allbursts_shape_lowfreq =  nan(size(data_conc,1),size(data.trial,2), length(data.time{1})*2);

for e = 1:size(data_conc,1) %loop electrodes
    
    %get time frequency of concatenated data
    [B,T,F]=BOSC_tf(data_conc(e,:),frequencies,data.fsample,6);
    %get average spectrum
    temp = mean(B,2);
    
    
    %%%%%%%%%%%% 1f correction %%%%%%%%%%%%
    
    %subtract main peak to improve fitting
    [pks,locs,w,p] = findpeaks(temp,frequencies);
    if ~isempty(pks) %if there is a peak
        %index main peak
        [~, index] = max(p);
        %get first and last frequency of the peak
        speak = locs(index) - w(index);
        endpeak = locs(index) + w(index);
        %define frequencies to include
        include = or(frequencies<speak , frequencies>endpeak);
    else
        include = 1:length(frequencies);
    end
    %remove peak from spectrum
    temp (~include) = nan;
    
    %fit trend
    try
        b = robustfit(log10(frequencies),log10(temp));
    catch %if not enough points after peak removal it wouold give you an error. In that case fit the raw spectrum.
        b = robustfit(log10(frequencies),log10(mean(B,2)));
    end
    
    %get 1/f power
    pv(1) = b(2);
    pv(2) = b(1);
    amplitude_1f=10.^(polyval(pv,log10(frequencies))); %get amplitude of 1/f trend
    
    %subtract 1/f from timefreq
    B =   B - repmat(amplitude_1f,size(B,2),1)';
    B(B<0) = 0; %get rid of power below 1/f
    
    
    
    %mark bursts above 1/f trend and below the noise power limit
    Bz = zscore(B,0,2);
    mask = B>0 & Bz>noiselimit(1) & Bz<noiselimit(2);
    
    %get length of each trial
    tlen =  cellfun(@length,data.trial);
    
    %cumulative sum is the end of each trial
    endtr = cumsum(tlen);
    
    %create temporal cells to fill
    bamp_temp = cell(size(data.trial,2), length(frequencies));
    bdur_temp = cell(size(data.trial,2), length(frequencies));
    brate_temp = cell(size(data.trial,2), length(frequencies));
    btime_temp = cell(size(data.trial,2), length(frequencies));
    burst_spectrum_temp= nan(size(data.trial,2), length(frequencies), 100,length(frequencies));
    bshape_temp= nan(size(data.trial,2), length(frequencies), 100,length(data.time{1})*2);
    
    bamp_temp2 = cell(size(data.trial,2), length(frequencies));
    bdur_temp2 = cell(size(data.trial,2), length(frequencies));
    brate_temp2 = cell(size(data.trial,2), length(frequencies));
    btime_temp2 = cell(size(data.trial,2), length(frequencies));
    burst_spectrum_temp2= nan(size(data.trial,2), length(frequencies), 100,length(frequencies));
    bshape_temp2= nan(size(data.trial,2), length(frequencies), 100,length(data.time{1})*2);
    
    %get time frequency matrices per trial and detect bursts
    for t = 1:size(data.trial,2) %loop trials
        %define beginning of trial by subtracting its length
        beg = endtr(t) - tlen(t) +1;
        %define trial
        tf_epoched = B(:,beg:endtr(t));
        mask_epoched = mask(:,beg:endtr(t));
        epoch = data_conc(e,beg:endtr(t));
        %sanity check epoching is working
        %sanity{t} = data_conc(:,beg:endtr(t));
        
        %exclude edges of time frequency matrix
        %wavelet is not reliable in the first / last 3 cycles of the time-frequency matrix
        edges = round(1./frequencies*3.*data.fsample);
        for fr = 1:size(mask_epoched,1)
            mask_epoched(fr,1:edges(fr)) = 0;
            mask_epoched(fr,size(mask_epoched,2)-edges(fr): size(mask_epoched,2)) = 0;
        end
        
        %add one more column of 0s so you can detect beg and end of burst
        %this way everything is shifted to the left by one sample
        mask_epoched2 = mask_epoched;
        mask_epoched2(:,size(mask_epoched,2)+1) = 0;
        
        %get temporal derivative to detect begin and end of bursts per trial
        %frequency and electrode
        Y = diff(mask_epoched2,1,2);
        
        for f = 1:size(Y,1) %loop frequencies
            
            if frequencies(f) >= lowerfreq % if this is a frequency of interest
                
                beginb = find(Y(f,:)==1)+1; %get burst beginning
                endb= find(Y(f,:)==-1);%get burst end
                
                if ~isempty(beginb) %if there are bursts
                    
                    %put nans in all possible bursts to be filled later
                    bamp_temp{t,f}(1:length(beginb))= nan;
                    bdur_temp{t,f}(1:length(beginb))= nan;
                    brate_temp{t,f}(1:length(beginb))= nan;
                    btime_temp{t,f}(1:length(beginb))= nan;
                    
                    bamp_temp2{t,f}(1:length(beginb))= nan;
                    bdur_temp2{t,f}(1:length(beginb))= nan;
                    brate_temp2{t,f}(1:length(beginb))= nan;
                    btime_temp2{t,f}(1:length(beginb))= nan;
                    
                    for bu = 1: length(beginb) %loop bursts
                        
                        %index begining and end of the burst in time points
                        index = beginb(bu):endb(bu);
                        
                        %get duration of burst in samples
                        samples = length(index);
                        %get 1 cycle duration at this frequency in samples
                        cycle_duration = 1/frequencies(f)*data.fsample;
                        
                        if samples>cycle_duration*nmbcycles %if the burst has at least X cycles

                            %%%% save time series of bursts
                            % we save a vector with the raw burst time series
                            % centered to their maximum peak
                            clearvars burst temp_burst
                            %get time series of burst
                            burst = epoch(index);
                            
                            %get maximum peak
                            [pks,locs,w,p] = findpeaks(burst);
                            [~, pos] = max(pks);
                            
                            %lock burst to maximum peak
                            center_bu =  locs(pos);
                            half2 = burst(center_bu:end);
                            half1 = burst(1:center_bu);
                            
                            center_temp = round(size(data.time{1},2));
                            temp_burst = nan(1,size(data.time{1},2)*2);
                            temp_burst(center_temp:center_temp+length(half2)-1) = half2;
                            temp_burst(center_temp-length(half1)+1:center_temp) = half1;
     
                            %get spectrum of burst to estimate peak amplitude
                            temp = mean(tf_epoched(:,index),2);
                            %find peaks
                            [pks,locs,w,p] = findpeaks(temp,frequencies);

                            %if there is a peak at that frequency
                            if ismember(frequencies(f),locs)
                                %save params when the burst peak has the maximum amplitude of spectrum
                                [~, pos] = max(pks); %get maximum peak
                                if locs(pos) == frequencies(f) %if the detected burst is the maximum peak
                                    bamp_temp{t,f}(bu)= p(locs==frequencies(f));%save prominence
                                    bdur_temp{t,f}(bu) = samples./data.fsample; % save duration in seconds
                                    brate_temp{t,f}(bu) = 1; % save unit for nmb bursts
                                    burst_spectrum_temp(t,f,bu,:)= temp;%save spectrum
                                    btime_temp{t,f}(bu)= timevec(index(1));%save start time
                                    bshape_temp(t,f,bu,:) = temp_burst;%save shape
                                    
                                end
                                
                                if locs(pos) < lowerfreq %if the detected burst is NOT the maximum peak
                                    bamp_temp2{t,f}(bu)= p(locs==frequencies(f));%save prominence
                                    bdur_temp2{t,f}(bu) = samples./data.fsample; % save duration in seconds
                                    brate_temp2{t,f}(bu) = 1; % save unit for nmb bursts
                                    burst_spectrum_temp2(t,f,bu,:)= temp;%save spectrum
                                    btime_temp2{t,f}(bu)= timevec(index(1)); %save start time
                                    bshape_temp2(t,f,bu,:) = temp_burst;%save waveform shape
                                end
                                
                                
                                
                                
                            end
                            
                            
                        end
                        
                        
                    end
                    
                end
                
            end
            
            
            
            
            
        end
    end

    %obtain genuine beta burst output by averaging burst dimension
    bamp(e,:,:) =  cellfun(@nanmean,bamp_temp);
    bdur(e,:,:) =  cellfun(@nanmean,bdur_temp);
    bnumb(e,:,:) = cellfun(@nansum,brate_temp);
    btime(e,:,:) = cellfun(@nanmean,btime_temp);
    allbursts_spectrum(e,:,:) = squeeze(nanmean(nanmean(burst_spectrum_temp,3),2));
    allbursts_shape(e,:,:) = squeeze(nanmean(nanmean(bshape_temp,3),2));
    
    %same for non-genuine beta bursts
    bamp_lowfreq(e,:,:) =  cellfun(@nanmean,bamp_temp2);
    bdur_lowfreq(e,:,:) =  cellfun(@nanmean,bdur_temp2);
    bnumb_lowfreq(e,:,:) = cellfun(@nansum,brate_temp2);
    btime_lowfreq(e,:,:) = cellfun(@nanmean,btime_temp2);
    allbursts_shape_lowfreq(e,:,:) = squeeze(nanmean(nanmean(bshape_temp2,3),2));
    allbursts_spectrum_lowfreq(e,:,:) = squeeze(nanmean(nanmean(burst_spectrum_temp2,3),2));
    
end




end



%%%%%%%%    This file is part of the Better OSCillation detection (BOSC) library %%%%%%%%%%%%%%%%.
%
%    The BOSC library is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The BOSC library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2010 Jeremy B. Caplan, Adam M. Hughes, Tara A. Whitten
%    and Clayton T. Dickson.

function [B,T,F]=BOSC_tf(eegsignal,F,Fsample,wavenumber)
% [B,T,F]=BOSC_tf(eegsignal,F,Fsample,wavenumber);
%
% This function computes a continuous wavelet (Morlet) transform on
% a segment of EEG signal; this can be used to estimate the
% background spectrum (BOSC_bgfit) or to apply the BOSC method to
% detect oscillatory episodes in signal of interest (BOSC_detect).
%
% parameters:
% eegsignal - a row vector containing a segment of EEG signal to be
%             transformed
% F - a set of frequencies to sample (Hz)
% Fsample - sampling rate of the time-domain signal (Hz)
% wavenumber is the size of the wavelet (typically, width=6)
%
% returns:
% B - time-frequency spectrogram: power as a function of frequency
%     (rows) and time (columns)
% T - vector of time values (based on sampling rate, Fsample)

st=1./(2*pi*(F/wavenumber));
A=1./sqrt(st*sqrt(pi));
B = zeros(length(F),length(eegsignal)); % initialize the time-frequency matrix
for f=1:length(F) % loop through sampled frequencies
    t=-3.6*st(f):(1/Fsample):3.6*st(f);
    m=A(f)*exp(-t.^2/(2*st(f)^2)).*exp(i*2*pi*F(f).*t); % Morlet wavelet
    y=conv(eegsignal,m); y=abs(y).^2;
    B(f,:)=y(ceil(length(m)/2):length(y)-floor(length(m)/2));
end
T=(1:size(eegsignal,2))/Fsample;
end