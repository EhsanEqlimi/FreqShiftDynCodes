 
function [FitParam,AperiodicFit,FlatPower,Alpha]=FnOneOverFLogLog(PowerMat,FreqVec,FOI,IsLogPower)
FitIdx=Freq >= FOI(1) & Freq <= FOI(2);
FreqSel=FreqVec(FitIdx);
LogFreq=log10(FreqSel);
if IsLogPower==1
    SpecSel=PowerMat(:,FitIdx)/10; % convert dB → log10(power)
else
    SpecSel=log10(PowerMat);
end
for Ch=1:size(SpecSel,1) % channel loop
    ChannelPower=SpecSel(Ch,:)'; %no need to log10 (spectopo is db)
    FitParam=polyfit(LogFreq,ChannelPower,1);
    AperiodicFit=polyval(FitParam,LogFreq);
    FlatPower(:,Ch)=ChannelPower-AperiodicFit;
    Alpha(Ch)=-FitParam(1); %kinda scaling exponent
end