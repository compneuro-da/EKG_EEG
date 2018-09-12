clear;clc;close all
sub=12; %subject index, you can loop on this
load(['C:\Users\dmarinaz\Documents\code\DEAP\preprocessed\s' num2str(sub,'%02.0f') '.mat']) %load the full dataset, change directory accordingly
% see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
% data	40 x 40 x 8064	video/trial x channel x data
% labels	40 x 4	video/trial x label (valence, arousal, dominance, liking)
for itrial=12%:40 % here you choose a trial, eventually you can loop on it
    ECG_srate=128;
    time_start=-.2; % time pre heartbeat (200 ms)
    time_end=.3; %time post heartbeat (300 ms)
    PS=squeeze(data(itrial,39,3*ECG_srate+1:end)); % plethysmograph; let's discard the 3 seconds baseline
    time_ECG=(0:length(PS)-1)/ECG_srate;
    %figure;plot(time_ECG,PS);title('plethysmograph signal');
    %% now let's detect R peaks from plethysmogram, we need to detrend the signal, and differentiate it twice
    [ddPS, PS]  = PS2H( PS, ECG_srate );
    %%
    [pks,locs] = findpeaks(ddPS,ECG_srate,'MinPeakDistance',.6,'MinPeakHeight',.1);% look for peaks at least .4 sec away
    IBI=diff(locs); %interbeat interval
    locs=locs-mean(IBI)/5; % shift locations back in time since the second derivative of the Plethysmogram sugnal comes later than the R peak
    figure;
    subplot(2,1,1);plot(time_ECG(2:end-1),ddPS);hold on;scatter(locs,pks,'r');title('R peaks on second derivative of PS')
    subplot(2,1,2);plot(time_ECG(2:end-1),PS(2:end-1));hold on;scatter(locs,PS(round(locs*ECG_srate)),'r');title('R peaks on PS')
    suptitle(['sub' ' ' num2str(sub,'%02.0f') ', trial' ' ' num2str(itrial)])
    
    figure;
    plot(time_ECG(2:end-1),ddPS);hold on;scatter(locs,pks,'r');title('R peaks based on second derivative of PS')
    hold on;plot(time_ECG(2:end-1),PS(2:end-1),'k');hold on;scatter(locs,PS(round(locs*ECG_srate)),'r');
    title(['sub' ' ' num2str(sub,'%02.0f') ', trial' ' ' num2str(itrial)])
end