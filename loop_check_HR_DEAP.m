clear;clc;close all
load DEAP_chnames
sub = 1; % subject index
load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(sub,'%02.0f') '.mat']) % load the full dataset, change directory accordingly
% see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
% data      40 x 40 x 8064	video/trial x channel (32 EEG + 8 phys) x data (63 s at 128 Hz)
% labels	40 x 4          video/trial x label (valence, arousal, dominance, liking)
%%
for itrial = 1%:40 % choose a trial
    %%
    ECG_srate = 128;
    PS = squeeze(data(itrial,39,3*ECG_srate+1:end)); % plethysmograph (chan 39); discard the 3 s baseline
    time_ECG = (0:length(PS)-1)/ECG_srate;
    time_start = -.2; % time pre heartbeat (200 ms)
    time_end = .6; % time post heartbeat (600 ms)
    figure; plot(time_ECG,PS); title('plethysmograph signal');
    
    %% now let's detect R peaks from plethysmogram, we need to detrend the signal, and differentiate it twice
    [ddPS,PS] = PS2H(PS,ECG_srate); % [WS,PS_dz] (see PS2H.m)
    figure; plot(time_ECG,PS); title('plethysmograph signal (detrended)');
    time_ddPS = (0:length(ddPS)-1)/ECG_srate;
    figure; plot(time_ddPS,ddPS); title('plethysmograph signal (second derivative)');
    %%
    [pks,locs] = findpeaks(ddPS,ECG_srate,'MinPeakDistance',.6,'MinPeakHeight',.1); % look for R peaks at least .6 sec away
    IBI = diff(locs);        % interbeat intervals
    locs = locs-mean(IBI)/5; % shift locations back in time since the second derivative of the plethysmogram signal comes later than the R peak
    
    figure;
    subplot(2,1,1); plot(time_ECG(2:end-1),ddPS); hold on; scatter(locs,pks,'r'); title('R peaks on second derivative of PS')
    subplot(2,1,2); plot(time_ECG(2:end-1),PS(2:end-1)); hold on; scatter(locs,PS(round(locs*ECG_srate)),'r'); title('R peaks on PS')
    suptitle(['sub' ' ' num2str(sub,'%02.0f') ', trial' ' ' num2str(itrial)])
    
    figure;
    plot(time_ECG(2:end-1),ddPS); hold on; scatter(locs,pks,'r'); title('R peaks based on second derivative of PS')
    hold on; plot(time_ECG(2:end-1),PS(2:end-1),'k'); hold on; scatter(locs,PS(round(locs*ECG_srate)),'r');
    title(['sub' ' ' num2str(sub,'%02.0f') ', trial' ' ' num2str(itrial)])
end
