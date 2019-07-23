clear;clc;close all
load DEAP_chnames
ECG_srate = 128;
time_start = -.2; % time pre heartbeat (200 ms)
time_end = .6; % time post heartbeat (600 ms)
EEG_srate = 128;
points_start = abs(floor(time_start*EEG_srate)); % number of points before HB
points_end = abs(floor(time_end*EEG_srate)); % number of points after HB
epoch_length = points_start + points_end;
%%
for isub = 1:32 % subject index
    disp(isub);
    load(['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '.mat']) % load the full dataset, change directory accordingly
    % see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
    % data      40 x 40 x 8064	video/trial x channel (32 EEG + 8 phys) x data (63 s at 128 Hz)
    % labels	40 x 4          video/trial x label (valence, arousal, dominance, liking)
    
    % set to zero for each subject
    avg_HEP = zeros(40,epoch_length,32);
    clear HB_trials
    clear ratings
    %%
    for itrial = 1:40 % choose a trial        
        %%
        PS = squeeze(data(itrial,39,3*ECG_srate+1:end)); % plethysmograph (chan 39); discard the 3 s baseline
        time_ECG = (0:length(PS)-1)/ECG_srate;
        
        %% now let's detect R peaks from plethysmogram, we need to detrend the signal, and differentiate it twice
        [ddPS,PS] = PS2H(PS,ECG_srate); % [WS,PS_dz] (see PS2H.m)
        %%
        [pks,locs] = findpeaks(ddPS,ECG_srate,'MinPeakDistance',.6,'MinPeakHeight',.1); % look for R peaks at least .6 s away
        
        IBI = diff(locs);             % interbeat intervals
        min_IBI = min(IBI);           % see what the shortest one is
        locs(1) = []; locs(end) = []; % discard the first and last epoch
        IBI = diff(locs);             % interbeat intervals
        locs = locs-mean(IBI)/5;      % shift locations back in time since the second derivative of the plethysmogram signal comes later than the R peak
        n_epochs = length(locs);       
        
        %% now let's load the EEG
        eeg = squeeze(data(itrial,1:32,3*EEG_srate+1:end))'; % discard the 3 s baseline
        eeg = zscore(eeg);
        [npoints,nchan] = size(eeg);
        time_EEG = (0:length(eeg)-1)/EEG_srate;
        
        %% extract the epochs
        clear HB_EEG_epochs % clear each iteration because every trial has a different number of epochs
        
        HB_points = ceil(locs*EEG_srate); % points in the EEG time series when the heartbeats happen
        for i_epochs=1:n_epochs
            HB_EEG_epochs(i_epochs,:,:) = eeg(HB_points(i_epochs)-points_start:HB_points(i_epochs)+points_end-1,:); % same format as in Ince et al.'s code (trials x time (x electrodes))
        end
        
        avg_HEP(itrial,:,:) = squeeze(mean(HB_EEG_epochs,1)); % average HEP
        
        %% save individual epochs
        if itrial == 1
           HB_trials = HB_EEG_epochs;
           ratings = labels(repmat(itrial,n_epochs,1),:);
        else
           HB_trials = cat(1,HB_trials,HB_EEG_epochs);
           ratings = cat(1,ratings,labels(repmat(itrial,n_epochs,1),:));
        end
        
    end
    %% save
    
    % save average HEP per music video
    savefile = ['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP.mat'];
    save(savefile,'avg_HEP','labels');
    % avg_HEP    40 x 102 x 32    video/trial x data x channel
    % labels     40 x 4		      video/trial x label (valence, arousal, dominance, liking)
    
    % you can also save the individual epochs
    %savefile = ['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_HEP.mat'];
    %save(savefile,'HB_trials','ratings');
    % HB_trials  .... x 102 x 32  HB/trial x data x channel
    % ratings    .... x 4		  HB/trial x label (valence, arousal, dominance, liking)
    
end
