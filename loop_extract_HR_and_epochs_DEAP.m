clear;clc;close all
load DEAP_chnames
ECG_srate=128;
time_start=-.2; % time pre heartbeat (200 ms)
time_end=.6; %time post heartbeat (300 ms)
EEG_srate=128;
points_start=abs(floor(time_start*EEG_srate)); %number of points before HB
points_end=abs(floor(time_end*EEG_srate)); %number of points after HB
epoch_length=points_start+points_end;
for isub=1:32; %subject index, you can loop on this
    disp(isub);
    load(['C:\Users\dmarinaz\Documents\code\DEAP\preprocessed\s' num2str(isub,'%02.0f') '.mat']) %load the full dataset, change directory accordingly
    % see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
    % data	40 x 40 x 8064	video/trial x channel x data
    % labels	40 x 4	video/trial x label (valence, arousal, dominance, liking)
    % channels have a different order in subs 23 to 32, here the index to
    % reorder them. We choose the Twente order and name.
    for itrial=1:40; % here you choose a trial, eventually you can loop on it
        
        PS=squeeze(data(itrial,39,:)); % plethysmograph;
        time_ECG=(0:length(PS)-1)/ECG_srate;
        %% now let's detect R peaks from plethysmogram, we need to detrend the signal, and differentiate it twice
        [ddPS, PS]  = PS2H( PS, ECG_srate );
        %%
        [pks,locs] = findpeaks(ddPS,ECG_srate,'MinPeakDistance',.6,'MinPeakHeight',.1);% look for peaks at least .4 sec away
        
        IBI=diff(locs); %interbeat interval
        min_IBI=min(IBI); %let's see what the shortest one is
        locs(1)=[];locs(end)=[]; % let's discard the first and last epoch
        IBI=diff(locs); %interbeat interval
        locs=locs-mean(IBI)/5; % shift locations back in time since the second derivative of the Plethysmogram sugnal comes later than the R peak
        n_epochs=length(locs);
        %% now let's load the EEG
        eeg=squeeze(data(itrial,1:32,:))';
        eeg=zscore(eeg);
        [npoints, nchan]=size(eeg);
        %epoch_length=floor(min_IBI*EEG_srate); % the length of the heartbeat evoked potential is equal to the shortest interbeat interval
        time_EEG=(0:length(eeg)-1)/EEG_srate;
        
        %% extract the epochs
        HB_points=ceil(locs*EEG_srate); %points in the EEG time series when the heartbeats happen
        HR_EEG_epochs=zeros(40,nchan,epoch_length,n_epochs); % vector of epochs, for 40 trials, same format as in Cohen's book
        
        for i_epochs=1:n_epochs
            HR_EEG_epochs(itrial,:,:,i_epochs)=eeg(HB_points(i_epochs)-points_start:HB_points(i_epochs)+points_end-1,:)'; % transpose since channels have to be the first dimension
        end
        
    end
    %% save
    savefile=['C:\Users\dmarinaz\Documents\code\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_HEP.mat'];
    save(savefile,'HR_EEG_epochs');
end
