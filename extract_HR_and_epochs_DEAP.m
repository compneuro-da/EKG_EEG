clear;clc;close all
load DEAP_chnames
sub=2; %subject index, you can loop on this
load(['C:\Users\dmarinaz\Documents\code\DEAP\preprocessed\s' num2str(sub,'%02.0f') '.mat']) %load the full dataset, change directory accordingly
% see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
% data	40 x 40 x 8064	video/trial x channel x data
% labels	40 x 4	video/trial x label (valence, arousal, dominance, liking)

itrial=1; % here you choose a trial, eventually you can loop on it
ECG_srate=128;
time_start=-.2; % time pre heartbeat (200 ms)
time_end=.6; %time post heartbeat (300 ms)
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
EEG_srate=128;
points_start=abs(floor(time_start*EEG_srate)); %number of points before HB
points_end=abs(floor(time_end*EEG_srate)); %number of points after HB
eeg=squeeze(data(itrial,1:32,:))'; 
eeg=zscore(eeg);
[npoints, nchan]=size(eeg);
%epoch_length=floor(min_IBI*EEG_srate); % the length of the heartbeat evoked potential is equal to the shortest interbeat interval
epoch_length=points_start+points_end;
time_EEG=(0:length(eeg)-1)/EEG_srate;

%% extract the epochs
HB_points=ceil(locs*EEG_srate); %points in the EEG time series when the heartbeats happen
HR_EEG_epochs=zeros(nchan,epoch_length,n_epochs); % vector of epochs, same format as in Cohen's book

for i_epochs=1:n_epochs
    HR_EEG_epochs(:,:,i_epochs)=eeg(HB_points(i_epochs)-points_start:HB_points(i_epochs)+points_end-1,:)'; % transpose since channels have to be the first dimension
end
avg_HRP=squeeze(mean(HR_EEG_epochs,3)); % average HRP
std_HRP=squeeze(std(HR_EEG_epochs,[],3)); % std of the the HRP
time_epoch=time_start:1/EEG_srate:time_end;
if length(time_epoch)>epoch_length
    time_epoch=time_start:1/EEG_srate:time_end-1/EEG_srate;
end

% the lines below are if you want to plot all the EEG channels
figure;
for i_chan=1:nchan
    subplot(4,8,i_chan);plot(time_epoch,avg_HRP(i_chan,:));
    xlim([time_start,time_end]);ylim([min(min(avg_HRP)) max(max(avg_HRP))]);
    title(ch_labels{i_chan})
end
