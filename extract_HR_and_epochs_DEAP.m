clear;clc;close all
load('C:\Users\dmarinaz\Documents\code\DEAP\preprocessed\s01.mat') %load the full dataset, change directory accordingly
% see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
% data	40 x 40 x 8064	video/trial x channel x data
% labels	40 x 4	video/trial x label (valence, arousal, dominance, liking)
itrial=1; % here you choose a trial, eventually you can loop on it
ECG_srate=128;
PS=squeeze(data(itrial,39,3*ECG_srate+1:end)); % plethysmograph; let's discard the 3 seconds baseline
time_ECG=(0:length(PS)-1)/ECG_srate;
figure;plot(time_ECG,PS);title('plethysmograph signal');
%% now let's detect R peaks, we need to detrend the signal
lambda_max = l1tf_lambdamax(PS);
[trend,status] = l1tf(PS, 0.001*lambda_max);
PS_d=PS-trend;
PS_dz=zscore(PS_d);
[pks,locs] = findpeaks(PS_dz,ECG_srate,'MinPeakDistance',.4,'MinPeakHeight',.01);% look for peaks at least .4 sec away
figure;plot(time_ECG,PS_dz);hold on;scatter(locs,pks,'r');title('see if we got the peaks right')
IBI=diff(locs); %interbeat interval
min_IBI=min(IBI); %let's see what the shortest one is
n_epochs=length(locs)-1; % let's discard the last epoch, since the last hearbeat could be too close to the end
%%now let's load the EEG
EEG_srate=128;
eeg=squeeze(data(itrial,1:32,3*EEG_srate+1:end))'; % let's discard the 3 seconds baseline
[npoints, nchan]=size(eeg);
epoch_length=floor(min_IBI*EEG_srate); % the length of the hartbeat evoked potential is equal to the shortest interbeat interval
time_EEG=(0:length(eeg)-1)/EEG_srate;

%% extract the epochs
HB_points=ceil(locs*EEG_srate); %points in the EEG time series when the heartbeats happen
HR_EEG_epochs=zeros(nchan,epoch_length,n_epochs); % vector of epochs, same format as in Cohen's book

for i_epochs=1:n_epochs
    HR_EEG_epochs(:,:,i_epochs)=eeg(HB_points(i_epochs):HB_points(i_epochs)+epoch_length-1,:)'; % transpose since channels have to be the first dimension
end
avg_HRP=squeeze(mean(HR_EEG_epochs,3)); % average to find the HRP
time_epoch=(0:epoch_length-1)/EEG_srate;
figure;
for i_chan=1:nchan
    subplot(6,6,i_chan);plot(time_epoch,avg_HRP(i_chan,:));% to add the channel name we should first define the labels
end