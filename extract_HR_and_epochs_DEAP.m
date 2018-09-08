clear;clc;close all
load('C:\Users\dmarinaz\Documents\code\DEAP\preprocessed\s01.mat') %load the full dataset, change directory accordingly
% see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
% data	40 x 40 x 8064	video/trial x channel x data
% labels	40 x 4	video/trial x label (valence, arousal, dominance, liking)
itrial=1; % here you choose a trial, eventually you can loop on it
ECG_srate=128;
time_start=-.2; % time pre heartbeat (200 ms)
time_end=.3; %time post heartbeat (200 ms)
PS=squeeze(data(itrial,39,3*ECG_srate+1:end)); % plethysmograph; let's discard the 3 seconds baseline
time_ECG=(0:length(PS)-1)/ECG_srate;
figure;plot(time_ECG,PS);title('plethysmograph signal');
%% now let's detect R peaks from plethysmogram, we need to detrend the signal, and differentiate it twice
lambda_max = l1tf_lambdamax(PS);
[trend,status] = l1tf(PS, 0.001*lambda_max);
PS_d=PS-trend;
PS_dz=zscore(PS_d);
ddPS=smooth(diff(PS_dz,2));
%%
[pks,locs] = findpeaks(ddPS,ECG_srate,'MinPeakDistance',.4,'MinPeakHeight',.005);% look for peaks at least .4 sec away
figure;plot(time_ECG(2:end-1),ddPS);hold on;scatter(locs,pks,'r');title('see if we got the peaks right')
IBI=diff(locs); %interbeat interval
min_IBI=min(IBI); %let's see what the shortest one is
locs(1)=[];locs(end)=[]; % let's discard the first and last epoch
n_epochs=length(locs);
%% now let's load the EEG
EEG_srate=128;
points_start=abs(floor(time_start*EEG_srate)); %number of points before HB
points_end=abs(floor(time_end*EEG_srate)); %number of points after HB
eeg=squeeze(data(itrial,1:32,3*EEG_srate+1:end))'; % let's discard the 3 seconds baseline
[npoints, nchan]=size(eeg);
%epoch_length=floor(min_IBI*EEG_srate); % the length of the heartbeat evoked potential is equal to the shortest interbeat interval
epoch_length=floor((time_end-time_start)*EEG_srate);
time_EEG=(0:length(eeg)-1)/EEG_srate;

%% extract the epochs
HB_points=ceil(locs*EEG_srate); %points in the EEG time series when the heartbeats happen
HR_EEG_epochs=zeros(nchan,epoch_length,n_epochs); % vector of epochs, same format as in Cohen's book

for i_epochs=1:n_epochs
    HR_EEG_epochs(:,:,i_epochs)=eeg(HB_points(i_epochs)-points_start:HB_points(i_epochs)+points_end-1,:)'; % transpose since channels have to be the first dimension
end
avg_HRP=squeeze(mean(HR_EEG_epochs,3)); % average to find the HRP
time_epoch=[-.2:1/EEG_srate:.3-1/EEG_srate];
figure;
for i_chan=1:nchan
    subplot(6,6,i_chan);plot(time_epoch,avg_HRP(i_chan,:));% to add the channel name we should first define the labels
    xlim([time_start,time_end]);
end