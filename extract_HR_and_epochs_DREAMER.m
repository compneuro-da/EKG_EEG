load('C:\Users\dmarinaz\Documents\code\DREAMER\DREAMER.mat') %load the full dataset, change directory accordingly
ecg=cell2mat(DREAMER.Data{1,1}.ECG.baseline(1)); %load the ecg for a given subject and condition
ecg(:,2)=[]; %use only one ecg channel
ECG_srate=DREAMER.ECG_SamplingRate;
time_ECG=(0:length(ecg)-1)/ECG_srate;
figure;plot(time_ECG,ecg)
%% now let's detect R peaks
x=zscore(ecg); %z score
figure;hist(x,100);
% you can see that the main peak is the subthreshold values, so the peak is
% at the extremities
[pks,locs] = findpeaks(x,ECG_srate,'MinPeakDistance',.4,'MinPeakHeight',3); % look for peaks at least .4 sec away, and 3 SD high
figure;plot(time_ECG,x);hold on;scatter(locs,pks,'r');title('see if we got the peaks right')
IBI=diff(locs); %interbeat interval
min_IBI=min(IBI); %let's see what the shortest one is
n_epochs=length(locs)-1; % let's discard the last epoch, since the last hearbeat could be too close to the end
%%now let's load the EEG
eeg=cell2mat(DREAMER.Data{1,1}.EEG.baseline(1));
[npoints, nchan]=size(eeg);
EEG_srate=DREAMER.EEG_SamplingRate;
time_start=-.2; % time pre heartbeat (200 ms)
time_end=.6; %time post heartbeat (300 ms)
%epoch_length=floor(min_IBI*EEG_srate); % the length of the hartbeat evoked potential is equal to the shortest interbeat interval
time_EEG=(0:length(eeg)-1)/EEG_srate;
points_start=abs(floor(time_start*EEG_srate)); %number of points before HB
points_end=abs(floor(time_end*EEG_srate)); %number of points after HB
epoch_length=points_start+points_end;
%% extract the epochs
HB_points=ceil(locs*EEG_srate); %points in the EEG time series when the heartbeats happen
HR_EEG_epochs=zeros(nchan,epoch_length,n_epochs); % vector of epochs, same format as in Cohen's book

for i_epochs=1:n_epochs
    HR_EEG_epochs(:,:,i_epochs)=eeg(HB_points(i_epochs):HB_points(i_epochs)+epoch_length-1,:)'; % transpose since channels have to be the first dimension
end
avg_HRP=squeeze(mean(HR_EEG_epochs,3)); % average to find the HRP on 3rd dimension
time_epoch=time_start:1/EEG_srate:time_end;
if length(time_epoch)>epoch_length
    time_epoch=time_start:1/EEG_srate:time_end-1/EEG_srate;
end
figure;
for i_chan=1:14
    subplot(2,7,i_chan);plot(time_epoch,avg_HRP(i_chan,:));title(DREAMER.EEG_Electrodes{i_chan});
    xlim([time_start,time_end]);%ylim([-8 8])
end
end
