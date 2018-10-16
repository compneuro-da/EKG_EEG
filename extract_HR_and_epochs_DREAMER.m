load('\\Client\C$\Users\Liesa\Desktop\EEG\data\DREAMER.mat') %load the full dataset, change directory accordingly
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
epoch_length=floor(min_IBI*EEG_srate); % the length of the hartbeat evoked potential is equal to the shortest interbeat interval
time_EEG=(0:length(eeg)-1)/EEG_srate;

%% extract the epochs
HB_points=ceil(locs*EEG_srate); %points in the EEG time series when the heartbeats happen
HR_EEG_epochs=zeros(nchan,epoch_length,n_epochs); % vector of epochs, same format as in Cohen's book

for i_epochs=1:n_epochs
    HR_EEG_epochs(:,:,i_epochs)=eeg(HB_points(i_epochs):HB_points(i_epochs)+epoch_length-1,:)'; % transpose since channels have to be the first dimension
end
avg_HRP=squeeze(mean(HR_EEG_epochs,3)); % average to find the HRP on 3rd dimension
time_epoch=(0:epoch_length-1)/EEG_srate;
figure;
for i_chan=1:14
    subplot(2,7,i_chan);plot(time_epoch,avg_HRP(i_chan,:));title(DREAMER.EEG_Electrodes{i_chan});
end
