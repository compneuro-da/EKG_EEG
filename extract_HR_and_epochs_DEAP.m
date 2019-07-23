clear;clc;close all
load DEAP_chnames
sub = 1; % subject index
load(['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(sub,'%02.0f') '.mat']) % load the full dataset, change directory accordingly
% see http://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
% data      40 x 40 x 8064	video/trial x channel (32 EEG + 8 phys) x data (63 s at 128 Hz)
% labels	40 x 4          video/trial x label (valence, arousal, dominance, liking)
itrial = 1; % choose a trial

%%
ECG_srate = 128;
PS = squeeze(data(itrial,39,3*ECG_srate+1:end)); % plethysmograph (chan 39); discard the 3 s baseline
time_ECG = (0:length(PS)-1)/ECG_srate;
time_start = -.2; % time pre heartbeat (200 ms)
time_end = .6; % time post heartbeat (600 ms)

%% now let's detect R peaks from plethysmogram, we need to detrend the signal, and differentiate it twice
[ddPS,PS] = PS2H(PS,ECG_srate); % [WS,PS_dz] (see PS2H.m)
%%
[pks,locs] = findpeaks(ddPS,ECG_srate,'MinPeakDistance',.6,'MinPeakHeight',.1); % look for peaks at least .6 sec away

IBI = diff(locs);             % interbeat intervals
min_IBI = min(IBI);           % see what the shortest one is
locs(1) = []; locs(end) = []; % discard the first and last epoch
IBI = diff(locs);             % interbeat intervals
locs = locs-mean(IBI)/5;      % shift locations back in time since the second derivative of the plethysmogram signal comes later than the R peak
n_epochs = length(locs);

%% now let's load the EEG
EEG_srate = 128;
eeg = squeeze(data(itrial,1:32,3*EEG_srate+1:end))'; % discard the 3 s baseline
eeg = zscore(eeg);
[npoints,nchan] = size(eeg);
time_EEG = (0:length(eeg)-1)/EEG_srate;

points_start = abs(floor(time_start*EEG_srate)); % number of points before HB
points_end = abs(floor(time_end*EEG_srate)); % number of points after HB
%epoch_length = floor(min_IBI*EEG_srate); % set the length of the heartbeat evoked potential equal to the shortest interbeat interval so that the epochs won't overlap
epoch_length = points_start + points_end;

%% extract the epochs
HB_points = ceil(locs*EEG_srate); % points in the EEG time series when the heartbeats happen
HB_EEG_epochs = zeros(nchan,epoch_length,n_epochs); % vector of epochs, same format as in Cohen's book (electrodes x time x trials)
for i_epochs = 1:n_epochs
    HB_EEG_epochs(:,:,i_epochs) = eeg(HB_points(i_epochs)-points_start:HB_points(i_epochs)+points_end-1,:)'; % transpose since channels have to be the first dimension
end

avg_HEP = squeeze(mean(HB_EEG_epochs,3)); % average HEP
std_HEP = squeeze(std(HB_EEG_epochs,[],3)); % std of the HEP

%% the lines below are if you want to plot all the EEG channels
time_epoch = time_start:1/EEG_srate:time_end;
if length(time_epoch)>epoch_length
   time_epoch = time_start:1/EEG_srate:time_end-1/EEG_srate;
end

figure;
for ichan = 1:nchan
    subplot(4,8,ichan); plot(time_epoch,avg_HEP(ichan,:));
    xlim([time_start,time_end]); ylim([min(min(avg_HEP)) max(max(avg_HEP))]);
    title(ch_labels{ichan})
end
