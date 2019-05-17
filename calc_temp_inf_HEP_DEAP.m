clear;clc;close all
load DEAP_chnames
isub = 1;
for isub = 1%:32 % subject index
    load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_HEP.mat']) % load data, change directory accordingly
    % HB_EEG_epochs  HB/trial x time x channel
    % ratings        HB/trial x label (valence, arousal, dominance, liking)
    chan = 24; % Cz
    dat = HB_EEG_epochs(:,:,chan); % HB epochs for Cz electrode
    [ntrls,ntime] = size(dat);
    time_start = -.2;
    time_end = .6;
    srate = 128;
    time = time_start:1/srate:time_end-1/srate;
    
    hep = squeeze(mean(dat,1));
    plot(time,hep) 
    
    %% copula transform and MI
    cdat = copnorm(dat);
    crat = copnorm(ratings);
    
    I = zeros(1,ntime);
    for t=1:ntime
        I(t) = mi_gg(cdat(:,t),crat(:,1),true,true);
    end
    plot(time,I)
    
    %% cross-temporal interaction information
    t1 = 1;
    t2 = 2;
    Ijoint = mi_gg([cdat(:,t1) cdat(:,t2)],crat(:,1),true,true);
   
    timewindow = (time>.2) & (time<.6);
    time_int = time(timewindow);
    ntime_int = length(time_int);
    
    cdat = copnorm(dat(:,timewindow));
    I = I(timewindow);
    
    I_int = zeros(ntime_int,ntime_int);
    for t1=1:ntime_int
        for t2=(t1+1):ntime_int
            Ijoint = mi_gg([cdat(:,t1) cdat(:,t2)],crat(:,1),true,true);
            I_int(t1,t2) = Ijoint - I(t1) - I(t2);
        end
    end
    I_int = I_int + I_int';
    
    figure
    imagesc(time_int,time_int,I_int)
    colormap parula
    colorbar
    axis square
    xlabel('Time (ms)')
    ylabel('Time (ms)')
    title('Interaction information (bits)')
    
end
