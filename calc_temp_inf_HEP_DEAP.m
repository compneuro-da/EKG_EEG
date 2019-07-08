clear;clc;close all
load DEAP_chnames
chan = 24; % choose an electrode (Cz)
time_start = -.2;
time_end = .6;
srate = 128;
time_epoch = time_start:1/srate:time_end-1/srate;
% time window we are interested in
timewindow = (time_epoch>.2) & (time_epoch<.6);
time_info = time_epoch(timewindow);
ntime_info = length(time_info);
%%
for isub = 1:32
    disp(isub);
    load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_HEP.mat']) % load data, change directory accordingly
    % HB_EEG_epochs  HBs x data points x EEG channels
    % ratings        HBs x labels (valence, arousal, dominance, liking)
    %% initialization
    dat = HB_EEG_epochs(:,:,chan);
    [ntrls,ntime_epoch] = size(dat);
    nlbls = size(ratings,2);
    
    hep = squeeze(mean(dat(:,timewindow),1));
    %plot(time_info,hep)

    %% copula transform and calculate mutual information (MI) at each time point
    cdat = copnorm(dat(:,timewindow));
    cratings = copnorm(ratings);
   
    MI = zeros(nlbls,ntime_info);
    for l=1:nlbls % loop over 4 labels
        for t=1:ntime_info
            MI(l,t) = mi_gg(cdat(:,t),cratings(:,l),true,true);
        end
    end

    %% cross-temporal interaction information (II)
    JMI = zeros(nlbls,1);
    II = zeros(ntime_info,ntime_info,nlbls);
    for l=1:nlbls
        for t1=1:ntime_info
            for t2=(t1+1):ntime_info
                JMI(l,:) = mi_gg([cdat(:,t1) cdat(:,t2)],cratings(:,l),true,true);
                II(t1,t2,l) = JMI(l,:) - MI(l,t1) - MI(l,t2);
            end
        end
        II(:,:,l) = II(:,:,l) + II(:,:,l)';
    end
    II(II==0)=NaN;
  
    %% partial information decomposition (PID)
    % see Faes et al., 2017; Ince, 2017
    % JMI = gcmi_cc(t,[d1 d2]); % t = target, d1 = source, d2 = source
    % MI1 = gcmi_cc(t,d1);
    % MI2 = gcmi_cc(t,d2);
    % RED = min(MI1,MI2);
    % U1 = MI1 - RED;
    % U2 = MI2 - RED;
    % SYN = JMI - RED - U1 - U2;
    % II = JMI - MI1 - MI2;
    
    JMI = zeros(nlbls,1);
    II = zeros(ntime_info,ntime_info,nlbls);
    RED = zeros(ntime_info,ntime_info,nlbls);
    U1 = zeros(nlbls,ntime_info);
    U2 = zeros(nlbls,ntime_info);
    SYN = zeros(ntime_info,ntime_info,nlbls);
    for l=1:nlbls
        for t1=1:ntime_info
            for t2=(t1+1):ntime_info
                JMI(l,:) = mi_gg([cdat(:,t1) cdat(:,t2)],cratings(:,l),true,true);
                II(t1,t2,l) = JMI(l,:) - MI(l,t1) - MI(l,t2);
         
                RED(t1,t2,l) = min(MI(l,t1),MI(l,t2));
                U1(l,t1) = MI(l,t1) - RED(t1,t2,l);
                U2(l,t2) = MI(l,t2) - RED(t1,t2,l);
                SYN(t1,t2,l) = JMI(l,:) - RED(t1,t2,l) - U1(l,t1) - U2(l,t2);
            end
        end
        II(:,:,l) = II(:,:,l) + II(:,:,l)';
        SYN(:,:,l) = SYN(:,:,l) + SYN(:,:,l)';
        RED(:,:,l) = RED(:,:,l) + RED(:,:,l)';
    end
    II(II==0)=NaN;
    SYN(SYN==0)=NaN;
    RED(RED==0)=NaN;
    
    %% save results
    %savefile = ['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_I.mat'];
    %save(savefile,'hep','MI','II','SYN','RED');
end
