clear;clc;close all
load DEAP_chnames
chan = 24; % choose an electrode
time_start = -.2;
time_end = .6;
srate = 128;
time_epoch = time_start:1/srate:time_end-1/srate;
ntime_epoch = length(time_epoch);
% time window we are interested in
timewindow = (time_epoch>.2) & (time_epoch<.6);
time_info = time_epoch(timewindow);
ntime_info = length(time_info);

%% calculate temporal information per subject
for isub = 1:32
    disp(isub);
    load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP.mat']) % load data, change directory accordingly
    %load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_HEP.mat'])
    
    if isub == 1
       avg_HEP_tot = avg_HEP;
       labels_tot = labels;
       %HB_trials_tot = HB_trials;
       %ratings_tot = ratings;
    else
       avg_HEP_tot = cat(1,avg_HEP_tot,avg_HEP);
       labels_tot = cat(1,labels_tot,labels);
       %HB_trials_tot = cat(1,HB_trials_tot,HB_trials);
       %ratings_tot = cat(1,ratings_tot,ratings);
    end
    
    dat = avg_HEP(:,timewindow,chan);
    %dat = HB_trials(:,timewindow,chan);
    [ntrls,~] = size(dat);
 
    HEP = squeeze(mean(dat,1));
    %plot(time_info,HEP)
 
    %% copula transform and calculate mutual information (MI) at each time point
    cdat = copnorm(dat);
    cratings = copnorm(labels);
    %cratings = copnorm(ratings);
    nlbls = size(cratings,2);
   
    MI = zeros(nlbls,ntime_info);
    for l=1:nlbls % loop over 4 labels
        for t=1:ntime_info
            MI(l,t) = mi_gg(cdat(:,t),cratings(:,l),false);
        end
    end

    %% cross-temporal interaction information (II)
    JMI = zeros(nlbls,1);
    II = zeros(ntime_info,ntime_info,nlbls);
    for l=1:nlbls
        for t1=1:ntime_info
            for t2=(t1+1):ntime_info
                JMI(l,:) = mi_gg([cdat(:,t1) cdat(:,t2)],cratings(:,l),false);
                II(t1,t2,l) = JMI(l,:) - MI(l,t1) - MI(l,t2);
            end
        end
        II(:,:,l) = II(:,:,l) + II(:,:,l)';
    end
    II(II==0)=NaN;
  
    %% partial information decomposition (PID)
    % see Faes, Marinazzo, & Stramaglia (2017)
    % JMI = gcmi_cc(t,[d1 d2]) with t = target, d1 = source, d2 = source
    % MI1 = gcmi_cc(t,d1)
    % MI2 = gcmi_cc(t,d2)
    % RED = min(MI1,MI2)
    % U1 = MI1 - RED
    % U2 = MI2 - RED
    % SYN = JMI - RED - U1 - U2
    % II = JMI - MI1 - MI2
    
    JMI = zeros(nlbls,1);
    II = zeros(ntime_info,ntime_info,nlbls);
    RED = zeros(ntime_info,ntime_info,nlbls);
    U1 = zeros(nlbls,ntime_info);
    U2 = zeros(nlbls,ntime_info);
    SYN = zeros(ntime_info,ntime_info,nlbls);
    for l=1:nlbls
        for t1=1:ntime_info
            for t2=(t1+1):ntime_info
                JMI(l,:) = mi_gg([cdat(:,t1) cdat(:,t2)],cratings(:,l),false);
                II(t1,t2,l) = JMI(l,:) - MI(l,t1) - MI(l,t2);
         
                RED(t1,t2,l) = min(MI(l,t1),MI(l,t2));
                U1(l,t1) = MI(l,t1) - RED(t1,t2,l);
                U2(l,t2) = MI(l,t2) - RED(t1,t2,l);
                SYN(t1,t2,l) = JMI(l,:) - U1(l,t1) - U2(l,t2) - RED(t1,t2,l);
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
    savefile = ['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP_I.mat'];
    %savefile = ['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_HEP_I.mat'];
    save(savefile,'HEP','MI','II','SYN','RED');
end

%% calculate temporal information over all subjects
dat_tot = avg_HEP_tot(:,timewindow,chan);
%dat_tot = HB_trials_tot(:,timewindow,chan);
[ntrls_tot,~] = size(dat_tot);

HEP_tot = squeeze(mean(dat_tot,1));
%plot(time_info,HEP_tot)

% copula transform and calculate mutual information (MI) at each time point
cdat_tot = copnorm(dat_tot);
cratings_tot = copnorm(labels_tot);
%cratings_tot = copnorm(ratings_tot);

MI_tot = zeros(nlbls,ntime_info);
for l=1:nlbls % loop over 4 labels
    for t=1:ntime_info
        MI_tot(l,t) = mi_gg(cdat_tot(:,t),cratings_tot(:,l),false);
    end
end
  
% cross-temporal interaction information (II) and PID
JMI_tot = zeros(nlbls,1);
II_tot = zeros(ntime_info,ntime_info,nlbls);
RED_tot = zeros(ntime_info,ntime_info,nlbls);
U1_tot = zeros(nlbls,ntime_info);
U2_tot = zeros(nlbls,ntime_info);
SYN_tot = zeros(ntime_info,ntime_info,nlbls);
for l=1:nlbls
    for t1=1:ntime_info
        for t2=(t1+1):ntime_info
            JMI_tot(l,:) = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,l),false);
            II_tot(t1,t2,l) = JMI_tot(l,:) - MI_tot(l,t1) - MI_tot(l,t2);
     
            RED_tot(t1,t2,l) = min(MI_tot(l,t1),MI_tot(l,t2));
            U1_tot(l,t1) = MI_tot(l,t1) - RED_tot(t1,t2,l);
            U2_tot(l,t2) = MI_tot(l,t2) - RED_tot(t1,t2,l);
            SYN_tot(t1,t2,l) = JMI_tot(l,:) - U1_tot(l,t1) - U2_tot(l,t2) - RED_tot(t1,t2,l);
        end
    end
    II_tot(:,:,l) = II_tot(:,:,l) + II_tot(:,:,l)';
    SYN_tot(:,:,l) = SYN_tot(:,:,l) + SYN_tot(:,:,l)';
    RED_tot(:,:,l) = RED_tot(:,:,l) + RED_tot(:,:,l)';
end
II_tot(II_tot==0)=NaN;
SYN_tot(SYN_tot==0)=NaN;
RED_tot(RED_tot==0)=NaN;

%% save results
savefile = '\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\tot_avgHEP_I.mat';
%savefile = '\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\tot_HEP_I.mat';
save(savefile,'HEP_tot','MI_tot','II_tot','SYN_tot','RED_tot');
