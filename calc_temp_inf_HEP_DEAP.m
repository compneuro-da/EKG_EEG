clear;clc;close all
load DEAP_chnames
isub = 1;

%for isub = 1:32 % subject index
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
    
    % time window we are interested in
    timewindow = (time>.2) & (time<.6);
    time_int = time(timewindow);
    ntime_int = length(time_int);
    
    %hep = squeeze(mean(dat,1));
    %plot(time,hep) 
    
    %% copula transform and MI
    cdat = copnorm(dat(:,timewindow));
    crat = copnorm(ratings);
    
    MI = zeros(4,ntime_int);
    for r=1:4 % loop over 4 labels
        for t=1:ntime_int
            MI(r,t) = mi_gg(cdat(:,t),crat(:,r),true,true);
        end
    end
    
    %% plot MI
    subplot(4,1,1)
    plot(time_int,MI(1,:))
    title('valence')
    subplot(4,1,2)
    plot(time_int,MI(2,:))
    title('arousal')
    subplot(4,1,3)
    plot(time_int,MI(3,:))
    title('dominance')
    subplot(4,1,4)
    plot(time_int,MI(4,:))
    title('liking')
    
    %% cross-temporal interaction information    
    II = zeros(ntime_int,ntime_int);
    %for r=1:4
        for t1=1:ntime_int
            MI1(t1) = mi_gg(cdat(:,t1),crat(:,1),true,true); % I have put this here because otherwise MI1 and MI2 don't have the same size
            for t2=(t1+1):ntime_int
                MI2(t2) = mi_gg(cdat(:,t2),crat(:,1),true,true);
                JMI = mi_gg([cdat(:,t1) cdat(:,t2)],crat(:,1),true,true);
                II(t1,t2) = JMI - MI(1,t1) - MI(1,t2);     
            end
        end
        
        II = II + II';
        RED = min(MI1,MI2);
        U1 = MI1 - RED;
        U2 = MI2 - RED; % null matrix?
        
        for t1=1:ntime_int 
            for t2=(t1+1):ntime_int  
                SYN(t1,t2) = JMI - RED - U1(t1) - U2(t2);
            end
        end
    
    %end
    
    %% plot II
    plottitle = sprintf('Interaction information (bits) - sub %02.0f', isub);
    suptitle(plottitle)
    for r=1:4
        subplot(2,2,r);
        imagesc(time_int,time_int,II(:,:,r))
        colormap parula
        colorbar
        axis square
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        if r==1
            title('valence')
        elseif r==2
            title('arousal')
        elseif r==3
            title('dominance')
        else
            title('liking')
        end
    end
   
%end
