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
    
    %hep = squeeze(mean(dat,1));
    %plot(time,hep) 
    
    % time window we are interested in
    timewindow = (time>.2) & (time<.6);
    time_info = time(timewindow);
    ntime_info = length(time_info);
    
    %% copula transform and MI
    cdat = copnorm(dat(:,timewindow));
    cratings = copnorm(ratings);
    
    MI = zeros(4,ntime_info);
    for r=1:4 % loop over 4 labels
        for t=1:ntime_info
            MI(r,t) = mi_gg(cdat(:,t),cratings(:,r),true,true);
        end
    end
    
    %% plot MI
    subplot(4,1,1)
    plot(time_info,MI(1,:))
    title('valence')
    subplot(4,1,2)
    plot(time_info,MI(2,:))
    title('arousal')
    subplot(4,1,3)
    plot(time_info,MI(3,:))
    title('dominance')
    subplot(4,1,4)
    plot(time_info,MI(4,:))
    title('liking')
    
    %% cross-temporal interaction information 
    JMI = zeros(4,1);
    II = zeros(ntime_info,ntime_info,4);
    for r=1:4
        for t1=1:ntime_info
            for t2=(t1+1):ntime_info
                JMI(r,:) = mi_gg([cdat(:,t1) cdat(:,t2)],cratings(:,r),true,true);
                II(t1,t2,r) = JMI(r,:) - MI(r,t1) - MI(r,t2);     
            end
        end
        II(:,:,r) = II(:,:,r) + II(:,:,r)';
    end
    
    %% plot II
    plottitle = sprintf('Interaction information (bits) - sub %02.0f', isub);
    suptitle(plottitle)
    for r=1:4
        subplot(2,2,r);
        imagesc(time_info,time_info,II(:,:,r))
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
    
    %% partial information decomposition (only valence for now)
    % see Faes et al., 2017; Ince, 2017
    % II = JMI - MI1 - MI2;
    % RED = min(MI1,MI2);
    % U1 = MI1 - RED;
    % U2 = MI2 - RED;
    % SYN = JMI - RED - U1 - U2;
    
    II = zeros(ntime_info,ntime_info);
    for t1=1:ntime_info
        for t2=(t1+1):ntime_info
            JMI = mi_gg([cdat(:,t1) cdat(:,t2)],cratings(:,1),true,true);
            II(t1,t2) = JMI - MI(1,t1) - MI(1,t2);
            
            RED(t1,t2) = min(MI(1,t1),MI(1,t2));
            U1(t1) = MI(1,t1) - RED(t1,t2);
            U2(t2) = MI(1,t2) - RED(t1,t2);
            SYN(t1,t2) = JMI - RED(t1,t2) - U1(t1) - U2(t2);
        end
    end
    II = II + II';
    extrarow = zeros(1,50);
    SYN = [SYN ; extrarow]; % to have a square matrix
    SYN = SYN + SYN';
    RED = [RED ; extrarow];
    RED = RED + RED';
    
    diagonal = NaN(1,50);
    n = size(SYN,1);
    SYN(1:(n+1):end) = diagonal;
    RED(1:(n+1):end) = diagonal;
    
    %% plot PI
    plottitle = sprintf('Partial information - sub %02.0f', isub);
    suptitle(plottitle)
    subplot(1,2,1);
    imagesc(time_info,time_info,SYN(:,:))
    colormap parula
    colorbar
    axis square
    xlabel('Time (ms)')
    ylabel('Time (ms)')
    title('Synergy (bits)')
    subplot(1,2,2);
    imagesc(time_info,time_info,RED(:,:))
    colormap parula
    colorbar
    axis square
    xlabel('Time (ms)')
    ylabel('Time (ms)')
    title('Redundancy (bits)')
    
%end
