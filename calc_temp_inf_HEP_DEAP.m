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
    nlbls = size(ratings,2);
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
    
    %% copula transform and calculate mutual information (MI) at each time point
    cdat = copnorm(dat(:,timewindow));
    cratings = copnorm(ratings);
    
    MI = zeros(nlbls,ntime_info);
    for l=1:nlbls % loop over 4 labels
        for t=1:ntime_info
            MI(l,t) = mi_gg(cdat(:,t),cratings(:,l),true,true);
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
    
    %% plot II
    plottitle = sprintf('Interaction information (bits) - sub %02.0f', isub);
    suptitle(plottitle)
    for l=1:nlbls
        subplot(2,2,l);
        imagesc(time_info,time_info,II(:,:,l))
        colormap parula
        colorbar
        axis square
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        if l==1
            title('valence')
        elseif l==2
            title('arousal')
        elseif l==3
            title('dominance')
        else
            title('liking')
        end
    end
    
    %% partial information decomposition (PID)
    % see Faes et al., 2017; Ince, 2017
    % II = JMI - MI1 - MI2;
    % RED = min(MI1,MI2);
    % U1 = MI1 - RED;
    % U2 = MI2 - RED;
    % SYN = JMI - RED - U1 - U2;
    
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
    
    %% plot PI
    SYN(SYN==0)=NaN;
    RED(RED==0)=NaN;
    
    plottitle = sprintf('Partial information (synergy) - sub %02.0f', isub);
    suptitle(plottitle)
    for l=1:nlbls
        subplot(2,2,l);
        imagesc(time_info,time_info,SYN(:,:,l))
        colormap parula
        colorbar
        axis square 
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        title('Synergy (bits)')
        if l==1
            title('valence')
        elseif l==2
            title('arousal')
        elseif l==3
            title('dominance')
        else
            title('liking')
        end
    end
    %%
    plottitle = sprintf('Partial information (redundancy) - sub %02.0f', isub);
    suptitle(plottitle)
    for l=1:nlbls
        subplot(2,2,l);
        imagesc(time_info,time_info,RED(:,:,l))
        colormap parula
        colorbar
        axis square 
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        title('Synergy (bits)')
        if l==1
            title('valence')
        elseif l==2
            title('arousal')
        elseif l==3
            title('dominance')
        else
            title('liking')
        end
    end
    
%end
