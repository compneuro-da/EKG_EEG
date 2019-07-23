clear;clc;close all
load DEAP_chnames
nchan = 3; % number of electrodes
chans = [24 19 16]; % choose electrodes of interest (Cz, Fz, Pz)
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
    load(['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP.mat']) % load data, change directory accordingly
    
    % plot example of HEP for one trial
    if isub == 1
       title('Heartbeat-evoked potential (one trial)')
       plot(time_epoch*1000,avg_HEP(40,:,24),'-o');
       caxis([-0.0025,0.0025])
       xlabel('Time (ms)')
    end
    
    if isub == 1
       avg_HEP_tot = avg_HEP;
       labels_tot = labels;
    else
       avg_HEP_tot = cat(1,avg_HEP_tot,avg_HEP);
       labels_tot = cat(1,labels_tot,labels);
    end
    
    for ichan = 1:nchan
        dat = avg_HEP(:,timewindow,chans(ichan));
        [ntrls,~] = size(dat);
        
        HEP = squeeze(mean(dat,1));
        %plot(time_info,HEP)
        
        %% copula transform and calculate mutual information (MI) at each time point
        cdat = copnorm(dat);
        cratings = copnorm(labels);
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
        
        %% save results
        savefile = ['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_I_' char(ch_labels(chans(ichan))) '.mat'];
        save(savefile,'HEP','MI','II','SYN','RED');
    end
end


%% calculate temporal information over all subjects
for ichan = 1:nchan
    dat_tot = avg_HEP_tot(:,timewindow,chans(ichan));
    [ntrls_tot,~] = size(dat_tot);
    
    HEP_tot = squeeze(mean(dat_tot,1));
    %plot(time_info,HEP_tot)
    
    % copula transform and calculate mutual information (MI) at each time point
    cdat_tot = copnorm(dat_tot);
    cratings_tot = copnorm(labels_tot);
    
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
    
    %% save results
    savefile = ['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\tot_I_' char(ch_labels(chans(ichan))) '.mat'];
    save(savefile,'HEP_tot','MI_tot','II_tot','SYN_tot','RED_tot');
end


%% plot the temporal information quantities
time_plot = time_info*1000;
for ichan = 1%:nchan
    load(['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\tot_I_' char(ch_labels(chans(ichan))) '.mat'])
    
    %% plot temporal II
    plottitle = sprintf('Temporal interaction information (bits)');
    suptitle(plottitle)
    for l=1:nlbls
        subplot(2,2,l);
        imagesc(time_plot,time_plot,II_tot(:,:,l))
        colormap(brewermap([],'*RdBu'));
        colorbar
        axis square
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        
        %max(max(II_tot(:,:,l))) % check max values
        %min(min(II_tot(:,:,l))) % check min values
        if l == 1
           title('Valence')
           if ichan == 1
              caxis([-0.0020,0.0020])
           elseif ichan == 2
              caxis([-0.0015,0.0015])
           else
              caxis([-0.0025,0.0025])
           end
        elseif l == 2
           title('Arousal')
           if ichan == 1
              caxis([-0.0035,0.0035])
           elseif ichan == 2
              caxis([-0.0030,0.0030])
           else
              caxis([-0.0020,0.0020])
           end
        elseif l == 3
           title('Dominance')
           if ichan == 1
              caxis([-0.0035,0.0035])
           elseif ichan == 2
              caxis([-0.0045,0.0045])
           else
              caxis([-0.0065,0.0065])
           end
        else
           title('Liking')
           if ichan == 1
              caxis([-0.0025,0.0025])
           elseif ichan == 2
              caxis([-0.0030,0.0030])
           else
              caxis([-0.0025,0.0025])
           end
        end
    end
    
    %% plot temporal synergy
    plottitle = sprintf('Temporal synergy (bits)');
    suptitle(plottitle)
    for l=1:nlbls
        subplot(2,2,l);
        imagesc(time_plot,time_plot,SYN_tot(:,:,l))
        colormap(brewermap([],'Reds'));
        colorbar
        axis square 
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        
        if l==1
           title('Valence')
           caxis([0,max(max(SYN_tot(:,:,l)))])
        elseif l==2
           title('Arousal')
           caxis([0,max(max(SYN_tot(:,:,l)))])
        elseif l==3
           title('Dominance')
           caxis([0,max(max(SYN_tot(:,:,l)))])
        else
           title('Liking')
           caxis([0,max(max(SYN_tot(:,:,l)))])
        end
    end
    
    %% plot temporal redundancy
    plottitle = sprintf('Temporal redundancy (bits)');
    suptitle(plottitle)
    for l=1:nlbls
        subplot(2,2,l);
        imagesc(time_plot,time_plot,RED_tot(:,:,l))
        colormap(brewermap([],'Blues'));
        colorbar
        axis square 
        xlabel('Time (ms)')
        ylabel('Time (ms)')
        
        if l==1
            title('Valence')
            caxis([0,max(max(RED_tot(:,:,l)))])
        elseif l==2
            title('Arousal')
            caxis([0,max(max(RED_tot(:,:,l)))])
        elseif l==3
            title('Dominance')
            caxis([0,max(max(RED_tot(:,:,l)))])
        else
            title('Liking')
            caxis([0,max(max(RED_tot(:,:,l)))])
        end
    end
end
