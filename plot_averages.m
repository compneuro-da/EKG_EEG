clear;clc;close all
load DEAP_chnames

nlbls = 4;
time_start = -.2;
time_end = .6;
srate = 128;
time_epoch = time_start:1/srate:time_end-1/srate;
% time window we are interested in
timewindow = (time_epoch>.2) & (time_epoch<.6);
time_info = time_epoch(timewindow);
ntime_info = length(time_info);

for isub = 1:32
    load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_I.mat'])
    if isub == 1
       hep_tot = hep;
       MI_tot = MI;
       II_tot = II;
       RED_tot = RED;
       SYN_tot = SYN;
    else
       hep_tot = cat(1,hep_tot,hep);
       MI_tot = cat(3,MI_tot,MI);
       II_tot = cat(4,II_tot,II);
       RED_tot = cat(4,RED_tot,RED);
       SYN_tot = cat(4,SYN_tot,SYN);
    end
end

II_avg = squeeze(mean(II_tot,4));
RED_avg = squeeze(mean(RED_tot,4));
SYN_avg = squeeze(mean(SYN_tot,4));

%% plot average temporal II
plottitle = sprintf('Temporal interaction information (bits) -  all subjects');
suptitle(plottitle)
for l=1:nlbls
     subplot(2,2,l);
     imagesc(time_info,time_info,II_avg(:,:,l))
     colormap parula
     colorbar
     axis square
     xlabel('Time (ms)')
     ylabel('Time (ms)')
     title(plottitle)
     
     if l==1
         title('Valence')
     elseif l==2
         title('Arousal')
     elseif l==3
         title('Dominance')
     else
         title('Liking')
     end
end

%% plot average synergy
plottitle = sprintf('Temporal synergy (bits) - all subjects');
suptitle(plottitle)
for l=1:nlbls
    subplot(2,2,l);
    imagesc(time_info,time_info,SYN_avg(:,:,l))
    colormap parula
    colorbar
    axis square 
    xlabel('Time (ms)')
    ylabel('Time (ms)')
    title(plottitle)
    
    if l==1
       title('Valence')
    elseif l==2
       title('Arousal')
    elseif l==3
       title('Dominance')
    else
       title('Liking')
    end
end

%% plot average redundancy
plottitle = sprintf('Temporal redundancy (bits) - all subjects');
suptitle(plottitle)
for l=1:nlbls
    subplot(2,2,l);
    imagesc(time_info,time_info,RED_avg(:,:,l))
    colormap parula
    colorbar
    axis square 
    xlabel('Time (ms)')
    ylabel('Time (ms)')
    title(plottitle)
    
    if l==1
        title('Valence')
    elseif l==2
        title('Arousal')
    elseif l==3
        title('Dominance')
    else
        title('Liking')
    end
end
