clear;clc;close all
load DEAP_chnames
isub = 1;
load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP_I.mat']) % load data, change directory accordingly

chan = 24; % choose an electrode
nlbls = 4;
time_start = -.2;
time_end = .6;
srate = 128;
time_epoch = time_start:1/srate:time_end-1/srate;
ntime_epoch = length(time_epoch);
% time window we are interested in
timewindow = (time_epoch>.2) & (time_epoch<.6);
time_info = time_epoch(timewindow);
ntime_info = length(time_info);

%% plot MI at each time point
suptitle('Mutual information (MI)');
subplot(4,1,1)
plot(time_info,MI(1,:))
title('Valence')
subplot(4,1,2)
plot(time_info,MI(2,:))
title('Arousal')
subplot(4,1,3)
plot(time_info,MI(3,:))
title('Dominance')
subplot(4,1,4)
plot(time_info,MI(4,:))
title('Liking')

%% plot temporal II
plottitle = sprintf('Temporal interaction information (bits) - sub %02.0f', isub);
suptitle(plottitle)
for l=1:nlbls
     subplot(2,2,l);
     imagesc(time_info,time_info,II(:,:,l))
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

%% plot MI (or HEP) together with previous plot (for each label individually)
l = 1; % choose label
xl = [.2 .6];
axm = subplot(5,5,[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
plottitle = sprintf('Temporal interaction information (bits) - sub %02.0f', isub);
suptitle(plottitle)
imagesc(time_info,time_info,II(:,:,l))
colormap parula
colorbar
axis square
xlabel('Time (ms)')
ylabel('Time (ms)')
if l==1
   title('Valence')
elseif l==2
   title('Arousal')
elseif l==3
   title('Dominance')
else
   title('Liking')
end
lw = 2;
subplot(5,5,[22 23 24 25])
plot(time_info,MI(l,:),'k','LineWidth',lw);
axis tight
xlim(xl)
intpos = get(axm,'Pos');
pos = get(gca,'Pos');
pos(1) = intpos(1);
pos(3) = intpos(3);
set(gca,'Pos',pos)
box off
subplot(5,5,[1 6 11 16])
plot(time_info,MI(l,:),'k','LineWidth',lw);
axis tight
box off
xlim(xl)
set(gca,'CameraUpVector',[-1 0 0])

%% plot temporal synergy
plottitle = sprintf('Temporal synergy (bits) - sub %02.0f', isub);
suptitle(plottitle)
for l=1:nlbls
    subplot(2,2,l);
    imagesc(time_info,time_info,SYN(:,:,l))
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
%%
l = 1; % choose label
xl = [.2 .6];
axm = subplot(5,5,[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
plottitle = sprintf('Temporal synergy (bits) - sub %02.0f', isub);
suptitle(plottitle)
imagesc(time_info,time_info,SYN(:,:,l))
colormap parula
colorbar
axis square
xlabel('Time (ms)')
ylabel('Time (ms)')
if l==1
   title('Valence')
elseif l==2
   title('Arousal')
elseif l==3
   title('Dominance')
else
   title('Liking')
end
lw = 2;
subplot(5,5,[22 23 24 25])
plot(time_info,MI(l,:),'k','LineWidth',lw);
axis tight
xlim(xl)
intpos = get(axm,'Pos');
pos = get(gca,'Pos');
pos(1) = intpos(1);
pos(3) = intpos(3);
set(gca,'Pos',pos)
box off
subplot(5,5,[1 6 11 16])
plot(time_info,MI(l,:),'k','LineWidth',lw);
axis tight
box off
xlim(xl)
set(gca,'CameraUpVector',[-1 0 0])

%% plot temporal redundancy
plottitle = sprintf('Temporal redundancy (bits) - sub %02.0f', isub);
suptitle(plottitle)
for l=1:nlbls
    subplot(2,2,l);
    imagesc(time_info,time_info,RED(:,:,l))
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
%%
l = 1; % choose label
xl = [.2 .6];
axm = subplot(5,5,[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
plottitle = sprintf('Temporal redundancy (bits) - sub %02.0f', isub);
suptitle(plottitle)
imagesc(time_info,time_info,RED(:,:,l))
colormap parula
colorbar
axis square
xlabel('Time (ms)')
ylabel('Time (ms)')
if l==1
   title('Valence')
elseif l==2
   title('Arousal')
elseif l==3
   title('Dominance')
else
   title('Liking')
end
lw = 2;
subplot(5,5,[22 23 24 25])
plot(time_info,MI(l,:),'k','LineWidth',lw);
axis tight
xlim(xl)
intpos = get(axm,'Pos');
pos = get(gca,'Pos');
pos(1) = intpos(1);
pos(3) = intpos(3);
set(gca,'Pos',pos)
box off
subplot(5,5,[1 6 11 16])
plot(time_info,MI(l,:),'k','LineWidth',lw);
axis tight
box off
xlim(xl)
set(gca,'CameraUpVector',[-1 0 0])