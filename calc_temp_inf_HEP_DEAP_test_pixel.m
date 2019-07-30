clear;clc;close all
load DEAP_chnames
chan = 24; % choose a channel
l = 1; % choose a label
time_start = -.2;
time_end = .6;
srate = 128;
time_epoch = time_start:1/srate:time_end-1/srate;
ntime_epoch = length(time_epoch);
% time window we are interested in
timewindow = (time_epoch>.2) & (time_epoch<.6);
time_info = time_epoch(timewindow);
ntime_info = length(time_info);

nperm = 5000; % number of permutations
pval = 0.05; % uncorrected threshold
%pcorr = pval*6/(ntime_info*(ntime_info-1)*(ntime_info-2)); % Bonferroni correction for triplets (from https://github.com/compneuro-da/interactioninformation/blob/master/main_interaction_information_pwp.m)

%% load the data
for isub = 1:32
    load(['D:\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP.mat'])
    if isub == 1
       avg_HEP_tot = avg_HEP;
       labels_tot = labels;
    else
       avg_HEP_tot = cat(1,avg_HEP_tot,avg_HEP);
       labels_tot = cat(1,labels_tot,labels);
    end
end
dat_tot = avg_HEP_tot(:,timewindow,chan);
[ntrls_tot,~] = size(dat_tot);

cdat_tot = copnorm(dat_tot);
cratings_tot = copnorm(labels_tot(:,l));
%%
h0 = zeros(ntime_info,ntime_info);
h = zeros(ntime_info,ntime_info);
for t1=1:ntime_info
    for t2=(t1+1):ntime_info
        JMI_tot = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,1),false);
        II_tot = JMI_tot - mi_gg(cdat_tot(:,t1),cratings_tot(:,1),false) - mi_gg(cdat_tot(:,t2),cratings_tot(:,1),false);
        
        h0(t1,t2) = II_tot;
        n = length(cratings_tot);
        ncont = 0;
        
        if h0(t1,t2) == 0
           h(t1,t2) = 0;
        else
           for i = 1:nperm
               cratings_perm(:,1) = cratings_tot(randperm(n));
               JMI_perm = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_perm(:,1),false);
               II_perm = JMI_perm - mi_gg(cdat_tot(:,t1),cratings_perm(:,1),false) - mi_gg(cdat_tot(:,t2),cratings_perm(:,1),false);
               
               hp = II_perm;
               
               if h0(t1,t2)>0 && hp>h0(t1,t2)
                  ncont=ncont+1;
               elseif h0(t1,t2)<0 && hp<h0(t1,t2)
                  ncont=ncont+1;
               end
           end
           if ncont>pval*nperm % p = ncont/nperm
              h(t1,t2) = 0;
           else
              h(t1,t2) = h0(t1,t2);
           end
        end
    end
end
h0 = h0 + h0';
h = h + h';
h = logical(h);

%% plot
time_plot = time_info*1000;
imagesc(time_plot,time_plot,h0(:,:))
hold on
contour(time_plot,time_plot,h,1,'linecolor','k','LineWidth',1)
colormap(brewermap([],'*RdBu'));
colorbar
axis square
lim = max(abs(min(min(h0(:,:)))),max(max(h0(:,:))));
caxis([-lim,lim]) % center zero
xlabel('Time (ms)'),ylabel('Time (ms)')
