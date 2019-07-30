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
preclust_pval = 0.05;
clust_pval = 0.05;

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

%% build null hypothesis maps
h0 = zeros(ntime_info,ntime_info);
hp = zeros(ntime_info,ntime_info,nperm);
clust_max = zeros(nperm,1);
for t1=1:ntime_info
    for t2=(t1+1):ntime_info
        JMI_tot = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,1),false);
        II_tot = JMI_tot - mi_gg(cdat_tot(:,t1),cratings_tot(:,1),false) - mi_gg(cdat_tot(:,t2),cratings_tot(:,1),false);
        
        h0(t1,t2) = II_tot;
        n = length(cratings_tot);
        
        for iperm = 1:nperm
            cratings_perm(:,1) = cratings_tot(randperm(n));
            JMI_perm = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_perm(:,1),false);
            II_perm = JMI_perm - mi_gg(cdat_tot(:,t1),cratings_perm(:,1),false) - mi_gg(cdat_tot(:,t2),cratings_perm(:,1),false);
            
            hp(t1,t2,iperm) = II_perm;
        end
    end
end

% fill in the other half of the matrices
h0 = h0 + h0';
for iperm = 1:nperm 
    hp(:,:,iperm) = hp(:,:,iperm) + hp(:,:,iperm)'; 
end 

%% collect the largest suprathreshold clusters
for iperm = 1:nperm
    perms = true(1,nperm);
    perms(iperm) = 0;
    zvals = squeeze((hp(:,:,iperm) - mean(hp(:,:,perms),3)) ./ std(hp(:,:,perms),[],3)); % see Cohen's book chapter 33 equation 33.1
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0; 
    
    clust_info = bwconncomp(zvals);
    clust_max(iperm) = max([0 cellfun(@numel,clust_info.PixelIdxList)]);
end

%% identify significant clusters in real data
zmap = squeeze((h0 - mean(hp,3)) ./ std(hp,[],3));
zmap(abs(zmap) < norminv(1 - preclust_pval)) = 0;

clust_info = bwconncomp(zmap);
clust_size = cellfun(@numel,clust_info.PixelIdxList);
clust_th = prctile(clust_max,100-clust_pval*100);
clust_rem = find(clust_size < clust_th);
for i=1:length(clust_rem)
    zmap(clust_info.PixelIdxList{clust_rem(i)})=0;
end
zmap(isnan(zmap)) = 0;
zmap = logical(zmap);

%% plot
time_plot = time_info*1000;
imagesc(time_plot,time_plot,h0(:,:))
hold on
contour(time_plot,time_plot,zmap,1,'linecolor','k','LineWidth',1)
colormap(brewermap([],'*RdBu'));
colorbar
axis square
lim = max(abs(min(min(h0(:,:)))),max(max(h0(:,:))));
caxis([-lim,lim]) % center zero
xlabel('Time (ms)'),ylabel('Time (ms)')
