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

nperm = 1000; % number of permutations
pval = 0.05;
%pcorr = pval*6/(3*(3-1)*(3-2)); % Bonferroni correction for triplets

%% load the data
for isub = 1:32
    load(['\\client\d$\Users\Liesa\Documents\Universiteit Gent\Theoretische en experimentele psychologie\MA05\05 J\5 Masterproef II\DEAP\preprocessed\s' num2str(isub,'%02.0f') '_avgHEP.mat'])
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

%% copula transform and calculate mutual information (MI) at each time point
cdat_tot = copnorm(dat_tot);
cratings_tot = copnorm(labels_tot(:,1));

MI_tot = zeros(1,ntime_info);
for t=1:ntime_info
    MI_tot(t) = mi_gg(cdat_tot(:,t),cratings_tot(:,1),false);
end

% cross-temporal interaction information (II)
II_tot = zeros(ntime_info,ntime_info);
for t1=1:ntime_info
     for t2=(t1+1):ntime_info
         JMI_tot = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,1),false);
         II_tot(t1,t2) = JMI_tot - MI_tot(t1) - MI_tot(t2);
     end
end
II_tot(:,:) = II_tot(:,:) + II_tot(:,:)';
II_tot(II_tot==0)=NaN;

%% test 
h0 = II_tot;
n = length(cratings_tot);
ncont = 0;
h = zeros(ntime_info,ntime_info);

if h0 == 0
   h = 0;
else
   for i = 1:nperm
       cratings_tot(:,1) = cratings_tot(randperm(n));
       MI_tot = zeros(1,ntime_info);
       for t=1:ntime_info
           MI_tot(t) = mi_gg(cdat_tot(:,t),cratings_tot(:,1),false);
       end
       II_tot = zeros(ntime_info,ntime_info);
       for t1=1:ntime_info
           for t2=(t1+1):ntime_info
               JMI_tot = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,1),false);
               II_tot(t1,t2) = JMI_tot - MI_tot(t1) - MI_tot(t2);
           end
       end
       II_tot(:,:) = II_tot(:,:) + II_tot(:,:)';
       II_tot(II_tot==0)=NaN;
       
       h = II_tot;
       if h>h0
          ncont = ncont+1;
       end
       if h<h0
          ncont = ncont+1;
       end
   end
   if ncont>pval*nperm
      h = 0;
   else
      h = h0;
   end
end
