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
%ntr = 50; % number of triplets
%pcorr = pval*6/(ntr*(ntr-1)*(ntr-2)); % Bonferroni correction for triplets

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

%%
cdat_tot = copnorm(dat_tot);
cratings_tot = copnorm(labels_tot(:,3));

MI_tot = zeros(1,ntime_info);
for t=1:ntime_info
    MI_tot(t) = mi_gg(cdat_tot(:,t),cratings_tot(:,1),false);
end

for t1=1:ntime_info
     for t2=(t1+1):ntime_info
         JMI_tot = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,1),false);
         II_tot = JMI_tot - MI_tot(t1) - MI_tot(t2);
         
         h0 = II_tot;
         n = length(cratings_tot);
         ncont = 0;
         
         if h0 == 0
            h(t1,t2) = 0;
         else
             for i = 1:nperm
                 cratings_tot(:,1) = cratings_tot(randperm(n));
                 JMI_tot = mi_gg([cdat_tot(:,t1) cdat_tot(:,t2)],cratings_tot(:,1),false);
                 II_tot = JMI_tot - mi_gg(cdat_tot(:,t1),cratings_tot(:,1),false) - mi_gg(cdat_tot(:,t2),cratings_tot(:,1),false);
                 hp = II_tot;
                 if hp>h0
                    ncont=ncont+1;
                 elseif hp<h0
                    ncont=ncont+1;
                 end
             end
             if ncont>pval*nperm
                h(t1,t2) = 0;
             else
                h(t1,t2) = h0;
             end
         end
     end
end
