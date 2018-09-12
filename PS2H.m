function [WS, PS_dz]  = PS2H( PS, ECG_srate )
b=fir1(48,[0.8 2.8]/(ECG_srate/2));
% PS = filtfilt(b,1,PS);
lambda_max = l1tf_lambdamax(PS);
[trend,~] = l1tf(PS, 0.0001*lambda_max);
PS_d=PS-trend;
PS_d = filtfilt(b,1,PS_d);
PS_dz=zscore(PS_d);
%PS_dz=PS;
%dPS=(diff(PS_dz,1));
%bdPS=dPS>0;
ddPS=(diff(PS_dz,2));
%ddPS(ddPS<0)=0;
%WS=bdPS(2:end).*ddPS;
WS=zscore(ddPS);
WS = zscore(filtfilt(b,1,WS));
% if min(WS)<0
%     WS=WS-min(WS);
% end
end

