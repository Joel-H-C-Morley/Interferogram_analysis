tic

Data12a=transpose(importdata('C:\Users\MOT lab\Documents\Data\12ms-1.csv'));
Data12=Data12a(2,:);
xData=importdata('C:\Users\MOT lab\Documents\Data\yData.csv')';
x = xData;
% trimx=linspace(-1e-3,1e-3,40);
% x=linspace(-5.3e-3,5.3e-3,916);

Res=length(x);

hold on
[pks,locs,widths,proms]=findpeaks(Data12,x,'MinPeakProminence',0.3);
plot(x,Data12)
plot(locs,pks,'+');

toc