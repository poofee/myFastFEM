clear all
close all

fname=['D:\Staff\-毕业设计\myFastFEM\matlab\读取Mesh并绘制\model_normal.mphtxt'];

[xy,TR,DM] = readComsol(fname);

% 绘制单元
 for i=1:length(TR)
plot(xy(TR(i,[1 2 3 1]),1),xy(TR(i,[1 2 3 1]),2));
     hold on;
 end
  axis equal
%  分域绘制单元
figure
minDM = min(DM);
maxDM = max(DM);
for i=minDM:maxDM
    p = find(DM == i);
    for j=1:length(p)
        fill(xy(TR(p(j),:),1),xy(TR(p(j),:),2),[i/maxDM i*0/maxDM i*0/maxDM]);
        hold on
    end
end
 axis equal