clc; clear all;
%first, run analysis on target interval
tempoLabels = ["fast","medium","slow"];
targetMean = [0.4 0.6 1.0];
targetSD = targetMean.*0.05;


% actualMean = [0.402, 0.570, 0.933];
% bias = actualMean - targetMean;
% actualSD = [0.014, 0.02, 0.037];


% plot(actualMean,actualSD,'ko'); hold on;
% plot(targetMean,targetSD,'ro');
% xlim([0.3 1.1]);
% ylim([0 0.05]); hold off;

rep = 1728;

%create Gaussian distribution and draw from this distribution
for i = 1:3
    gaussian_noisy(i,:) = normrnd(targetMean(i),targetSD(i),rep,1);
    rt_adjustment(i,:) = [gaussian_noisy(i,:) - [gaussian_noisy(i,2:rep) 0]];
    rt_adjustment(i,2:rep) = rt_adjustment(i,1:(rep-1)) + normrnd(0,targetSD(i),1,(rep-1)); %add in Gaussian noise
    rt_adjustment(i,1) = NaN;
end

for i = 1:3
    for j = 10:rep
        running_mean(i,j) = mean(gaussian_noisy(i,(j-9):j));
        running_mean(i,1:9) = NaN;
    end
end

for i = 1:3
 histogram(running_mean(i,:)); hold on;
end


deviation_from_the_mean = gaussian_noisy - running_mean;

%%
% for i = 1:3
%  subplot(1,3,i);
%  histogram(deviation_from_the_mean(i,:));
%  xlim([-0.2 0.2]);
% end

for i = 1:3
 subplot(1,3,i);
 histogram(rt_adjustment(i,:));
end

%%
for i = 1:3
 plot(deviation_from_the_mean(i,:),rt_adjustment(i,:),'o'); 
 lsline; hold on;
 xlim([-0.2 0.2]);
 ylim([-0.2 0.2]);
end
