function [r,a,cpa]=sampling(pa,ns)

a=zeros(1,ns);

sample=mnrnd(ns,pa,1); %model abundance, aim too high for numbers so can get right number after removing zeros
sample=sample(sample>0); %don't know which ones were missed

a(1:length(sample))=-sort(-sample)/sum(sample); %make sure sorted, and normalize
r=1:length(a);

cpa=cumsum(a);

end