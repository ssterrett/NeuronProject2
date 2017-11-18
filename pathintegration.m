function pathintegration(dataset)
load(dataset);
close all;
steps = 50;
for k = 1:10;
    
ind = find(spikes_binned(:,k) == 1);
figure(k); hold on;
for i = 1:length(ind)
    plot(xN(ind(i)-steps:ind(i))-xN(ind(i)),yN(ind(i)-steps:ind(i))-yN(ind(i)));
end
end