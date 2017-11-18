function [x_at_spiketimes y_at_spiketimes] = spatialplots(dataset)
load(dataset);
[rows neurons] = size(spikes_binned);
x_at_spiketimes = cell(neurons);
y_at_spiketimes = cell(neurons);
for i = 1:10
x_at_spiketimes{i} = xN(find(spikes_binned(:,i) == 1));
y_at_spiketimes{i} = yN(find(spikes_binned(:,i) == 1));
figure(i);
plot(xN,yN,x_at_spiketimes{i}, y_at_spiketimes{i},'r.');
end
end
