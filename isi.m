function isi(datafile)

load(datafile); 
[rows neurons] = size(spikes_binned);
spike_times = cell(neurons);
isi = cell(neurons);
for i = 1:neurons;
spike_times{i} = find(spikes_binned(:,i) == 1);
isi{i} = diff(spike_times{i});
figure(i); 
histogram(isi{i},0:10:10000);
xlim([0 1000]);
xlabel('ISI (ms)');
ylabel('Count');
title(['ISI Histogram for Neuron ' num2str(i)]);
end
