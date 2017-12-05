function isi(datafile)
close all

if ~isa(datafile,'char')
    disp('Input argument must be in form ''[filename].mat''');
    return
end

load(datafile); 
[rows neurons] = size(spikes_binned);
spike_times = cell(neurons);
isi = cell(neurons);

for i = 1:neurons;
spike_times{i} = find(spikes_binned(:,i) == 1);
isi{i} = diff(spike_times{i});

if i==3
    % Separate figures
    % figure(i); 

    % One figure with 10 subplots
    subplot(2,5,i);

    histogram(isi{i},15000);
    set(gca,'FontSize',14);
    xlim([0 300]);
    xlabel('ISI (ms)','FontSize',15);
    ylabel('Count','FontSize',15);
    title(['ISI Histogram for Neuron ' num2str(i)]);
elseif i==4
    % Separate figures
    % figure(i); 

    % One figure with 10 subplots
    subplot(2,5,i);

    histogram(isi{i},8000);
    set(gca,'FontSize',14);
    xlim([0 300]);
    xlabel('ISI (ms)','FontSize',15);
    ylabel('Count','FontSize',15);
    title(['ISI Histogram for Neuron ' num2str(i)]);
elseif i==6
    % Separate figures
    % figure(i); 

    % One figure with 10 subplots
    subplot(2,5,i);

    histogram(isi{i}, 6000);
    set(gca,'FontSize',14);
    xlim([0 300]);
    xlabel('ISI (ms)','FontSize',15);
    ylabel('Count','FontSize',15);
    title(['ISI Histogram for Neuron ' num2str(i)]);
elseif i==7
    % Separate figures
    % figure(i); 

    % One figure with 10 subplots
    subplot(2,5,i);

    histogram(isi{i},15000);
    set(gca,'FontSize',14);
    xlim([0 300]);
    xlabel('ISI (ms)','FontSize',15);
    ylabel('Count','FontSize',15);
    title(['ISI Histogram for Neuron ' num2str(i)]);
else
    % Separate figures
    % figure(i); 

    % One figure with 10 subplots
    subplot(2,5,i);

    histogram(isi{i},10000);
    set(gca,'FontSize',14);
    xlim([0 300]);
    xlabel('ISI (ms)','FontSize',15);
    ylabel('Count','FontSize',15);
    title(['ISI Histogram for Neuron ' num2str(i)]);
end

end
