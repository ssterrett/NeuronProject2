Fs = 1000;
n = 1:length(spikes_binned(:,1));
for F = 6:10;
    w = 2*pi*F/Fs;
    figure(F);
    title(['Rhythm = ' num2str(F) ' Hz']);

for i = 1:10
    phase = asin(sin(w.*n));
    phase_spikes = phase(find(spikes_binned(:,i) == 1));
    subplot(2,5,i);
    histogram(phase_spikes);
    title(['Neuron ' num2str(i)]);
end
end
