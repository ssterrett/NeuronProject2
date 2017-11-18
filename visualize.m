function visualize(file)
% Visualizes spike positions over time

% Error condition
if ~isa(file,'char')
    disp('Input argument must be in form ''[filename].mat''');
    return
end


load(file);

figure;
for i = 1:length(spikes_binned(1,:))
    % Define x,y position where spikes occur
    x_at_spiketimes = xN(find(spikes_binned(:,i)==1));
    y_at_spiketimes = yN(find(spikes_binned(:,i)==1));
    % NEED TO SAVE X/Y AT SPIKETIMES
    
    % Plot
    subplot(2,5,i)
    plot(xN,yN,x_at_spiketimes,y_at_spiketimes,'r.');
    axis tight square;
    xlabel('x position (m)'); ylabel('y position (m)');
    title(['Spike positions, Neuron ' num2str(i)]);
    hold on;
end

end
