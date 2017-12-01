function [x_at_spiketimes, y_at_spiketimes]= visualize(file)
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
    x_at_spiketimes{i} = xN(find(spikes_binned(:,i)==1));
    y_at_spiketimes{i} = yN(find(spikes_binned(:,i)==1));
    
    % Plot
    subplot(2,5,i)
    plot(xN,yN,x_at_spiketimes{i},y_at_spiketimes{i},'r.');
    set(gca,'FontSize',14);
    axis tight square;
    xlabel('x position (m)','FontSize',18); ylabel('y position (m)','FontSize',18);
    title(['Spike positions, Neuron ' num2str(i)]);
    hold on;
end

end
