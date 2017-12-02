function pathintegration(dataset)
load(dataset);
close all;
steps = 25;
for k = 1:10;
    
ind = find(spikes_binned(:,k) == 1);
subplot(2,5,k); hold on;
% axis square;
store = zeros(length(ind),1);
for i = 1:length(ind)
    pos = [(xN(ind(i)-steps) - xN(ind(i))) (yN(ind(i) - steps) - yN(ind(i))) 0];
    store(i) = angle((xN(ind(i)-steps) - xN(ind(i))) -1j*(yN(ind(i) - steps) - yN(ind(i))));
%     plot([xN(ind(i)-steps)-xN(ind(i)) 0]./(sqrt((xN(ind(i)-steps)-xN(ind(i))).^2 + (yN(ind(i)-steps)-yN(ind(i))).^2)),[yN(ind(i)-steps)-yN(ind(i)) 0]./(sqrt((xN(ind(i)-steps)-xN(ind(i))).^2 + (yN(ind(i)-steps)-yN(ind(i))).^2)));
end
histogram(store);
xticks([-pi 0 pi]);
xticklabels({'-\pi' '0' '\pi'});
ylabel('Spikes'); xlabel('Head Direction');
title(['Neuron ' num2str(k)]);

xNd = xN(50:end,k) - xN(1:(end-50),k);
yNd = yN(50:end,k) - yN(1:(end-50),k);

end