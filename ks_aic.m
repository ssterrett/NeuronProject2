function [ks_stat, KSSorted, AIC] = ks_aic(file,b,dev,stats)
% Creates KS plot for data in file

% Error condition
if ~isa(file,'char')
    disp('Input argument must be in form ''[filename].mat''');
    return
end

load(file);

ks_stat = {};
KSSorted = {};
AIC = {};

%---------------- Model 1a, Covariates: X, Y, X^2, Y^2 --------------------

% Confidence intervals for each covariate
figure(1);
for i = 1:length(b{1}(1,:))
    subplot(2,5,i)
    errorbar(b{1}(:,i),stats{1}(i).se);
    hold on;
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 1, Neuron ' num2str(i)]);
    plot(1:length(b{1}(:,i)),zeros(length(b{1}(:,i)),1),'r:');
end

% P-values
for i = 1:length(b{1}(1,:))
    disp(['P-values for model 1, Neuron ' num2str(i) ': ' num2str(stats{1}(i).p')]);
end

% AIC
for i = 1:length(b{1}(1,:))
    aic1(i) = dev{1}(i) + 2*length(b{1}(:,i));
    disp(['AIC for Model 1, Neuron ' num2str(i) ': ' num2str(aic1(i))]);
end
    
% KS-plot
figure(101);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))
%     [b_lin(:,i),dev_lin(:,i),stats_lin(i)] = glmfit([xN yN xN.^2 yN.^2],spikes_binned(:,i),'poisson');
%     lambdaEst = exp(b_lin(1,i) + b_lin(2,i)*xN + b_lin(3,i)*yN + b_lin(4,i)*xN.^2 + b_lin(5,i)*yN.^2); 

    lambdaEst = exp(b{1}(1,i) + b{1}(2,i)*xN + b{1}(3,i)*yN + b{1}(4,i)*xN.^2 + b{1}(5,i)*yN.^2);
    
    timestep = 1;
    lambdaInt = 0;
    j=0;

    for t=1:length(spikes_binned(:,i))
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes_binned(t,i))
            j = j + 1;
            KS(j,i) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted1 = sort( KS(:,i) );
    N = length( KSSorted1);
    figure(101);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted1, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 1: Neuron ' num2str(i) ': KS Plot']);

    ks_stat1(i) = max(abs(KSSorted1' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS Stats: ' num2str(ks_stat1)]);

AIC{1} = aic1;
ks_stat{1} = ks_stat1;
KSSorted{1} = KSSorted1;

%---------------- Model 1b, Covariates: X, Y, X^2, Y^2, wall dist, speed, refractory, short, long --------------------
% clearvars -except b dev lambda lambdaAll stats AIC ks_stat KSSorted
% Computing speed
vN = zeros(length(xN),1);
vN(1) = 0;
for j = 2:(length(xN)-1)
    vN(j) = sqrt((xN(j)-xN(j-1)).^2 +(yN(j)-yN(j-1)).^2);
end
vN(length(xN)) = vN(length(xN)-1);

% Confidence intervals for each covariate
figure(2);
for i = 1:length(b{2}(1,:))
    subplot(2,5,i)
    e = errorbar(b{2}(:,i),stats{2}(i).se);
    e.Color = 'b';
    hold on;
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 1 (full), Neuron ' num2str(i)]);
    plot(1:length(b{2}(:,i)),zeros(length(b{2}(:,i)),1),'r:');
end

% P-values
for i = 1:length(b{2}(1,:))
    disp(['P-values for model 1 (full), Neuron ' num2str(i) ': ' num2str(stats{2}(i).p')]);
end

% AIC
for i = 1:length(b{2}(1,:))
    aic2(i) = dev{2}(i) + 2*length(b{2}(:,i));
    disp(['AIC for Model 1 (full), Neuron ' num2str(i) ': ' num2str(aic2(i))]);
end
    
% KS-plot
figure(102);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))
    % Defining history
    lambhr = zeros(length(spikes_binned(:,i))-260,length(spikes_binned(1:5,i)));
    for k = 1:length(spikes_binned(1:5,i))
        lambhr(:,k) = spikes_binned((261-k):(end-k),i);
    end
    
    lambhs = zeros(length(spikes_binned(:,i))-260,length(spikes_binned(80:110,i)));
    for k = 1:length(spikes_binned(80:110,i))
        lambhs(:,k) = spikes_binned((261-(k+79)):(end-(k+79)),i);
    end
    
    lambhl = zeros(length(spikes_binned(:,i))-260,length(spikes_binned(240:260,i)));
    for k = 1:length(spikes_binned(240:260,i))
        lambhl(:,k) = spikes_binned((261-(k+239)):(end-(k+239)),i);
    end

    % Defining lambdaEst
    % History terms
    lambhrBeta = 0;
    for k = 8:12
        lambhrBeta = lambhrBeta + b{2}(k,i).*lambhr(:,k-7);
    end
    
    lambhsBeta = 0;
    for k = 13:43
        lambhsBeta = lambhsBeta + b{2}(k,i).*lambhs(:,k-12);
    end
    
    lambhlBeta = 0;
    for k = 44:64
        lambhlBeta = lambhlBeta + b{2}(k,i).*lambhl(:,k-43);
    end
    
    % Splitting lambdaEst due to memory constraints
    lambdaEst = exp(b{2}(1,i) + b{2}(2,i)*xN(261:end) + b{2}(3,i)*yN(261:end) + b{2}(4,i)*xN(261:end).^2 + b{2}(5,i)*yN(261:end).^2 + b{2}(6,i)*abs(1-sqrt(xN(261:end).^2 + yN(261:end).^2)));
    lambdaEst2 = exp(b{2}(7,i)*vN(261:end) + lambhrBeta);
    lambdaEst3 = exp(lambhsBeta + lambhlBeta);
    
    timestep = 1;
    lambdaInt = 0;
    j=0;

    for t=1:length(spikes_binned(:,i))-261
        lambdaInt = lambdaInt + (lambdaEst(t).*lambdaEst2(t))*timestep;
        if (spikes_binned(t,i))
            j = j + 1;
            KS(j,i) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted2 = sort( KS(:,i) );
    N = length( KSSorted2);
    figure(102);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted2, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 1: Neuron ' num2str(i) ': KS Plot']);

    ks_stat2(i) = max(abs(KSSorted2' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS Stats: ' num2str(ks_stat2)]);

AIC{2} = aic2;
ks_stat{2} = ks_stat2;
KSSorted{2} = KSSorted2;




%---------------- Model 2, Covariates: Sinusoids -------------------- 
% Computing speed
vN = zeros(length(xN),1);
vN(1) = 0;
for j = 2:(length(xN)-1)
    vN(j) = sqrt((xN(j)-xN(j-1)).^2 +(yN(j)-yN(j-1)).^2);
end
vN(length(xN)) = vN(length(xN)-1);

% Confidence intervals for each covariate
figure(3);
for i = 1:length(b{3}(1,:))
    subplot(2,5,i)
    errorbar(b{3}(:,i),stats{3}(i).se);
    hold on;
    plot(1:length(b{3}(:,i)),zeros(length(b{3}(:,i)),1),'r:');
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 2, Neuron ' num2str(i)]);
end

% P-values
for i = 1:length(b{3}(1,:))
    disp(['P-values for model 2, Neuron ' num2str(i) ': ' num2str(stats{3}(i).p')]);
end

% AIC
for i = 1:length(b{3}(1,:))
    aic3(i) = dev{3}(i) + 2*length(b{3}(:,i));
    disp(['AIC for Model 2, Neuron ' num2str(i) ': ' num2str(aic3(i))]);
end

% KS-plot
figure(103);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))

    lambdaEst = exp(b{3}(1,i) + b{3}(2,i)*sin(2.2*pi*(yN+0.8*xN)) + b{3}(3,i)*sin(3*pi*(yN-0.35*xN)-9) + b{3}(4,i)*vN);
    
    timestep = 1;
    lambdaInt = 0;
    j=0;

    for t=1:length(spikes_binned(:,i))
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes_binned(t,i))
            j = j + 1;
            KS(j,i) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted3 = sort( KS(:,i) );
    N = length( KSSorted3);
    figure(103);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted3, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 2: Neuron ' num2str(i) ': KS Plot']);

    ks_stat3(i) = max(abs(KSSorted3' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS Stats: ' num2str(ks_stat3)]);

AIC{3} = aic3;
ks_stat{3} = ks_stat3;
KSSorted{3} = KSSorted3;



end



