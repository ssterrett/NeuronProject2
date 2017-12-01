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

%---------------- Model 1, Covariates: X, Y, X^2, Y^2 --------------------

% Confidence intervals for each covariate
figure(1);
for i = 1:length(b{1}(1,:))
    subplot(2,5,i)
    errorbar(b{1}(:,i),stats{1}(i).se);
    hold on;
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 1, Neuron ' num2str(i)]);
    plot(1:length(b{1}(:,i)),zeros(length(b{1}(:,i)),1),'r-');
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

%---------------- Model 2, Covariates: Sinusoids -------------------- 

% Confidence intervals for each covariate
figure(2);
for i = 1:length(b{2}(1,:))
    subplot(2,5,i)
    errorbar(b{2}(:,i),stats{2}(i).se);
    hold on;
    plot(1:length(b{2}(:,i)),zeros(length(b{2}(:,i)),1),'r-');
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 2, Neuron ' num2str(i)]);
end

% P-values
for i = 1:length(b{2}(1,:))
    disp(['P-values for model 2, Neuron ' num2str(i) ': ' num2str(stats{2}(i).p')]);
end

% AIC
for i = 1:length(b{2}(1,:))
    aic2(i) = dev{2}(i) + 2*length(b{2}(:,i));
    disp(['AIC for Model 2, Neuron ' num2str(i) ': ' num2str(aic2(i))]);
end

% KS-plot
figure(102);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))

    lambdaEst = exp(b{2}(1,i) + b{2}(2,i)*sin(2.5*pi*xN) + b{2}(3,i)*sin(2.5*pi*yN));
    
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
    figure(102);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted1, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 2: Neuron ' num2str(i) ': KS Plot']);

    ks_stat1(i) = max(abs(KSSorted1' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS Stats: ' num2str(ks_stat1)]);

AIC{2} = aic2;
ks_stat{2} = ks_stat1;
KSSorted{2} = KSSorted1;

%---------------- Model 3, Covariates: Refractory period history-dependence ---------------
% Confidence intervals for each covariate
figure(3);
for i = 1:length(b{3}(1,:))
    subplot(2,5,i)
    errorbar(b{3}(:,i),stats{3}(i).se);
    hold on;
    plot(1:length(b{3}(:,i)),zeros(length(b{3}(:,i)),1),'r-');
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 3, Neuron ' num2str(i)]);
end

% P-values
for i = 1:length(b{3}(1,:))
    disp(['P-values for model 3, Neuron ' num2str(i) ': ' num2str(stats{3}(i).p')]);
end

% AIC
for i = 1:length(b{3}(1,:))
    aic3(i) = dev{3}(i) + 2*length(b{3}(:,i));
    disp(['AIC for Model 3, Neuron ' num2str(i) ': ' num2str(aic3(i))]);
end

% KS-plot
figure(103);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))
    % Define history
    lambh = zeros(length(spikes_binned(:,i))-5,length(spikes_binned(1:5,i)));
    for k = 1:length(spikes_binned(1:5,i))
        lambh(:,k) = spikes_binned((6-k):(end-k),i);
    end
    
    % Lambda without history
    lambdaEst = exp(b{3}(1,i));
    % Add in history
    for k = 2:length(b{3}(:,i))
     lambdaEst = lambdaEst + b{3}(k,i)*lambh(:,k-1);
    end
    
    timestep = 1;
    lambdaInt = 0;
    j=0;

    for t=1:length(spikes_binned(:,i))-5
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes_binned(t,i))
            j = j + 1;
            KS(j,i) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted1 = sort( KS(:,i) );
    N = length( KSSorted1);
    figure(103);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted1, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 3: Neuron ' num2str(i) ': KS Plot']);

    ks_stat1(i) = max(abs(KSSorted1' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS Stats: ' num2str(ks_stat1)]);

AIC{3} = aic3;
ks_stat{3} = ks_stat1;
KSSorted{3} = KSSorted1;

%---------------- Model 4, Covariates: Short-term history-dependence ---------------
% Confidence intervals for each covariate
figure(4);
for i = 1:length(b{4}(1,:))
    subplot(2,5,i)
    errorbar(b{4}(:,i),stats{4}(i).se);
    hold on;
    plot(1:length(b{4}(:,i)),zeros(length(b{4}(:,i)),1),'r-');
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 4, Neuron ' num2str(i)]);
end

% P-values
for i = 1:length(b{4}(1,:))
    disp(['P-values for model 4, Neuron ' num2str(i) ': ' num2str(stats{4}(i).p')]);
end

% AIC
for i = 1:length(b{4}(1,:))
    aic4(i) = dev{4}(i) + 2*length(b{4}(:,i));
    disp(['AIC for Model 4, Neuron ' num2str(i) ': ' num2str(aic4(i))]);
end

% KS-plot
figure(104);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))
    % Define history
    lambh = zeros(length(spikes_binned(:,i))-125,length(spikes_binned(75:125,i)));
    for k = 1:length(spikes_binned(75:125,i))
        lambh(:,k) = spikes_binned((126-k):(end-k),i);
    end
    
    % Lambda without history
    lambdaEst = exp(b{4}(1,i));
    % Add in history
    for k = 2:length(b{4}(:,i))
     lambdaEst = lambdaEst + b{4}(k,i)*lambh(:,k-1);
    end
    
    timestep = 1;
    lambdaInt = 0;
    j=0;

    for t=1:length(spikes_binned(:,i))-125
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes_binned(t,i))
            j = j + 1;
            KS(j,i) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted1 = sort( KS(:,i) );
    N = length( KSSorted1);
    figure(104);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted1, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 4: Neuron ' num2str(i) ': KS Plot']);

    ks_stat1(i) = max(abs(KSSorted1' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS stats: ' num2str(ks_stat1)]);

AIC{4} = aic4;
ks_stat{4} = ks_stat1;
KSSorted{4} = KSSorted1;


%---------------- Model 5, Covariates: Long-term history-dependence ---------------
% Confidence intervals for each covariate
figure(5);
for i = 1:length(b{5}(1,:))
    subplot(2,5,i)
    errorbar(b{5}(:,i),stats{5}(i).se);
    hold on;
    plot(1:length(b{5}(:,i)),zeros(length(b{5}(:,i)),1),'r-');
    xlabel('\beta_i');
    ylabel('\beta value');
    title(['Model 5, Neuron ' num2str(i)]);
end

% P-values
for i = 1:length(b{5}(1,:))
    disp(['P-values for model 5, Neuron ' num2str(i) ': ' num2str(stats{5}(i).p')]);
end

% AIC
for i = 1:length(b{5}(1,:))
    aic5(i) = dev{5}(i) + 2*length(b{5}(:,i));
    disp(['AIC for Model 5, Neuron ' num2str(i) ': ' num2str(aic5(i))]);
end

% KS-plot
figure(105);
clear b_lin dev_lin stats_lin lamdaEst KS
for i = 1:length(spikes_binned(1,:))
    % Define history
    lambh = zeros(length(spikes_binned(:,i))-275,length(spikes_binned(200:275,i)));
    for k = 1:length(spikes_binned(200:275,i))
        lambh(:,k) = spikes_binned((276-k):(end-k),i);
    end
    
    % Lambda without history
    lambdaEst = exp(b{5}(1,i));
    % Add in history
    for k = 2:length(b{5}(:,i))
     lambdaEst = lambdaEst + b{5}(k,i)*lambh(:,k-1);
    end
    
    timestep = 1;
    lambdaInt = 0;
    j=0;

    for t=1:length(spikes_binned(:,i))-275
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes_binned(t,i))
            j = j + 1;
            KS(j,i) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted1 = sort( KS(:,i) );
    N = length( KSSorted1);
    figure(105);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted1, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Model 5: Neuron ' num2str(i) ': KS Plot']);

    ks_stat1(i) = max(abs(KSSorted1' - ((1:N)-.5)/N));
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(['KS Stats: ' num2str(ks_stat1)]);

AIC{5} = aic5;
ks_stat{5} = ks_stat1;
KSSorted{5} = KSSorted1;

end



