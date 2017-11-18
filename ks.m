function ks(file)
% Creates KS plot for data in file

% Error condition
if ~isa(file,'char')
    disp('Input argument must be in form ''[filename].mat''');
    return
end

load(file);




%---------------- Model 1, Covariates: X, Y, X^2, Y^2 --------------------

figure;
for i = 1:length(spikes_binned(1,:))
    [b_lin(:,i),dev_lin(:,i),stats_lin(:,i)] = glmfit([xN yN xN.^2 yN.^2],spikes_binned(:,i),'poisson');

    lambdaEst = exp(b_lin(1,i) + b_lin(2,i)*xN + b_lin(3,i)*yN + b_lin(4,i)*xN.^2 + b_lin(5,i)*yN.^2); 

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
%     disp(['sorted ' num2str(i) ' is ' num2str(sort(KS(:,i)'))]);
    KSSorted = sort( KS(:,i) );
    N = length( KSSorted);
    subplot(2,5,i);
    plot( ([1:N]-.5)/N, KSSorted, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
    axis( [0 1 0 1] );
    xlabel('Uniform CDF');
    ylabel('Empirical CDF of Rescaled ISIs');
    title(['Neuron ' num2str(i) ': KS Plot, 95% Confidence']);

    ks_stat(i) = max(abs(KSSorted' - ((1:N)-.5)/N));
    % NEED TO SAVE KS STATISTIC, MAYBE KS MATRIX TOO?
    
    disp(['Completed Neuron ' num2str(i) '.']);

end


disp(ks_stat);



%---------------- Model 2, Covariates: ? -------------------- 



%---------------- Model 3, Covariates: ? --------------------













end