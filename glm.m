function [b1, dev1, stats1, lambda1, b2, dev2, stats2, lambda2, b3, dev3, stats3, lambda3, b4, dev4, stats4, lambda4, b5, dev5, stats5, lambda5] = glm(file)
% Fits GLM for 3 specified models to data in file

% Error condition
if ~isa(file,'char')
    disp('Input argument must be in form ''[filename].mat''');
    return
end

load(file);

%---------------- Model 1, Covariates: X, Y, X^2, Y^2 --------------------

% Positions to plot
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b1(:,i),dev1(:,i),stats1(i)] = glmfit([xN yN xN.^2 yN.^2],spikes_binned(:,i),'poisson');

    % Computing lambda
    lamb = exp(b1(1,i) + b1(2,i)*x_new + b1(3,i)*y_new + b1(4,i).*x_new.^2 + b1(5,i).*y_new.^2);
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda1(:,:,i) = lamb;
    % NEED TO SAVE LAMBDAS

    % Plotting lambda and circle defining position limits
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda1(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    
    disp(['Completed Neuron ' num2str(i) '.']);
end




%---------------- Model 2, Covariates: Sinusoids ---------------

% Positions to plot
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b2(:,i),dev2(:,i),stats2(i)] = glmfit([sin(2.5*pi*xN) sin(2.5*pi*yN)],spikes_binned(:,i),'poisson');

    % Computing lambda
    lamb = exp(b2(1,i) + b2(2,i)*sin(2.5*pi*x_new) + b2(3,i)*sin(2.5*pi*y_new));
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda2(:,:,i) = lamb;

    % Plotting lambda and circle defining position limits
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda2(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    
    disp(['Completed Neuron ' num2str(i) '.']);
end



%---------------- Model 3, Covariates: Refractory period history-dependence ---------------

% Positions to plot
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b3(:,i),dev3(:,i),stats3(i)] = glmfit([xN(6:end) yN(6:end) xN(6:end).^2 yN(6:end).^2 ones(size(xN(6:end))).*spikes_binned(1:5,i)'],spikes_binned(6:end,i),'poisson');

    % Computing lambda
    % History-dependent portion
    lambh = 0;
    for k = 1:length(spikes_binned(1:5,i))
    lambh = lambh + b3(k+5,i)*spikes_binned(k,i)';
    end
    lamb = exp(b3(1,i) + b3(2,i)*x_new + b3(3,i)*y_new + b3(4,i).*x_new.^2 + b3(5,i).*y_new.^2 + lambh);
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda3(:,:,i) = lamb;

    % Plotting lambda and circle defining position limits
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda3(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    
    disp(['Completed Neuron ' num2str(i) '.']);
end

%---------------- Model 4, Covariates: Short-term history-dependence ---------------
% Using history from 75 ms to 175 ms
% Positions to plot
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b4(:,i),dev4(:,i),stats4(i)] = glmfit([xN(176:end) yN(176:end) xN(176:end).^2 yN(176:end).^2 ones(size(xN(176:end))).*spikes_binned(75:175,i)'],spikes_binned(176:end,i),'poisson');

    % Computing lambda
    % History-dependent portion
    lambh = 0;
    for k = 1:length(spikes_binned(75:175,i))
    lambh = lambh + b4(k+5,i)*spikes_binned(k,i)';
    end
    lamb = exp(b4(1,i) + b4(2,i)*x_new + b4(3,i)*y_new + b4(4,i).*x_new.^2 + b4(5,i).*y_new.^2 + lambh);
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda4(:,:,i) = lamb;

    % Plotting lambda and circle defining position limits
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda4(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    
    disp(['Completed Neuron ' num2str(i) '.']);
end

%---------------- Model 5, Covariates: Long-term history-dependence ---------------
% Using history from 200 ms to 275 ms
% Positions to plot
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b5(:,i),dev5(:,i),stats5(i)] = glmfit([xN(201:end) yN(201:end) xN(201:end).^2 yN(201:end).^2 ones(size(xN(201:end))).*spikes_binned(200:275,i)'],spikes_binned(201:end,i),'poisson');

    % Computing lambda
    % History-dependent portion
    lambh = 0;
    for k = 1:length(spikes_binned(200:275,i))
    lambh = lambh + b5(k+5,i)*spikes_binned(k,i)';
    end
    lamb = exp(b5(1,i) + b5(2,i)*x_new + b5(3,i)*y_new + b5(4,i).*x_new.^2 + b5(5,i).*y_new.^2 + lambh);
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda5(:,:,i) = lamb;

    % Plotting lambda and circle defining position limits
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda5(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    
    disp(['Completed Neuron ' num2str(i) '.']);
end






%---------------- Model ?, Covariates: Velocity? -------------------- 





end
