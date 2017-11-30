function [b, dev, stats, lambda] = glm(file)
% Fits GLM for 3 specified models to data in file

% Error condition
if ~isa(file,'char')
    disp('Input argument must be in form ''[filename].mat''');
    return
end

load(file);

b = {};
dev = {};
stats = {};
lambda = {};

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

b{1} = b1;
dev{1} = dev1;
stats{1} = stats1;
lambda{1} = lambda1;


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

b{2} = b2;
dev{2} = dev2;
stats{2} = stats2;
lambda{2} = lambda2;

%---------------- Model 3, Covariates: Refractory period history-dependence ---------------
% Using history from 1 ms to 5 ms

for i = 1:length(spikes_binned(1,:))
    % Defining history
    lambh = zeros(length(spikes_binned(:,i))-5,length(spikes_binned(1:5,i)));
    for k = 1:length(spikes_binned(1:5,i))
        lambh(:,k) = spikes_binned((6-k):(end-k),i);
    end
    
    % GLM coefficients
    [b3(:,i),dev3(:,i),stats3(i)] = glmfit([lambh],spikes_binned(6:end,i),'poisson');
    % Computing lambda
    % Lambda without history
    lamb = exp(b3(1,i));
    % Adding history components
    for k = 2:length(b3(:,i))
     lamb = lamb + b3(k,i)*lambh(:,k-1);
    end
    
    lambda3(:,:,i) = lamb;
    
    disp(['Completed Neuron ' num2str(i) '.']);
end

b{3} = b3;
dev{3} = dev3;
stats{3} = stats3;
lambda{3} = lambda3;

% POSITION AND HISTORY DEPENDENCE
% % Positions to plot
% figure;
% [x_new,y_new]=meshgrid(-1:.1:1);
% y_new = flipud(y_new);
% x_new = fliplr(x_new);
% 
% for i = 1:length(spikes_binned(1,:))
%     % Defining history components
%     lambh = zeros(length(spikes_binned(:,i))-5,length(spikes_binned(1:5,i)));
%     for k = 1:length(spikes_binned(1:5,i))
%         lambh(:,k) = spikes_binned((6-k):(end-k),i);
%     end
%     
%     % GLM coefficients
%     [b3(:,i),dev3(:,i),stats3(i)] = glmfit([xN(6:end) yN(6:end) xN(6:end).^2 yN(6:end).^2 lambh],spikes_binned(6:end,i),'poisson');
%     disp(b3);
%     % Computing lambda
%     % Lambda without history
%     lamb = exp(b3(1,i) + b3(2,i)*x_new + b3(3,i)*y_new + b3(4,i).*x_new.^2 + b3(5,i).*y_new.^2);
%     % Adding history components
%     for k = 6:length(b3(:,i))
%      lamb = lamb + b3(k,i)*lambh(:,k-5);
%     end
% 
%     lamb(find(x_new.^2+y_new.^2>1))=nan;
%     lambda3(:,:,i) = lamb;
% 
%     % Plotting lambda and circle defining position limits
%     subplot(2,5,i)
%     h_mesh = mesh(x_new,y_new,lambda3(:,:,i),'AlphaData',0);
%     get(h_mesh,'AlphaData');
%     set(h_mesh,'AlphaData',0);
%     xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
%     hold on;
%     plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
% 
%     
%     disp(['Completed Neuron ' num2str(i) '.']);
% end

%---------------- Model 4, Covariates: Short-term history-dependence ---------------
% Using history from 75 ms to 125 ms

for i = 1:length(spikes_binned(1,:))
    % Defining history
    lambh = zeros(length(spikes_binned(:,i))-125,length(spikes_binned(75:125,i)));
    for k = 1:length(spikes_binned(75:125,i))
        lambh(:,k) = spikes_binned((126-k):(end-k),i);
    end
    
    % GLM coefficients
    [b4(:,i),dev4(:,i),stats4(i)] = glmfit([lambh],spikes_binned(126:end,i),'poisson');
    % Computing lambda
    % Lambda without history
    lamb = exp(b4(1,i));
    % Adding history components
    for k = 2:length(b4(:,i))
     lamb = lamb + b4(k,i)*lambh(:,k-1);
    end
    
    lambda4(:,:,i) = lamb;
    
    disp(['Completed Neuron ' num2str(i) '.']);
end

b{4} = b4;
dev{4} = dev4;
stats{4} = stats4;
lambda{4} = lambda4;

% POSITION AND HISTORY DEPENDENCE
% % Using history from 75 ms to 175 ms
% % Positions to plot
% figure;
% [x_new,y_new]=meshgrid(-1:.1:1);
% y_new = flipud(y_new);
% x_new = fliplr(x_new);
% 
% 
% for i = 1:length(spikes_binned(1,:))
%     % Defining history components
%     lambh = zeros(length(spikes_binned(:,i))-175,length(spikes_binned(75:175,i)));
%     for k = 1:length(spikes_binned(75:175,i))
%         lambh(:,k) = spikes_binned((176-k):(end-k),i);
%     end
%     
%     % GLM coefficients
%     [b4(:,i),dev4(:,i),stats4(i)] = glmfit([xN(176:end) yN(176:end) xN(176:end).^2 yN(176:end).^2 lambh],spikes_binned(176:end,i),'poisson');
% 
%     % Computing lambda
%     % Lambda without history
%     lamb = exp(b4(1,i) + b4(2,i)*x_new + b4(3,i)*y_new + b4(4,i).*x_new.^2 + b4(5,i).*y_new.^2);
%     % Adding history components
%     for k = 6:length(b4(:,i))
%      lamb = lamb + b4(k,i).*lambh(:,k-5);
%     end
%     
%     lamb(find(x_new.^2+y_new.^2>1))=nan;
%     lambda4(:,:,i) = lamb;
% 
%     % Plotting lambda and circle defining position limits
%     subplot(2,5,i)
%     h_mesh = mesh(x_new,y_new,lambda4(:,:,i),'AlphaData',0);
%     get(h_mesh,'AlphaData');
%     set(h_mesh,'AlphaData',0);
%     xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
%     hold on;
%     plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
% 
%     
%     disp(['Completed Neuron ' num2str(i) '.']);
% end


%---------------- Model 5, Covariates: Long-term history-dependence ---------------
% Using history from 200 ms to 275 ms

for i = 1:length(spikes_binned(1,:))
    % Defining history
    lambh = zeros(length(spikes_binned(:,i))-275,length(spikes_binned(200:275,i)));
    for k = 1:length(spikes_binned(200:275,i))
        lambh(:,k) = spikes_binned((276-k):(end-k),i);
    end
    
    % GLM coefficients
    [b5(:,i),dev5(:,i),stats5(i)] = glmfit([lambh],spikes_binned(276:end,i),'poisson');
    % Computing lambda
    % Lambda without history
    lamb = exp(b5(1,i));
    % Adding history components
    for k = 2:length(b5(:,i))
     lamb = lamb + b5(k,i)*lambh(:,k-1);
    end
    
    lambda5(:,:,i) = lamb;
    
    disp(['Completed Neuron ' num2str(i) '.']);
end

b{5} = b5;
dev{5} = dev5;
stats{5} = stats5;
lambda{5} = lambda5;

% % Using history from 200 ms to 275 ms
% % Positions to plot
% figure;
% [x_new,y_new]=meshgrid(-1:.1:1);
% y_new = flipud(y_new);
% x_new = fliplr(x_new);
% 
% for i = 1:length(spikes_binned(1,:))
%     % Defining history components
%     lambh = zeros(length(spikes_binned(:,i))-275,length(spikes_binned(200:275,i)));
%     for k = 1:length(spikes_binned(200:275,i))
%         lambh(:,k) = spikes_binned((276-k):(end-k),i);
%     end
%     
%     % GLM coefficients
%     [b5(:,i),dev5(:,i),stats5(i)] = glmfit([xN(201:end) yN(201:end) xN(201:end).^2 yN(201:end).^2 lambh],spikes_binned(201:end,i),'poisson');
% 
%     % Computing lambda
%     % Lambda without history
%     lamb = exp(b5(1,i) + b5(2,i)*x_new + b5(3,i)*y_new + b5(4,i).*x_new.^2 + b5(5,i).*y_new.^2);
%     % Adding history components
%     for k = 6:length(b5(:,i))
%      lamb = lamb + b5(k,i).*lambh(:,k-5);
%     end
% 
%     lamb(find(x_new.^2+y_new.^2>1))=nan;
%     lambda5(:,:,i) = lamb;
% 
%     % Plotting lambda and circle defining position limits
%     subplot(2,5,i)
%     h_mesh = mesh(x_new,y_new,lambda5(:,:,i),'AlphaData',0);
%     get(h_mesh,'AlphaData');
%     set(h_mesh,'AlphaData',0);
%     xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
%     hold on;
%     plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
% 
%     
%     disp(['Completed Neuron ' num2str(i) '.']);
% end






%---------------- Model ?, Covariates: Velocity? -------------------- 





end
