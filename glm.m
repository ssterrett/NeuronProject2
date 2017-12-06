function [b, dev, stats, lambda, lambdaAll] = glm(file)
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
lambdaAll = {};

%---------------- Model 1a, Covariates: X, Y, X^2, Y^2 --------------------

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

%---------------- Model 1b, Covariates: X, Y, X^2, Y^2, wall dist, speed, refractory, short, long --------------------
% Refractory history 1 to 5 ms
% Shorter term history 80 to 110 ms
% Longer term history 240 to 260 ms

% Positions to plot
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

% Computing speed
vN = zeros(length(xN),1);
vN(1) = 0;
for j = 2:(length(xN)-1)
    vN(j) = sqrt((xN(j)-xN(j-1)).^2 +(yN(j)-yN(j-1)).^2);
end
vN(length(xN)) = vN(length(xN)-1);

figure(1000); figure(1001);
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
    
    
    % GLM coefficients
    [b2(:,i),dev2(:,i),stats2(i)] = glmfit([xN(261:end) yN(261:end) xN(261:end).^2 yN(261:end).^2 abs(1-sqrt(xN(261:end).^2 + yN(261:end).^2)) vN(261:end) lambhr lambhs lambhl],spikes_binned(261:end,i),'poisson');

    % Computing lambda
    % Without history
    lamb = exp(b2(1,i) + b2(2,i)*x_new + b2(3,i)*y_new + b2(4,i).*x_new.^2 + b2(5,i).*y_new.^2 + b2(6,i).*abs(1-sqrt(x_new.^2 + y_new.^2)));
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda2(:,:,i) = lamb;
    % All terms
    for k = 1:(length(b2(:,i)))
        lambdaAll_1(k,i) = exp(b2(k,i));
    end

    % Plotting positional lambda and circle defining position limits
    figure(1000);
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda2(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
    
    % Plotting all individual exp(beta) values
    figure(1001);
    subplot(2,5,i)
    plot(1:length(b2(:,i)),lambdaAll_1(:,i),'r-')
    xlabel('Covariate number'); ylabel('e^{\beta_i}'); 
    title(['Neuron ' num2str(i) ': e^{\beta_i}']);
    
    disp(['Completed Neuron ' num2str(i) '.']);
end

b{2} = b2;
dev{2} = dev2;
stats{2} = stats2;
lambda{2} = lambda2;
lambdaAll{2} = lambdaAll_1;

%---------------- Model 2, Covariates: Sinusoids ---------------

% Positions to plot
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);
figure(1002); figure(1003);

% Computing speed
vN = zeros(length(xN),1);
vN(1) = 0;
for j = 2:(length(xN)-1)
    vN(j) = sqrt((xN(j)-xN(j-1)).^2 +(yN(j)-yN(j-1)).^2);
end
vN(length(xN)) = vN(length(xN)-1);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b3(:,i),dev3(:,i),stats3(i)] = glmfit([sin(2.2*pi*(yN+0.8*xN)) sin(3*pi*(yN-0.35*xN)-9) vN],spikes_binned(:,i),'poisson');

    % Computing lambda
    lamb = exp(b3(1,i) + b3(2,i)*sin(2.2*pi*(y_new+0.8*x_new)) + b3(3,i)*sin(3*pi*(y_new-0.35*x_new)-9));
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda3(:,:,i) = lamb;
    
    % All terms
    for k = 1:(length(b3(:,i)))
        lambdaAll_2(k,i) = exp(b3(k,i));
    end

    % Plotting lambda and circle defining position limits
    figure(1002);
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda3(:,:,i),'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    % Plotting all individual exp(beta) values
    figure(1003);
    subplot(2,5,i)
    plot(1:length(b3(:,i)),lambdaAll_2(:,i),'r-')
    xlabel('Covariate number'); ylabel('e^{\beta_i}'); 
    title(['Neuron ' num2str(i) ': e^{\beta_i}']);
    
    
    disp(['Completed Neuron ' num2str(i) '.']);
end

b{3} = b3;
dev{3} = dev3;
stats{3} = stats3;
lambda{3} = lambda3;
lambdaAll{3} = lambdaAll_2;







%---------------- Model 3, Covariates: Head direction, Speed -------------------- 





end
