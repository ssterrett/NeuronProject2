function [b1, dev1, stats1, lambda1, b4, dev4, stats4, lambda4] = glm(file)
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






%---------------- Model 2, Covariates: Velocity? -------------------- 



%---------------- Model 3, Covariates: History dependence? ---------------





%---------------- Model 4, Covariates: Sinusoids ---------------

% Positions to plot
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

for i = 1:length(spikes_binned(1,:))
    
    % GLM coefficients
    [b4(:,i),dev4(:,i),stats4(i)] = glmfit([sin(2.5*pi*xN) sin(2.5*pi*yN)],spikes_binned(:,i),'poisson');

    % Computing lambda
    lamb = exp(b4(1,i) + b4(2,i)*sin(2.5*pi*x_new) + b4(3,i)*sin(2.5*pi*y_new));
    lamb(find(x_new.^2+y_new.^2>1))=nan;
    lambda4(:,:,i) = lamb;
    % NEED TO SAVE LAMBDAS

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



end
