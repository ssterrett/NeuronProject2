function glm(file)
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
    [b(:,i),dev(:,i),stats(:,1)] = glmfit([xN yN xN.^2 yN.^2],spikes_binned(:,i),'poisson');
    % NEED TO SAVE THESE VARIABLES

    % Computing lambda
    lambda = exp(b(1,i) + b(2,i)*x_new + b(3,i)*y_new + b(4,i).*x_new.^2 + b(5,i).*y_new.^2);
    lambda(find(x_new.^2+y_new.^2>1))=nan;
    % NEED TO SAVE LAMBDAS

    % Plotting lambda and circle defining position limits
    subplot(2,5,i)
    h_mesh = mesh(x_new,y_new,lambda,'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    xlabel('x'); ylabel('y'); title(['Lambda, Neuron ' num2str(i)]);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

    
    disp(['Completed Neuron ' num2str(i) '.']);
end






%---------------- Model 2, Covariates: ? -------------------- 



%---------------- Model 3, Covariates: ? --------------------







end
