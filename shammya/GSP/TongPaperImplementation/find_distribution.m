data_files = ["Armstrong2019","Barrett2019","Santan2019","Cedar2019","GWC2019","ISTB32019","MercadoA2019","Peralta2019","SP2019","Wrigley2019"];

all_distributions = containers.Map();
all_data = containers.Map();
for i = 1:length(data_files)
    data = csvread(strcat(data_files(i),".csv"));
    data = data(1:8640);
    loc= find(data<=0);
    data(loc) = [];
    % [~,score] = pca(data,"NumComponents",1);
    % GMModels = cell(3,1); % Preallocation
    % options = statset("MaxIter",1000);
    % rng(1); % For reproducibility
    % 
    % for j = 1:4
    %     GMModels{j} = fitgmdist(score,j,"Options",options);
    %     fprintf("\n GM Mean for %i Component(s)\n",j)
    %     Mu = GMModels{j}.mu
    % end
    % pd_kernel = fitdist(data,"Kernel");
    % pd_normal = fitdist(data,"Normal");
    % x = min(data):0.1:max(data);
    % 
    % pdf_kernel = pdf(pd_kernel,x);
    % pdf_normal = pdf(pd_normal,x);
    BIC = zeros(1,10);
    GMModels = cell(1,10);
    options = statset("MaxIter",5000);
    for k = 1:10
        GMModels{k} = fitgmdist(data,k,"Options",options);
        BIC(k)= GMModels{k}.BIC;
    end

    [minBIC,numComponents] = min(BIC);
    BestModel = GMModels{numComponents};
    all_distributions(data_files(i)) = BestModel;
    all_data (data_files(i)) = data;
end
save('distribution.mat','all_distributions');
save('all_data.mat','all_data');
save('data_file_list.mat','data_files');


% plot(x,pdf_kernel,"Color","b","LineWidth",2);
% hold on;
% plot(x,pdf_normal,"Color","r","LineStyle",":","LineWidth",2);
% hold on 
% plot(x,pdf(BestModel,x"),"LineWidth",2)
% hold on 
% histogram(data,"Normalization","pdf");
% legend("Kernel Distribution","Normal Distribution","GM","Location","NorthEast");
% hold off;
