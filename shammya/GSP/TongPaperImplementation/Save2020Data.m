new_data_files = ["Armstrong2020","Barrett2020","Santan2020","Cedar2020","GWC2020","ISTB32020","MercadoA2020","Peralta2020","SP2020","Wrigley2020"];

all_data_2020 = containers.Map();
for i = 1:length(new_data_files)
    data = csvread(strcat(new_data_files(i),".csv"));
    data = data(1:8640);
    loc= find(data<=0);
    data(loc) = [];
    all_data_2020 (new_data_files(i)) = data;
end


save('all_data_2020.mat','all_data_2020');
save('new_data_file_list.mat','new_data_files');

