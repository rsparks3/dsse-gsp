[~,Yf,~] = makeYbus(loadcase('case85.m'));

ammeter_location = load('ammeter.mat').ammeter_location;
power_flow = load('power_flow.mat').power_flow;

real_injection = zeros(length(power_flow),84+length(ammeter_location));
reactive_injection = zeros(length(power_flow),84+length(ammeter_location));

input_data = zeros(length(power_flow),84*2+length(ammeter_location));
output_data = zeros(length(power_flow),85*2);
define_constants;
% output_data = [];
for i = 1:length(power_flow)
    mpc = power_flow{i};
    pf = mpc.branch(:,PF)+random_distribution(84);
    qf = mpc.branch(:,QF)+random_distribution(84);
    current = Yf*mpc.bus(:,VM);
    current = current(ammeter_location) ;
    current = abs(current) + random_distribution(length(current),0.001);
    input_data(i,:) = [pf' qf' current'];
    vmagnitude = mpc.bus(:,VM);
    vangle = mpc.bus(:,VA);
    output_data(i,:) = [vmagnitude' vangle'];
%     output_data = [output_data;mpc.bus(2:end,VM);mpc.bus(2:end,VA)];
   
end

save('input_data.mat','input_data');
save('output_data.mat','output_data');
 


































