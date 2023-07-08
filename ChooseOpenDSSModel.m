function [file_directory,slack_bus_name] = ChooseOpenDSSModel(BusSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if BusSize == 10
    file_directory = 'Compile C:\models\test10bus\10bus.dss';
    slack_bus_name = upper('n1');
elseif BusSize == 123
    file_directory = 'Compile C:\models\123Busm\IEEE123Master.DSS';
    slack_bus_name = upper('150');
elseif BusSize == 37
    file_directory = 'Compile C:\models\37Bus\ieee37.DSS';
    slack_bus_name = upper('sourcebus');
elseif BusSize == 13
    file_directory = 'Compile C:\models\13Bus\IEEE13Nodeckt.dss';
    slack_bus_name = upper('SourceBus');
elseif BusSize == 8500
    file_directory = 'Compile C:\models\8500-Node\Master-unbal.dss';
elseif BusSize == 342
    file_directory = 'Compile C:\models\LVTestCaseNorthAmerican\Master.DSS';
    slack_bus_name = upper('P1');
elseif BusSize == 34 
    file_directory = 'Compile C:\models\34Bus\ieee34Mod2.dss';
    slack_bus_name = upper('sourcebus');
elseif BusSize == 1102 
    file_directory = 'Compile C:\models\circuit2\circuit2.dss';
    slack_bus_name = upper('sourcebus');
elseif BusSize == 1105 
    file_directory = 'Compile C:\models\circuit5\circuit5.dss';
    slack_bus_name = upper('sourcebus');
elseif BusSize == 1103 
    file_directory = 'Compile C:\models\circuit3\circuit3.dss';
    slack_bus_name = upper('sourcebus');
else
    file_directory = '';
    slack_bus_name = '';
end
end

