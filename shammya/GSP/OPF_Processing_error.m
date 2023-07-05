clc;
clear all;
close all; 
define_constants;
opf_results_gsp = load('OPF_GSP_600.mat');
matpower_opf = opf_results_gsp.opf_case_results;
total_demand = sum(matpower_opf.bus(:,PD))
total_power_matpower = sum((matpower_opf.gen(:,PG)))
pg = opf_results_gsp.pg;
total_power_gsp = matpower_opf.baseMVA*sum(pg)
Cs = opf_results_gsp.Cs;
matpower_cost= sum(Cs.*(matpower_opf.gen(:,PG)).*(matpower_opf.gen(:,PG)))

generation_bus = opf_results_gsp.generation_bus;
cost = matpower_opf.baseMVA^2*sum(Cs.*pg(generation_bus).*pg(generation_bus))
telapsed_gsp = opf_results_gsp.telapsed_gsp


opf_results_ngsp = load('OPF_NOGSP_600.mat');
pg_ngsp = opf_results_ngsp.pg_convex;
total_power_ngsp = matpower_opf.baseMVA*sum(pg_ngsp)
Cs_ngsp = opf_results_ngsp.Cs;
generation_bus_ngsp = opf_results_ngsp.generation_bus;
cost_ngsp = matpower_opf.baseMVA^2* sum(Cs_ngsp.*pg_ngsp(generation_bus_ngsp).*pg_ngsp(generation_bus_ngsp))
telapsed_ngsp = opf_results_ngsp.telapsed