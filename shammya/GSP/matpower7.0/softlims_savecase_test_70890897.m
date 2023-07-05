function mpc = softlims_savecase_test_70890897
%SOFTLIMS_SAVECASE_TEST_70890897

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin	lam_P	lam_Q	mu_Vmax	mu_Vmin
mpc.bus = [
	10	3	0	0	0	0	1	0.9	0	345	1	1.1	0.9	45.0000	-10.0000	0.0000	2075.4687;
	20	2	0	0	0	0	1	1.1	5.53422214	345	1	1.1	0.9	40.0000	20.0000	2207.2707	0.0000;
	30	2	0	0	0	0	1	1.1	17.0346969	345	1	1.1	0.9	28.0000	0.0000	6822.1173	0.0000;
	40	1	0	0	0	0	1	0.934104949	-0.795292694	345	1	0.9	0.9	44.7148	-10.5476	75000.0000	0.0000;
	50	1	90	30	0	0	1	0.960199579	0.255941117	345	1	1.1	0.9	44.4574	-3.9386	0.0000	0.0000;
	60	1	0	0	0	0	1	1.07435325	10.2605553	345	1	1.1	0.9	29.0032	3.6281	0.0000	0.0000;
	70	1	100	35	0	0	1	1.05630719	4.21056184	345	1	1.1	0.9	40.9359	15.2209	0.0000	0.0000;
	80	1	0	0	0	0	1	1.0660722	3.46199522	345	1	1.1	0.9	41.5392	22.5392	0.0000	0.0000;
	90	1	125	50	0	0	1	0.947380968	-3.24216966	345	1	1.1	1.05	55.1813	48.9425	0.0000	75000.0000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf	mu_Pmax	mu_Pmin	mu_Qmax	mu_Qmin
mpc.gen = [
	30	237.881973	62.2212912	300	-300	1.1	100	1	225	10	0	0	0	0	0	0	0	0	0	0	0	3.0000	0.0000	0.0000	0.0000;
	20	0	0	300	-300	1.025	100	0	300	10	0	0	0	0	0	0	0	0	0	0	0	0.0000	0.0000	0.0000	0.0000;
	10	20.2584453	-53.148382	300	-45	0.9	100	1	250	25	0	0	0	0	0	0	0	0	0	0	0	0.0000	5.0000	0.0000	10.0000;
	20	67.8452237	60.9399468	50	-300	1.1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0	0.0000	0.0000	20.0000	0.0000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax	Pf	Qf	Pt	Qt	mu_Sf	mu_St	mu_angmin	mu_angmax
mpc.branch = [
	10	40	0	0.0576	0	250	250	250	0	0	1	-60	60	20.2584	-53.1484	-20.2584	55.4489	0.0000	0.0000	0.0000	0.0000;
	40	50	0.017	0.092	0.158	0	250	250	0	0	1	-60	60	-22.0006	-29.1584	22.1915	16.0147	0.0000	0.0000	0.0000	0.0000;
	50	60	0.039	0.17	0.358	120	150	150	0	0	1	-5	60	-112.1915	-46.0147	117.8842	33.6646	0.0000	10.0000	20.0000	0.0000;
	30	60	0	0.0586	0	0	300	300	0	0	1	-60	5	237.8820	62.2213	-237.8820	-32.9410	0.0000	0.0000	0.0000	20.0000;
	30	60	0	0.0586	0	0	300	300	0	0	0	-60	60	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000;
	60	70	0.0119	0.1008	0.209	120	150	150	0	0	1	-60	60	119.9978	-0.7236	-118.5000	-10.3107	8.8116	0.0000	0.0000	0.0000;
	70	80	0.0085	0.072	0.149	0	250	250	0	0	1	-60	60	18.5000	-24.6893	-18.4535	8.3036	0.0000	0.0000	0.0000	0.0000;
	80	20	0	0.0625	0	250	250	250	0	0	1	-60	60	-67.8452	-56.6442	67.8452	60.9399	0.0000	0.0000	0.0000	0.0000;
	80	90	0.032	0.161	0.306	250	250	250	0	0	1	-60	60	86.2987	48.3405	-82.9853	-62.7909	0.0000	0.0000	0.0000	0.0000;
	90	40	0.01	0.085	0.176	250	250	250	0	0	1	-60	60	-42.0147	12.7909	42.2590	-26.2905	0.0000	0.0000	0.0000	0.0000;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	3000	0	2	25	0;
	2	2000	0	2	40	0;
	2	1500	0	2	50	0;
	2	2000	0	2	40	0;
];

%%-----  OPF Soft Limit Data  -----%%
%% VMAX soft limit data
mpc.softlims.VMAX.hl_mod = 'scale';     %% type of hard limit modification
mpc.softlims.VMAX.hl_val = [            %% value(s) used to set new hard limit
	1.5;
];
mpc.softlims.VMAX.idx = [               %% bus numbers for soft voltage limit
	1;
	2;
	3;
	4;
	5;
	6;
	7;
	8;
	9;
];
mpc.softlims.VMAX.cost = [              %% violation cost coefficient
	75000;
	75000;
	75000;
	75000;
	75000;
	75000;
	75000;
	75000;
	75000;
];
mpc.softlims.VMAX.overload = [          %% violation of original limit
	0;
	0;
	0;
	0.03410494894;
	0;
	0;
	0;
	0;
	0;
];
mpc.softlims.VMAX.ovl_cost = [          %% overload cost
	0;
	0;
	0;
	2557.871171;
	0;
	0;
	0;
	0;
	0;
];

%% VMIN soft limit data
mpc.softlims.VMIN.hl_mod = 'scale';     %% type of hard limit modification
mpc.softlims.VMIN.hl_val = [            %% value(s) used to set new hard limit
	0.5;
];
mpc.softlims.VMIN.idx = [               %% bus numbers for soft voltage limit
	1;
	2;
	8;
	9;
];
mpc.softlims.VMIN.busnum = [            %% bus numbers for soft voltage limit
	10;
	20;
	80;
	90;
];
mpc.softlims.VMIN.cost = [              %% violation cost coefficient
	75000;
	75000;
	75000;
	75000;
];
mpc.softlims.VMIN.overload = [          %% violation of original limit
	0;
	0;
	0;
	0;
	0;
	0;
	0;
	0;
	0.1026190319;
];
mpc.softlims.VMIN.ovl_cost = [          %% overload cost
	0;
	0;
	0;
	0;
	0;
	0;
	0;
	0;
	7696.42739;
];

%% PMAX soft limit data
mpc.softlims.PMAX.hl_mod = 'remove';    %% type of hard limit modification
mpc.softlims.PMAX.idx = [               %% gen matrix row indices
	1;
	3;
	4;
];
mpc.softlims.PMAX.cost = [              %% violation cost coefficient
	3;
	1000;
	1000;
];
mpc.softlims.PMAX.overload = [          %% violation of original limit
	12.88197308;
	0;
	0;
	0;
];
mpc.softlims.PMAX.ovl_cost = [          %% overload cost
	38.64591923;
	0;
	0;
	0;
];

%% PMIN soft limit data
mpc.softlims.PMIN.hl_mod = 'scale';     %% type of hard limit modification
mpc.softlims.PMIN.hl_val = [            %% value(s) used to set new hard limit
	0;
];
mpc.softlims.PMIN.idx = [               %% gen matrix row indices
	1;
	3;
	4;
];
mpc.softlims.PMIN.cost = [              %% violation cost coefficient
	5;
	5;
	5;
];
mpc.softlims.PMIN.overload = [          %% violation of original limit
	0;
	0;
	4.74155467;
	0;
];
mpc.softlims.PMIN.ovl_cost = [          %% overload cost
	0;
	0;
	23.70777335;
	0;
];

%% QMAX soft limit data
mpc.softlims.QMAX.hl_mod = 'remove';    %% type of hard limit modification
mpc.softlims.QMAX.idx = [               %% gen matrix row indices
	1;
	3;
	4;
];
mpc.softlims.QMAX.cost = [              %% violation cost coefficient
	1000;
	1000;
	20;
];
mpc.softlims.QMAX.overload = [          %% violation of original limit
	0;
	0;
	0;
	10.93994678;
];
mpc.softlims.QMAX.ovl_cost = [          %% overload cost
	0;
	0;
	0;
	218.7989356;
];

%% QMIN soft limit data
mpc.softlims.QMIN.hl_mod = 'remove';    %% type of hard limit modification
mpc.softlims.QMIN.idx = [               %% gen matrix row indices
	1;
	3;
	4;
];
mpc.softlims.QMIN.cost = [              %% violation cost coefficient
	1000;
	10;
	1000;
];
mpc.softlims.QMIN.overload = [          %% violation of original limit
	0;
	0;
	8.14838203;
	0;
];
mpc.softlims.QMIN.ovl_cost = [          %% overload cost
	0;
	0;
	81.4838203;
	0;
];

%% RATE_A soft limit data
mpc.softlims.RATE_A.hl_mod = 'scale';   %% type of hard limit modification
mpc.softlims.RATE_A.hl_val = [          %% value(s) used to set new hard limit
	1.5;
];
mpc.softlims.RATE_A.idx = [             %% branch matrix row indices
	1;
	3;
	6;
	8;
	9;
	10;
];
mpc.softlims.RATE_A.cost = [            %% violation cost coefficient
	10;
	10;
	10;
	10;
	10;
	10;
];
mpc.softlims.RATE_A.overload = [        %% violation of original limit
	0;
	0;
	2.596817619;
	0;
	0;
	0;
	0;
	0;
	0;
	0;
];
mpc.softlims.RATE_A.ovl_cost = [        %% overload cost
	0;
	0;
	25.96817619;
	0;
	0;
	0;
	0;
	0;
	0;
	0;
];

%% ANGMAX soft limit data
mpc.softlims.ANGMAX.hl_mod = 'replace'; %% type of hard limit modification
mpc.softlims.ANGMAX.hl_val = [          %% value(s) used to set new hard limit
	360;
];
mpc.softlims.ANGMAX.idx = [             %% branch matrix row indices
	1;
	2;
	3;
	4;
	6;
	7;
	8;
	9;
	10;
];
mpc.softlims.ANGMAX.cost = [            %% violation cost coefficient
	1000;
	1000;
	1000;
	20;
	1000;
	1000;
	1000;
	1000;
	1000;
];
mpc.softlims.ANGMAX.overload = [        %% violation of original limit
	0;
	0;
	0;
	1.774141588;
	0;
	0;
	0;
	0;
	0;
	0;
];
mpc.softlims.ANGMAX.ovl_cost = [        %% overload cost
	0;
	0;
	0;
	35.48283177;
	0;
	0;
	0;
	0;
	0;
	0;
];

%% ANGMIN soft limit data
mpc.softlims.ANGMIN.hl_mod = 'replace'; %% type of hard limit modification
mpc.softlims.ANGMIN.hl_val = [          %% value(s) used to set new hard limit
	-360;
];
mpc.softlims.ANGMIN.idx = [             %% branch matrix row indices
	1;
	2;
	3;
	4;
	6;
	7;
	8;
	9;
	10;
];
mpc.softlims.ANGMIN.cost = [            %% violation cost coefficient
	1000;
	1000;
	20;
	1000;
	1000;
	1000;
	1000;
	1000;
	1000;
];
mpc.softlims.ANGMIN.overload = [        %% violation of original limit
	0;
	0;
	5.004614232;
	0;
	0;
	0;
	0;
	0;
	0;
	0;
];
mpc.softlims.ANGMIN.ovl_cost = [        %% overload cost
	0;
	0;
	100.0922846;
	0;
	0;
	0;
	0;
	0;
	0;
	0;
];
