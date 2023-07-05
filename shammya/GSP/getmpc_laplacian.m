function out = getmpc_laplacian(fname)
%%% returns a structure with fields L for the laplacian and bus_name wich
%%% is a cell array of bus names as well as busnumber which is the Matpower
%%% EXTERNAL bus number
%%% the buses are reordered so that the Laplacian has structure
%%%          Lgen-gen     Lgen-nongen
%%%          Lnongen-gen  Lnongen-nongen

define_constants; % this creates a bunch of variables to index into mpc matrices.
mpc = loadcase(fname);
% for some of the commands to work you need to have all bus numbers
% consecutive starting at 1, etc. ext2int does that for you.
mpcint = ext2int(mpc);
%% catch some possible issues that are currently not handled
if size(mpc.bus,1) ~= size(mpcint.bus,1)
    error('getmpc_laplacian: some buses were removed in conversion to internal indexing. This is currently not handled.')
end
%% some useful numbers to have around
nb = size(mpc.bus,1); %this is the number of buses
%% This is the laplacian
B = makeBdc(mpcint);
L = speye(nb) - logical(B);
L = L + sparse(1:nb,1:nb, -sum(L,2),nb,nb);
%% let's make a boolean mask for buses
genbus = mpcint.gen(:,GEN_BUS); % this is the bus number for each generator
%check that the generator is on: (shouldn't be an issue because of internal
%conversion...bus still)
gen_status = mpcint.gen(:,GEN_STATUS) > 0;
% now you have only active generators:
genbus = genbus(gen_status);

% now make sure there are no repeat indices:
genbus = unique(genbus);

genmask = false(nb,1);
genmask(genbus) = true;
%% Order laplacian matrix and bus vectors
% Bord = [B(bool_mask, bool_mask), B(bool_mask, ~bool_mask);
%         B(~bool_mask, bool_mask), B(~bool_mask, ~bool_mask)];
out  = struct();
out.L        = [L( genmask, genmask), L( genmask, ~genmask);
                L(~genmask, genmask), L(~genmask, ~genmask)];
out.B = B;
out.busnum   = [mpc.bus(genmask,BUS_I); mpc.bus(~genmask,BUS_I)];
out.genmask  = genmask;
if isfield(mpc, 'bus_name')
  out.bus_name = vertcat(mpc.bus_name(genmask), mpc.bus_name(~genmask));
end