% create a network object, loading in structure data from a network file
filename = 'example.net';
NT = NetworkObj(filename);

%% plot network structure
plotopts = struct('labelnodes',1);
NT.plotNetwork(plotopts)

%% calculate analytic MFPTs to a set of targets
% which nodes are the targets?
targets = [3,6];

% MFPTs from each node to the first target hit
D = 1; % diffusivity
MFPTs = networkMFPTanalytic(NT,targets,D);

plotopts = struct();
plotopts.colornodes = MFPTs;
NT.plotNetwork(plotopts)
colorbar 