function MFPTs = networkMFPTanalytic(NT,targets,D)
% get the mfpt for diffusion on a network
% using analytic matrix approaches
 
% NT is a network object (created with NetworkObj)
% targets = list of target nodes (gives MFPT to any of the targets)
% D = diffusivity (default is 1)
% returns first passage from each initial node

%%
istarget = zeros(1,NT.nnode);
istarget(targets) = 1;
nottargets = find(~istarget);

Pmat = zeros(NT.nnode,NT.nnode);
Qvec = zeros(NT.nnode,1);
for nc = nottargets
    
    deg = NT.degrees(nc);
    
    % all edge lengths from this node
    lens = NT.edgelens(NT.nodeedges(nc,1:deg));
    lensinv = 1./lens;
    suminv = sum(lensinv);
    
    % get Pij at s=0 for this node
    Pmat(nc,NT.nodenodes(nc,1:deg)) = lensinv/ suminv;
    
    % get Qi at s=0 for this node
    Qvec(nc) = sum(lens)/suminv;
end
Pmat = Pmat(nottargets,nottargets);
Qvec = Qvec(nottargets);

% calculate MFPT to target nodes from each initial position
MFPTs = zeros(NT.nnode,1);
MFPTs(nottargets) = inv(eye(length(nottargets))-Pmat)*Qvec/2;

if (exist('D','var'))
    % scale by diffusivity if desired
    MFPTs = MFPTs/D;
end
end