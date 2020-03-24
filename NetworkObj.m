% class definition for a network object
% containing data on network structure and basic utilities for
% working with networks

classdef NetworkObj < handle
    
properties                               
    nnode 
    nedge 
    dim 
    degrees % degree of each node
    maxdeg % max degree in network
    nodepos % position of each node
    edgenodes % nodes connected by each edge
    nodenodes % nodes connected to each node
    nodeedges % edges connected to each node
    edgelens % edge lengths
    
end

methods
    function NT = NetworkObj(fname)
        % create a network object
        % optionally, load data from file
        
        NT.nnode = 0; % number of nodes
        NT.nedge = 0; % number of edges
        NT.dim = 2; % spatial dimension        
        NT.maxdeg = 10; % maximum allowed degree              
        
        if (nargin>0)
            NT.loadNetwork(fname);
        end
    end
    
    function loadNetwork(NT,fname)
        % load a network from file
        
        fid = fopen(fname);
        
        tline = fgetl(fid);
        clear nodepos
        while (ischar(tline))
            % ignore comments and blank lines
            if (length(tline)==0 || tline(1)=='#'); tline = fgetl(fid); continue; end
            
            %disp(tline);
            
            % split line by whitespace
            words= strsplit(tline);
            
            if(strcmpi(words(1),'NODE'))
                % read in node position
                dim = nnz(cellfun(@(i) length(i), words(2:end))>0)-1;
                nums= cellfun(@(i) str2num(i), words(2:dim+2));
                nodepos(nums(1),:) = nums(2:end);
            elseif(strcmpi(words(1),'EDGE'))
                nums= cellfun(@(i) str2num(i), words(2:end));
                edgenodes(nums(1),:) = nums(2:3);                
                if (length(nums)>3)
                    % read in edge length if supplied
                    edgelens(nums(1)) = nums(4);
                end
            end
            tline = fgetl(fid);
        end
        if (exist('edgelens','var'))            
            edgelens(edgelens==0) = NaN;
            if (length(edgelens)<size(edgenodes,1))
                edgelens(end+1: size(edgenodes,1)) = NaN;
            end
        else
            edgelens = NaN*ones(size(edgenodes,1),1);
        end
        
        fclose(fid);
        
        NT.setupNetwork(nodepos,edgenodes,edgelens);
    end
    
    function setupNetwork(NT,nodepos,edgenodes,edgelens)
        % set up network structure starting with nodepos and edgenodes arrays
        % only
        % optionally, also supply list of edge lengths
        NT.nodepos = nodepos;
        NT.edgenodes = edgenodes;
        
        NT.dim = size(NT.nodepos,2);
        NT.nnode = size(NT.nodepos,1);
        NT.nedge = size(NT.edgenodes,1);
        
        % get degree and nodes attached to each node
        NT.degrees = zeros(NT.nnode,1);
        NT.nodenodes = zeros(NT.nnode,NT.maxdeg);
        NT.nodeedges = zeros(NT.nnode,NT.maxdeg);        
        if (exist('edgelens','var'))
            NT.edgelens = edgelens;
        else
            NT.edgelens = NaN*ones(NT.nedge,1);
        end
        
        for ec = 1:NT.nedge
            n1 = NT.edgenodes(ec,1); n2 = NT.edgenodes(ec,2);
            NT.degrees(n1) = NT.degrees(n1)+1;
            NT.nodenodes(n1,NT.degrees(n1)) = n2;
            NT.nodeedges(n1,NT.degrees(n1)) = ec;
            
            NT.degrees(n2) = NT.degrees(n2)+1;
            NT.nodenodes(n2,NT.degrees(n2)) = n1;
            NT.nodeedges(n2,NT.degrees(n2)) = ec;
            
            % calculate edge length if not supplied
            if isnan(NT.edgelens(ec))
                NT.edgelens(ec) = norm(NT.nodepos(n1,:) - NT.nodepos(n2,:));
            end
        end
        
        NT.maxdeg = max(NT.degrees);        
        NT.nodenodes = NT.nodenodes(:,1:NT.maxdeg);
        NT.nodeedges = NT.nodeedges(:,1:NT.maxdeg);        
        
    end
    
    function outputNetwork(NT,outfile,edgepaths,nodelabels)
        % output network structure to file
        % nodevals: arbitrary numerical value for each node
        % nodelabels: string associated with each node
        
        of = fopen(outfile,'w')
        fprintf(of,'%s\n','# file defining network structure')
        fprintf(of,'%s\n\n','# made with matlab NetworkObj output')
        
        fprintf(of,'%s \n','# list of node indices and xy positions and values')        
        nodefmtstring = ['NODE %d ' repmat(['%20.10f '],1,NT.dim+1) '\n'];
        nodelblfmtstring = ['NODE %d ' repmat(['%20.10f '],1,NT.dim+1) '%s \n'];
        %
        for pc = 1:size(NT.nodepos)
            
            if (~isfield(NT,'nodevals') || isempty(NT.nodevals))
                val = 0.0;              
            else
                val = NT.nodevals(pc);                
            end
            
            if (exist('nodelabels','var'))
                fprintf(of,nodelblfmtstring,pc, NT.nodepos(pc,:), val,nodelabels{pc});
            else
                fprintf(of,nodefmtstring,pc, NT.nodepos(pc,:), val);
            end
        end
        % edge information
        edgefmtstring = 'EDGE %d %d %d %20.10f\n';
        for ec = 1:size(NT.edgenodes)
            fprintf(of,edgefmtstring,ec, NT.edgenodes(ec,:), NT.edgelens(ec));
        end
               
        fclose(of)
    end
    
    function removeDoubleEdges(NT)
        %% get rid of all duplicate edges, so you only have one edge
        % between each pair of nodes
        badind = [];
        newedgenodes = [];
        edgenodes= NT.edgenodes;
        for ec1 = 1:NT.nedge
            duplicatefound = 0;
            for ec2 = ec1+1:NT.nedge
                if (edgenodes(ec1,1)==edgenodes(ec2,2) & edgenodes(ec1,2)==edgenodes(ec2,1))
                    %disp([ec1 ec2 edgenodes(ec1,:) edgenodes(ec2,:)])
                    duplicatefound = 1;
                elseif (edgenodes(ec1,1)==edgenodes(ec2,1) & edgenodes(ec1,2)==edgenodes(ec2,2))
                    %disp([ec1 ec2 edgenodes(ec1,:) edgenodes(ec2,:)])
                    duplicatefound=1;
                end
            end
            if ~duplicatefound
                newedgenodes = [newedgenodes; edgenodes(ec1,:)];
            end
        end
        NT.edgenodes = newedgenodes;
    end
    
    
    
     function [mapold2new] = keepNodes(NT,keepind)        
        % truncate network to only keep the node indices given by keepind
        % remove any edges not between nodes in the new network.
        % keep only largest connected component
        % mapnew2old(i) = index of new node i in the old node list
        % mapold2new(i) = index of old node i in new node list
        
        [newnodepos,newedgenodes,mapold2new] = truncateNetworkNodes(keepind,NT.nodepos,NT.edgenodes);
        
        NT.nodepos = newnodepos;
        NT.edgenodes = newedgenodes;
        
        NT.setupNetwork()
                    
        %if (~isempty(NT.nodevals)); NT.nodevals = NT.nodevals(keepind); end
        %if (~isempty(NT.edgevals)); NT.edgevals = NT.edgevals(keepind); end
     
     end
    
     function plotNetwork(NT,options)
         
         % default options
         opt.labelnodes = 0; % do not label nodes
         opt.nodeplotopt = {'MarkerSize',20}; % options to pass for plotting nodes
         opt.scatterplotopt = {'filled'}; % options for plotting color-coded nodes
         opt.markersize = 50; % marker size for colored scatter plot
         
         opt.edgeplotopt = {'LineWidth',1}; % options to pass for plotting edges
         
         opt.plotnodes = 1; % do plot the nodes
         opt.plotedges = 1; % do plot the edges
         opt.colornodes = NaN; % supply a list of values and color nodes according to the values
         
         if (exist('options','var')) % input options
             opt =copyStruct(options,opt,1);
         end
         
         if (opt.plotedges)
             p1 = NT.edgenodes(:,1); p2 = NT.edgenodes(:,2);
             edgepts = zeros(2,NT.nedge,NT.dim);
             for dc = 1:NT.dim
                 edgepts(:,:,dc) = [NT.nodepos(p1,dc)'; NT.nodepos(p2,dc)'];
             end
             if (NT.dim==2)
                 plot(edgepts(:,:,1),edgepts(:,:,2),'k',opt.edgeplotopt{:})
             else
                 plot3(edgepts(:,:,1),edgepts(:,:,2),edgepts(:,:,3),'k',opt.edgeplotopt{:})
             end
             
         end
         hold all
         
         if (opt.plotnodes)
             if (isnan(opt.colornodes))
                 if (NT.dim==3)
                     plot3(NT.nodepos(:,1),NT.nodepos(:,2),NT.nodepos(:,3),'b.',opt.nodeplotopt{:})
                 else
                     plot(NT.nodepos(:,1),NT.nodepos(:,2),'b.',opt.nodeplotopt{:})
                 end
             else
                 % colored nodes
                 if (NT.dim==3)
                     scatter3(NT.nodepos(:,1),NT.nodepos(:,2),NT.nodepos(:,3),opt.markersize,opt.colornodes,opt.scatterplotopt{:})
                 else
                     scatter(NT.nodepos(:,1),NT.nodepos(:,2),opt.markersize,opt.colornodes,opt.scatterplotopt{:})
                 end
             end
             axis equal
             hold all
         end
         
         if (opt.labelnodes)
             for pc = 1:length(NT.nodepos)
                 text(NT.nodepos(pc,1),NT.nodepos(pc,2),sprintf('%d',pc))
             end
         end
         axis equal
         hold off
     end
end

end