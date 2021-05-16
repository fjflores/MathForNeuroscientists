%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% makemd.m
% 
% Generates the morphology data structure based on the neuron
% characteristics structure, and sets some parameters.
%
% Modified 9/3/08 by Anthony Kellems:
%          Includes simParams input to take in default parameters
% Modified 12/17/07 by Anthony Kellems:
%          Combines swc and asc checks into this file (better portability,
%          easier code changing)
% Modified 9/26/07 by Anthony Kellems:
%          md.dist_to_soma : cumulative distance from branch to soma
% Modified 7/16/09 by S Cox: removed simParams argument
%
% INPUTS:   nc        = neuron characteristics structure
%           numcomps  = desired compartment size (um)
%           %simParams = structure of parameters (from init_default_sim_params)
%
% OUTPUT:   md = morphology data structure

function md = makemd(nc,numcomps)

if strcmp(nc.filename(end-3:end),'.asc')
    fileType = 'asc';
else
    fileType = 'swc';
end

md.count.N = length(nc.celldata);   % total number of branches

md.grid.h_desired = numcomps;

md.parents = zeros(md.count.N,1);
md.par = [];

p = 1;
paryes = 0;

% traverses through all of the branches to interpolate radii

for j = 1:md.count.N         
    
    % finds physical length of the current branch
    md.grid.segL(j) = nc.lendata{j}(end)*1e-4;        % convert from microm to cm
    
    % md.interpolates radii for the current segment
    
    % number of compartments needed
    md.grid.comps(j) = max(5, round(nc.lendata{j}(end)/numcomps) );       
    
    % compartment length calculated
    md.grid.hstep(j) = md.grid.segL(j)/md.grid.comps(j);            
    
    % check for proper radius interpolation (asc uses diameter, swc uses radius)
    if strcmp(fileType,'asc')
        md.interp.radii{j} = interp1(1e-4*nc.lendata{j},...
            1e-4*(nc.celldata{j}(1:nc.segsize(j),4))',...
            0:md.grid.hstep(j):md.grid.segL(j))/2;
    else
        % no divide by 2, because SWC format specifies radius, not diameter
        md.interp.radii{j} = interp1(1e-4*nc.lendata{j},...
            1e-4*(nc.celldata{j}(1:nc.segsize(j),4))',...
            0:md.grid.hstep(j):md.grid.segL(j));   % convert from microm to cm
    end
    
    md.interp.size(j) = length(md.interp.radii{j});
    
    % regular spatial mesh 
    md.grid.x{j} = linspace(0,md.grid.segL(j),md.interp.size(j));
                             
    % find parents
    if size(nc.treedata{j},2) > 0
        md.par(p) = j;
        p = p + 1;
        paryes = 1;
        
        for k = 1:length(nc.treedata{j})
            md.parents(nc.treedata{j}(k)) = j;
        end
    end    
    
    % get those leaves
    leaves = ones(1,md.count.N);
    leaves(md.par) = 0;
    md.leaf = find(leaves == 1);

    clear leaves

    % find largest radius
    md.maxrad(j) = max(nc.celldata{j}(1:nc.segsize(j),4));
    
end
disp(['Total compartments: ' num2str(1+sum(md.grid.comps))])

md.count.L = sum(md.interp.size)-size(md.interp.size,2)+1;  % total number of compartments

if paryes 
    md.count.P = size(md.par,2);    % total number of parents
else
    md.count.P = 0;
end

% find the largest radius
md.maxrad = max(md.maxrad);

