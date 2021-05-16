%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% makeH.m
%
% Construct the Hine's Matrix for the cell with morphology encoded by md
%

function [H,SA,psa] = makeH(md,Ra,As)

% preallocate memory for connection strength matrix
S = spalloc(md.count.L,md.count.L,md.count.L*3-4);

soma_index = md.count.L;     % soma is indexed last

sacnt = 0;

psa = [];

for ii = 1:md.count.N    % fill in S for each branch
    
    % save low_index and high index for block for current branch 
    low_index  = sum(md.grid.comps(1:ii-1))+1;
    high_index = sum(md.grid.comps(1:ii));
       
    dx = md.grid.hstep(ii);     % save dx for this branch
       
    a = md.interp.radii{ii}(:);     % save radii for current branch
    
    %plot(a)
    %keyboard
    
    psa = [psa; 2*pi*dx.*a(1:end-1)];
    
    main_diag  = zeros(length(a),1);
    sub_diag   = zeros(length(a),1);
    super_diag = zeros(length(a),1);
     
    % spatial finite difference at interior branch nodes
    for jj = 2:length(a)-2   
        aterm = a(jj)*(a(jj+1)-a(jj-1))/(2*dx^2);
        main_diag(jj)    = -2*a(jj)^2/dx^2 / (2*Ra*a(jj));
        super_diag(jj+1) = (a(jj)^2/dx^2 + aterm) / (2*Ra*a(jj));
        sub_diag(jj-1)   = (a(jj)^2/dx^2 - aterm) / (2*Ra*a(jj));
    end
    
    % build chunck of S matrix that describes all Inner Branch
    % Compartment connections along branch ii
    tridiag_S = spdiags([sub_diag(:) main_diag(:) super_diag(:)],[-1 0 1],...
        md.grid.comps(ii),md.grid.comps(ii));

    % insert
    S(low_index:high_index,low_index:high_index) = tridiag_S;
    
    % do most distal node
    if ~isempty(find(md.leaf == ii))                            % if branch is a leaf
        S(low_index,low_index)   = -a(1)^2 / (2*Ra*dx^2*a(1));
        S(low_index,low_index+1) = a(1)^2 / (2*Ra*dx^2*a(1));
    else                                                        % if branch is a parent 
        parent = ii;
        
        children = find(md.parents == ii);      % get children branch indices
        
        % gives index where parent attaches to daughter branch
        parent_index = sum(md.grid.comps(1:parent-1))+1;
        
        % get nodes, radii, and step size of children
        child_index = zeros(length(children),1);
        for c = 1:length(children)
            child_index(c)  = sum(md.grid.comps(1:children(c)));
            child_dx(c)     = md.grid.hstep(children(c));
            child_radius(c) = md.interp.radii{children(c)}(end-1);  % get at next to last node
        end      
        
        parent_dx = md.grid.hstep(parent);      % step dx in parent branch
        
        % radius of connecting parent branch: get one node beyond last!
        parent_radii = md.interp.radii{parent}(2);
        
        % set S matrix entries
        for c = 1:length(children)
            cfactor = child_radius(c)^2 /(2*Ra*child_dx(c)*parent_dx*parent_radii);           
            S(parent_index,parent_index)   = S(parent_index,parent_index) - cfactor;
            S(parent_index,child_index(c)) = cfactor;
        end
        pfactor = parent_radii^2 /(2*Ra*parent_dx^2*parent_radii);
        
        S(parent_index,parent_index)   = S(parent_index,parent_index) - pfactor;
        S(parent_index,parent_index+1) = pfactor;
        
    end
        
    parent = md.parents(ii);  % store parent of current branch 
    
    % now handle connection of branch ii to it's predeccesor
    if(parent > 0)
        parent_index = sum(md.grid.comps(1:parent-1))+1;
        
        % radius of connecting parent branch
        parent_radius = md.interp.radii{parent}(1);
        
        rad_idx = max(length(a)-2,1);
     
        aterm = a(end-1)*(parent_radius-a(rad_idx))/(2*dx^2);
                
        S(high_index,high_index)   = -2*a(end-1)^2/dx^2 / (2*Ra*a(end-1));
        S(high_index,high_index-1) = (a(end-1)^2/dx^2 - aterm) / (2*Ra*a(end-1));
        S(high_index,parent_index) = (a(end-1)^2/dx^2 + aterm) / (2*Ra*a(end-1));
        
    else       % for soma ordering is the last, soma_index

        rad_idx = max(length(a)-2,1);

        aterm = a(end-1)*(a(end)-a(rad_idx))/(2*dx^2);

        S(high_index,high_index)   = -2*a(end-1)^2/dx^2 / (2*Ra*a(end-1));
        S(high_index,high_index-1) = (a(end-1)^2/dx^2 - aterm) / (2*Ra*a(end-1));
        S(high_index,soma_index)   = (a(end-1)^2/dx^2 + aterm) / (2*Ra*a(end-1));
        
        sacnt = sacnt + 1;
        SA(sacnt,1) = 2*pi*a(end)*dx;

        % compartment touching soma
        S(soma_index,high_index) = pi*a(end-1)^2/(Ra*dx*As);

        % soma compartment
        S(soma_index,soma_index) = S(soma_index,soma_index) - pi*a(end-1)^2 /(Ra*dx*As);
    end
end

psa = [psa; 2*pi*dx*As];

H = -S;