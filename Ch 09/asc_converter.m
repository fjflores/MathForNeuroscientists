%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% Converts Neurolucida .asc files into raw data for assembler
% Processes data for special cases
%
% usage: neur_converter(filename)  (ex. test.asc)
%
% ex. neur_coverter(test.asc)  
%
% output: 
%   nc.treedata   : cell array containing the connectivity data
%   nc.celldata   : cell array containing cell morphology
%   nc.somadata   : struct array containing information about the soma
%   nc.lendata    : cell array containing cumulative lengths of each segment
%   nc.roots      : array containing the segment indices of the nc.roots
%   OPENF         : handle of the opened file
%
% contains code from convert.m by Steve Cox and Jay Raol
%
% Nan Xiao
% neur_converter.m

% updated 6/29/05

% 'Hines' ordering version

function [nc, OPENF] = asc_converter(filename)
OPENF = fopen(filename,'rt');

if OPENF == -1 
    disp(['failed to open file ' filename]);
    return
end

fprintf(1,'\n');
disp(['Reading .asc data for ' filename])  

fprintf(1,'\n');
fprintf(1,'storing coordinates');

% Loops through the lines and Find Soma Data
word = ''; 		            % word on a line

while ~strcmp(word,'CellBody')  
    line = fgetl(OPENF);
    word = getword(line);
end

line = fgetl(OPENF);
word = getword(line);

isNumber = 0;
nc.somadata.x = [];
nc.somadata.y = [];
nc.somadata.z = [];

% loops through cellbody data
while ~strcmp(word,'Endofcontour'),     
    
    num = 0;
    number = '';

    for j=1:length(line),

        % goes through each character in a line
        
        if line(j) ~= ' ' & line(j) ~= '(' & num <= 3,
            
            % the character is not a space or '(' then it is a number
            % and stored in number
            
            number = [number line(j)];
            isNumber = 1;
        
        elseif isNumber,
            
            % if a number has just been read in and the character
            % is not a number, then reset
            
            num = num + 1;
            isNumber = 0;
            storeword = 1;
        
        end

        if num == 1 & storeword,          %store x,y,z information
            
            nc.somadata.x = [nc.somadata.x str2num(number)];
            number = '';
            storeword = 0;
        
        elseif num == 2 & storeword,
            
            nc.somadata.y = [nc.somadata.y str2num(number)];
            number = '';
            storeword = 0;
        
        elseif num == 3 & storeword,
            
            nc.somadata.z = [nc.somadata.z str2num(number)];
            number = '';
            storeword = 0;
        
        end

    end	% end for j

    line = fgetl(OPENF);
    word = getword(line);

end  % end while

% Calculate Soma Properties

pcount = 1;

for j = 1:length(nc.somadata.x),
    
    for i = 1:length(nc.somadata.y),
    
        diams(pcount) = sqrt((nc.somadata.y(j)-nc.somadata.y(i))^2 ...
                        +(nc.somadata.x(j)-nc.somadata.x(i))^2);   
        pcount = pcount + 1;    
    
    end

end

nc.somadata.d = max(diams);
nc.somadata.l = nc.somadata.d;

% calculates area of soma
nc.somadata.As = abs(area(nc.somadata.x,nc.somadata.y));	

% begin to read in dendrite information

line = fgetl(OPENF);           
word = getword(line);

while ~strcmp(word,'Dendrite')
    
   line = fgetl(OPENF);
   word = getword(line);
   
end

segcount = 0;                   

while ~feof(OPENF),            % loop until the end of file
    line = fgetl(OPENF);
    isNumber = 0;
    goloop = 1;
    
    % check segment data
    
    % if the line of data has 'R', as in the R-1-2, etc...
    
    if findstr(line,'; R') > 0, 
        
        % then we have found a segment chunk
        
        goloop = 1;             
        segcount = segcount + 1;

        segx = [];
        segy = [];
        segz = [];
        segd = [];
        connstr = [];
        
    else

        goloop = 0;
        
    end

    coordcount = 1;

    % loops through each segment
    
    while goloop == 1,          

        num = 0;
        number = '';

        % goes through each character in a line

        for j=1:length(line),	

            if line(j) ~= ' ' & line(j) ~= '(' & ...
               line(j) ~= ')' & num <= 4,
                
                % the character is not a space or '(' 
                % then it is a number
                % and stored in number
                
                number = [number line(j)];
                isNumber = 1;
            
            elseif isNumber
                
                % if a number has just been read in and the character
                % is not a number, then reset
                
                num = num + 1;
                isNumber = 0;
                storeword = 1;
            
            end

            %store x coordinates in x vector, y in y vector

            if num == 1 & storeword,
                
                segx = [segx str2num(number)];
                number = '';
                storeword = 0;
                
            elseif num == 2 & storeword,
                
                segy = [segy str2num(number)];
                number = '';
                storeword = 0;
                
            elseif num == 3 & storeword,
                
                segz = [segz str2num(number)];
                number = '';
                storeword = 0;
                
            elseif num == 4 & storeword,
                
                segd = [segd str2num(number)];
                number = '';
                storeword = 0;
                
            end

        end

        % get the connectivity pattern
        
        if coordcount == 1,          
            
            j = 3;
            tempstr = line(length(line)-j);
            
            while tempstr ~= 'R',
                
                tempstr = line(length(line)-j);
                j = j+1;
                
            end

            connstr = line(length(line)-j:length(line)-3);
        end

        line = fgetl(OPENF);

        coordcount = coordcount + 1;

        % check if we are still inside a segment chunk

        if findstr(line,'R') > 0
        
            goloop = 1;
        
        else

            goloop = 0;
        
        end

    end  % end while

    segdata = [segx; segy; segz; segd]';
    nc.celldata{segcount} = segdata;
    conndata{segcount} = connstr;

    if mod(segcount,5) == 0
        fprintf(1,'.');
    end
end

% close file
fclose(OPENF);

for j =1:length(conndata)

    if conndata{j}(1) == ' ',
    
        conndata{j}(1) = [];
    
    end
end

% connections: nc.roots to nc.somadata

treecount = 1;

for j = 1:segcount,

    if conndata{j} == 'R'
 
        nc.roots(treecount) = j;
        treecount = treecount + 1;
        
    end
    
end

% write the branch connections 
% (using string search, not base 3 conversion)

nc.roots(length(nc.roots)+1) = segcount;
a = 1;
nc.treedata{segcount} = [];

fprintf(1,'\n');
fprintf(1,'Writing connections');

while a < length(nc.roots),
    
    j = nc.roots(a)+1;
    
    if a == length(nc.roots)-1,
        cond = segcount+1;
    else    
        cond = nc.roots(a+1);
    end
    
    while j < cond
    
        b = 1;
        lastdash = 0;
        
        while ~lastdash,
    
            if conndata{j}(length(conndata{j})-b) == '-',
                
                lastdash = length(conndata{j})-b;
                
            else

                b = b+1;
                
            end
            
        end
       
        k = nc.roots(a);
        
        while k < nc.roots(a+1),
        
            if length(conndata{k}) == ...
                      length(conndata{j}(1:lastdash-1)) &...
                      conndata{k} == conndata{j}(1:lastdash-1),
                  
                nc.treedata{k} = [nc.treedata{k} j];
                
                nc.celldata{j}(2:size(nc.celldata{j},1)+1,:) = ...
                  nc.celldata{j}(1:size(nc.celldata{j},1),:);
              
                nc.celldata{j}(1,:) = ...
                  nc.celldata{k}(size(nc.celldata{k},1),:);
            end

            k = k + 1;
            
        end
        
        j = j + 1;
    end
    
    a = a + 1;
    
    fprintf(1,'.');
end    

% Remove possible redundant points

for j = 1:segcount
    
    numslink = 0;
    k = 1;
    %disp(['j = ' num2str(j)])
    
    while k < size(nc.celldata{j},1)-numslink 
        
        if (nc.celldata{j}(k,1) == nc.celldata{j}(k+1,1) & ...
            nc.celldata{j}(k,2) == nc.celldata{j}(k+1,2)) | ...
           (nc.celldata{j}(k,1) == nc.celldata{j}(k+1,1) & ...
            nc.celldata{j}(k,2) == nc.celldata{j}(k+1,2) & ...
            nc.celldata{j}(k,3) == nc.celldata{j}(k+1,3)),
            
            nc.celldata{j}(k,:) = [];
            numslink = numslink + 1;
            
        end

        k = k + 1;
        
    end
   
end

% Check for "single-forked" branches 
% (weird case, but it's in some of the traces)

numsingle = 0;

% index shows whether or not the segment is a lone child

singleindex = zeros(1,segcount);        

% combines child and parent into one segment

for j = 1:segcount                      
    
    % checks for lone children
    
    if length(nc.treedata{j}) == 1         

        % if so, store number of points in parent
        
        parentlength = size(nc.celldata{j},1);                 
        
        % store number of points in child
        
        childlength = size(nc.celldata{nc.treedata{j}(1)},1);     
        
        % merge parent with child
        
        nc.celldata{j}(parentlength+1:parentlength+childlength-1,:) ...    
                     = nc.celldata{nc.treedata{j}(1)}(2:childlength,:);
                 
        % memorize that this segment is a lone child
                     
        singleindex(nc.treedata{j}(1)) = 1;

        % increment number of lone children
        
        numsingle = numsingle + 1;              

    end
end

% re-adjust non-single child segments

if numsingle > 0                      

    newlength = size(nc.celldata,2) - numsingle;
    
    % initialize new cell array
    newnc.celldata{newlength} = [];                
    
    % initialize new tree array
    newnc.treedata{newlength} = [];                
    
    % For each lone child (that will be deleted from the cell array) ...
    % every index in tree data that is greater than the index of the lone...
    % child must be decremented by one.
   
    % traverse the segments
    
    for j = 1:length(singleindex)         
        
        % check if the branch is a lone child
        
        if singleindex(j) == 1                  
            
            % traverse the tree data
            
            for k = 1:size(nc.treedata,2)          
                
                % check if segment has 2 children (normal case)
                
                if length(nc.treedata{k}) == 2     
                    
                    % if so, check if children index are larger than ... 
                    % that of the current segment
                    
                    if nc.treedata{k}(1) > j       
                        
                        nc.treedata{k}(1) = nc.treedata{k}(1) - 1;    
                        
                    end
                    
                    if nc.treedata{k}(2) > j
                        
                        nc.treedata{k}(2) = nc.treedata{k}(2) - 1;
                        
                    end
                    
                end
                
            end
            
            % traverse nc.roots
            
            for i = 1:length(nc.roots)             
                
                % decrement root indices along the same lines
                
                if nc.roots(i) > j                 
                    
                    nc.roots(i) = nc.roots(i) - 1;    
                    
                end

            end

        end
        
    end
    
    k = 1;
    i = 1;
    
    % Remove the lone child from the original arrays by copying ...
    % the normal segments into new arrays
    
    % traverse segments
    
    for j = 1:segcount                      
    
        % copy only if segment is not a lone child
        
        if singleindex(j) == 0              
            
            newnc.celldata{k} = nc.celldata{j};
            k = k + 1;
            
        end

        % copy only if tree data does not point to lone child
        
        if length(nc.treedata{j}) ~= 1         
            
            newnc.treedata{i} = nc.treedata{j};
            i = i + 1;
            
        end
        
    end    
    
    % new fixed arrays
    
    nc.celldata = newnc.celldata;                    
    nc.treedata = newnc.treedata;     
    
    % new fixed segcount
    
    segcount = size(nc.celldata,2);            
    
    clear newnc.celldata;
    clear newnc.treedata;
    
end

% fudge with trifurcating branches
% Add a minute segment (one point from a trifurcation child)
% Change connections so that new segment connects to two of the...
% parent's children; the parent still is connected to its first ...
% child but is now also connected to the new segment.

for j = 1:segcount
    
    if length(nc.treedata{j}) == 3
    
        newindex = segcount + 1;
        nc.celldata{newindex}(1:2,:) = nc.celldata{j}(1:2,:);
        nc.celldata{j}(1,:) = [];
        nc.treedata{newindex}(1) = nc.treedata{j}(2);
        nc.treedata{newindex}(2) = nc.treedata{j}(3);
        nc.treedata{j}(2) = newindex;
        nc.treedata{j}(3) = [];
        segcount = segcount + 1;
        
    end
    
end

% reorder

nc = neur_reorder(nc);

% Write down branch lengths

fprintf(1,'\n');
fprintf(1,'Writing branch lengths')

nc.lendata{segcount} = [];

for j = 1:segcount 
    
    nc.lendata{j} = ...
         getlength(nc.celldata{j}(1:size(nc.celldata{j},1),1),...
                   nc.celldata{j}(1:size(nc.celldata{j},1),2));
               
    if mod(j,5) == 0
        
        fprintf(1,'.');
        
    end
    
    nc.segsize(j) = size(nc.celldata{j},1);
end

nc.type = 'asc';

treecount = treecount - 1;
fprintf(1,'\n');
fprintf(1,'\n');
disp(['Trees: ' num2str(treecount)]);
disp(['Segments: ' num2str(segcount)]);

% nc.roots

nc.roots = nc.roots(1:length(nc.roots)-1);

% write filename

nc.filename = filename;

function arc = getlength(x,y)

arc = 0;
sum = 0;

for i=2:length(x),
    
   sum = sum + sqrt((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2);
   arc = [arc sum];
   
end
%pause

return

function word = getword(line)

word = '';

for j=1:length(line),
    
   if line(j) == '"',
       
      word = '';
      return;
      
   end

   if isletter(line(j)),
       
      word = [word line(j)];
      
   end
   
end

return

function  a = area(x,y)

% AREA(X,Y)  Calculates the area of the 2-dimensional polygon
%      formed by vertices with coordinates vectors  X and Y.
%   Kirill K. Pankratov, kirill@plume.mit.edu.
%   April 20, 1994

x = [x(:); x(1)]; y = [y(:); y(1)];
lx = length(x);
a = abs( (x(2:lx)-x(1:lx-1))' * (y(1:lx-1)+y(2:lx)) )/2;

return

% perform a Hines reordering of the data

function newnc = neur_reorder(nc)

k = 1;
parent = 0;
children = 1;

for ii = 1:length(nc.treedata)
 
    if size(nc.treedata{ii},2) > 0
        
    % find parents
        
        parent(k) = ii;
        children(k) = nc.treedata{ii}(1);
        k = k + 1;
        parent(k) = ii;
        children(k) = nc.treedata{ii}(2);
        k = k + 1;
        
    end
    
end

pardex = zeros(1,length(nc.treedata));
pardex(children) = parent;

% loop through branches to find depth

for ii = 1:length(nc.treedata)
    
    jj = ii;
    depth(ii) = 1;
    jj = pardex(jj);
    
    while jj
        
        jj = pardex(jj);
        depth(ii) = depth(ii) + 1;
        
    end
    
end

% rearrange the branches by depth (highest to lowest)

m = 1;

maxdepth = max(depth);

newnc.somadata = nc.somadata;
%newnc.filename = nc.filename;

%nc.roots = find(depth == 1);

for jj = maxdepth:-1:1

    cdex = find(depth == jj);
    
    for kk = 1:length(cdex)
    
        ii = cdex(kk);
        
        newnc.celldata{m} = nc.celldata{ii};
        newnc.treedata{m} = nc.treedata{ii};
        %newnc.lendata{m} = nc.lendata{ii};
        %newnc.segsize(m) = nc.segsize(ii);
        
        newdex(ii) = m;
        
        m = m + 1;
        
    end    
    
end

% take care of new connections
% and reorder the roots

for ii = 1:length(nc.treedata)
    
    if size(newnc.treedata{ii},2) > 0
    
        newnc.treedata{ii}(1) = newdex(newnc.treedata{ii}(1));
        newnc.treedata{ii}(2) = newdex(newnc.treedata{ii}(2));
        
    end
    
end

for ii = 1:length(nc.roots)
   
    newnc.roots(ii) = newdex(nc.roots(ii));
    
end

% reversal of fortune

for ii = 1:length(nc.treedata)
    
    newnc.celldata{ii} = flipud(newnc.celldata{ii});
    
end   

m = 1;

% store depths,leaves, and parents

newnc.parent = zeros(1,length(newnc.treedata));

for ii = 1:length(depth)

    newnc.depth(newdex(ii)) = depth(ii);

    if(length(newnc.treedata{ii}) == 0)

        newnc.leafs(m) = ii;
        m = m + 1;

    else

        newnc.parent(newnc.treedata{ii}) = ii;

    end

end