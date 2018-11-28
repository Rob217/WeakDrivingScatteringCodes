
%% LatticePositions - calculate bare lattice positions
function [nDips, rDips] = LatticePositions(nX,nY,nZ,shape)
% this calculates the bare lattice positions for a given lattice type
% i.e. sep = 1, every other parameter = 0


% if no inputs, create a custom lattice
if nargin == 0
    
    disp('plotting custom lattice')
    
    nX = 5; 
    nY = 1;
    nZ = 1;
    shape = 'h'; % equivalent to a=0.8 10x10 square lattice
    
end


%% different lattice types

if strcmp(shape,'s')==1;
    %% square lattice
        
    % total number of atoms
    nDips = nX*nY*nZ;

    rDips = zeros(nDips,3); % positions matrix
    nXs = -nX/2+0.5:1:nX/2-0.5; % indices of positions in x,y,z
    nYs = -nY/2+0.5:1:nY/2-0.5;
    nZs = -nZ/2+0.5:1:nZ/2-0.5;

    % shift at an angle
%     shift_atom = -(nDips-1)/2:1:(nDips-1)/2;
    shift_atom = zeros(nDips);
    
    Nth = 1; % Nth atom
    for ll = 1:nZ % loop over atoms
        for jj = 1:nY 
            for kk = 1:nX
                rDips(Nth,1) = nXs(kk)+shift_atom(Nth); % x coordinate 
                rDips(Nth,2) = nYs(jj); % y coordinate
                rDips(Nth,3) = nZs(ll); % z coordinate
                Nth = Nth+1; % goes to next atom
            end
        end
    end
    
elseif strfind(shape, 'PenroseUnit1') == 1
    % unit cell of Penrose tiling
    %
    %
    
    
elseif strfind(shape, 'Penrose1') == 1
    % quasicrystal based on Penrose tiling
    %
    % 'Penrose1_nIter' where nIter is the number of iterations in lattice construction algorithm
    
    underscoreInd = strfind(shape, '_');
    nSteps = str2num(shape(underscoreInd+1:end));
    
    % create ten simple triangles in a ring
    triangleList = zeros([10, 4]);
    triangleList(:,1) = 0;
    triangleList(:,2) = 0j;
    triangleList(:,3) = exp(1j*[0.5:9.5].*2.*pi/10);
    triangleList(:,4) = exp(1j*[1.5:10.5].*2.*pi/10);
    triangleList(2:2:end, [3,4]) = triangleList(2:2:end, [4,3]);

    goldenRatio = (1+sqrt(5))/2;
%     nSteps = 5;
    iStep = 1;    
    while iStep <= nSteps
        
        triangleListOld = triangleList;
        iTriangle = 1;
        
        for ii = 1:size(triangleListOld, 1)
            A = triangleListOld(ii, 2);
            B = triangleListOld(ii, 3);
            C = triangleListOld(ii, 4);
            if triangleListOld(ii, 1) == 0 % type 0
                P = A + (B - A)./goldenRatio;
                triangleListNew(iTriangle, :) = [0, C, P, B];
                triangleListNew(iTriangle + 1, :) = [1, P, C, A];
                iTriangle = iTriangle + 2;
            elseif triangleListOld(ii, 1) == 1 % type 1
                R = B + (C-B)./goldenRatio;
                Q = B + (A-B)./goldenRatio;
                triangleListNew(iTriangle, :) = [1, R, C, A];
                triangleListNew(iTriangle + 1, :) = [1, Q, R, B];
                triangleListNew(iTriangle + 2, :) = [0, R, Q, A];
                iTriangle = iTriangle + 3;
            end
        end
        triangleList = triangleListNew;
        iStep = iStep + 1;
    end
    
    rDips = triangleList(:, 2:4);
    rDips = [real(rDips(:)), imag(rDips(:)), zeros([3*size(rDips,1), 1])];
    ndips = size(rDips, 1);
    
    xdips = real(triangleList(:,2:4));
    ydips = imag(triangleList(:,2:4));
    
    % remove duplicates and scale by smallest separation such that nearest neighbor spacing = 1
    xdipsRep = repmat(rDips(:,1), [1, ndips]);
    ydipsRep = repmat(rDips(:,2), [1, ndips]);
    sep = sqrt((xdipsRep-xdipsRep').^2 + (ydipsRep-ydipsRep').^2);
    
    % find minimum sep not including duplicates
    minSep = min(sep(find(sep>1e-9)));
    
    % find groups of atoms that are on the same site
    sep = sep + 1e10.*eye(size(sep,1));
    [iDuplicate, jDuplicate] = find(sep<1e-9);
%     indsToInclude = [];
%     for ii = 1:length(iDuplicate);
%         iDuplicateInds = find(iDuplicate == iDuplicate(ii));
%         iDuplicateToAdd = max([iDuplicate(ii); jDuplicate(iDuplicateInds)]);
%         if any(find(iDuplicateToAdd==indsToInclude))==0
%             indsToInclude = [indsToInclude, iDuplicateToAdd];
%         end
%     end    
    
    indsMat = zeros([size(sep)]);
    indsDuplicate = [iDuplicate, jDuplicate];
    sub2indDuplicate = sub2ind(size(sep), iDuplicate, jDuplicate);
    indsMat(sub2indDuplicate) = sub2indDuplicate;
    allInds = 1:ndips;
    indsMat1 = repmat(allInds, [ndips, 1]);
    indsMat2 = permute(indsMat1, [2,1]);
    indsMat1Sub = indsMat1;
    indsMat2Sub = indsMat2;
    indsMat1Sub(sub2indDuplicate) = 0;
    indsMat2Sub(sub2indDuplicate) = 0;
    indsMat1 = indsMat1 - indsMat1Sub;
    indsMat2 = indsMat2 - indsMat2Sub;
    max1 = max(indsMat1, [], 1)';
    max2 = max(indsMat1, [], 2);    
    indsToInclude = unique(max(max1, max2))';
    
    nonDuplicatedInds = setdiff(allInds, unique(iDuplicate));
    indsToInclude = [nonDuplicatedInds; indsToInclude];
    
    % add atoms, only including first atom among duplicates
    rDipsDuplicatesRemoved = rDips(indsToInclude,:);
    rDips = rDipsDuplicatesRemoved./minSep; % scale s.t. nearest neighbor spacing = 1 
    ndips = size(rDips, 1);
    
    plotPenrose = 0;
    if plotPenrose == 1
        figure(1001)
        clf
        hold on, box on
        shape0Inds = find(triangleList(:,1) == 0);
        shape1Inds = find(triangleList(:,1) == 1);
        for ii = 1:length(shape0Inds)
            ind0 = shape0Inds(ii);
            fill(xdips(ind0,[1,2,3,1]), ydips(ind0,[1,2,3,1]), 'r')
        end
        for jj = 1:length(shape1Inds)
            ind1 = shape1Inds(jj);
            fill(xdips(ind1,:), ydips(ind1,:), 'b')
        end
%         scatter(rDips(:,1), rDips(:,2), [], 'k', 'filled')
        daspect([1,1,1])
    end  
    
    
elseif strfind(shape, 'Penrose2') == 1
    % quasicrystal based on Penrose tiling
    %
    % 'Penrose2_nIter' where nIter is the number of iterations in lattice construction algorithm
    
    underscoreInd = strfind(shape, '_');
    nSteps = str2num(shape(underscoreInd+1:end));
    
    % create ten simple triangles in a ring
    triangleList = zeros([10, 4]);
    triangleList(:,1) = 0;
    triangleList(:,2) = exp(1j*[0.5:9.5].*2.*pi/10);
    triangleList(:,3) = 0j;
    triangleList(:,4) = exp(1j*[1.5:10.5].*2.*pi/10);
    triangleList(2:2:end, [2,4]) = triangleList(2:2:end, [4,2]);

    goldenRatio = (1+sqrt(5))/2;
    iStep = 1;    
    while iStep <= nSteps
        
        triangleListOld = triangleList;
        iTriangle = 1;
        
        for ii = 1:size(triangleListOld, 1)            

            A = triangleListOld(ii, 2);
            B = triangleListOld(ii, 3);
            C = triangleListOld(ii, 4);
            if triangleListOld(ii, 1) == 0 % type 0
                Q = A + (B - A)./goldenRatio;
                R = B + (C - B)./goldenRatio;
                triangleListNew(iTriangle, :) = [0, C, A, R];
                triangleListNew(iTriangle + 1, :) = [1, R, Q, B];
                triangleListNew(iTriangle + 2, :) = [0, Q, A, R];
                iTriangle = iTriangle + 3;
            elseif triangleListOld(ii, 1) == 1 % type 1
                P = C + (A - C)./goldenRatio;
                triangleListNew(iTriangle, :) = [1, B, P, A];
                triangleListNew(iTriangle + 1, :) = [0, P, C, B];
                iTriangle = iTriangle + 2;
            end
        end
        triangleList = triangleListNew;
        iStep = iStep + 1;
    end
    
    rDips = triangleList(:, 2:4);
    rDips = [real(rDips(:)), imag(rDips(:)), zeros([3*size(rDips,1), 1])];
    ndips = size(rDips, 1);
    
    xdips = real(triangleList(:,2:4));
    ydips = imag(triangleList(:,2:4));
    
    % remove duplicates and scale by smallest separation such that nearest neighbor spacing = 1
    xdipsRep = repmat(rDips(:,1), [1, ndips]);
    ydipsRep = repmat(rDips(:,2), [1, ndips]);
    sep = sqrt((xdipsRep-xdipsRep').^2 + (ydipsRep-ydipsRep').^2);
    
    % find minimum sep not including duplicates
    minSep = min(sep(find(sep>1e-9)));
    
    % find groups of atoms that are on the same site
    sep = sep + 1e10.*eye(size(sep,1));
    [iDuplicate, jDuplicate] = find(sep<1e-9);  
    
    indsMat = zeros([size(sep)]);
    indsDuplicate = [iDuplicate, jDuplicate];
    sub2indDuplicate = sub2ind(size(sep), iDuplicate, jDuplicate);
    indsMat(sub2indDuplicate) = sub2indDuplicate;
    allInds = 1:ndips;
    indsMat1 = repmat(allInds, [ndips, 1]);
    indsMat2 = permute(indsMat1, [2,1]);
    indsMat1Sub = indsMat1;
    indsMat2Sub = indsMat2;
    indsMat1Sub(sub2indDuplicate) = 0;
    indsMat2Sub(sub2indDuplicate) = 0;
    indsMat1 = indsMat1 - indsMat1Sub;
    indsMat2 = indsMat2 - indsMat2Sub;
    max1 = max(indsMat1, [], 1)';
    max2 = max(indsMat1, [], 2);    
    indsToInclude = unique(max(max1, max2))';
    
    nonDuplicatedInds = setdiff(allInds, unique(iDuplicate));
    indsToInclude = [nonDuplicatedInds; indsToInclude];
    
    % add atoms, only including first atom among duplicates
    rDipsDuplicatesRemoved = rDips(indsToInclude,:);
    rDips = rDipsDuplicatesRemoved./minSep; % scale s.t. nearest neighbor spacing = 1 
    ndips = size(rDips, 1);
    
    plotPenrose = 1;
    if plotPenrose == 1 && nargin == 0
        figure(1001)
        clf
        hold on, box on
        shape0Inds = find(triangleList(:,1) == 0);
        shape1Inds = find(triangleList(:,1) == 1);
        for ii = 1:length(shape0Inds)
            ind0 = shape0Inds(ii);
            fill(xdips(ind0,[1,2,3,1]), ydips(ind0,[1,2,3,1]), 'r')
        end
        for jj = 1:length(shape1Inds)
            ind1 = shape1Inds(jj);
            fill(xdips(ind1,:), ydips(ind1,:), 'b')
        end
%         scatter(rDips(:,1), rDips(:,2), [], 'k', 'filled')
        daspect([1,1,1])
    end 
    
    

elseif strfind(shape, 'q1') == 1
    % quasicrystal based on 4 overlapping optical lattice beams forming two square lattices rotated with respect to each other
    %
    % q1_theta_radius_Vmin_boundary where theta = angle between square lattices and radius = relative radius of inclusion and Vmin = minimum potential to include lattice site and boundary = 'circ' or 'square'
    underscoreInds = strfind(shape, '_');
    theta = str2num(shape(underscoreInds(1)+1:underscoreInds(2)-1));
    radius = str2num(shape(underscoreInds(2)+1:underscoreInds(3)-1));
    Vmin = str2num(shape(underscoreInds(3)+1:underscoreInds(4)-1));
    boundary = shape(underscoreInds(4)+1:end);
    
    theta = theta*pi/180; % convert to radians
    
    x = 1.1*radius.*linspace(-1,1,20*radius);
    y = 1.1*radius.*linspace(-1,1,20*radius);
    [X,Y] = meshgrid(x,y);

%     k = 4/sqrt(0.55);
    k = 5.35;
    E1 = sin(k.*X) + sin(k.*Y);
    E2 = sin(cos(theta).*k.*X-sin(theta).*k.*Y) + sin(sin(theta).*k.*X + cos(theta).*k.*Y);
    Vtot = abs(E1+E2).^2;

%     minIndsX = [];
%     minIndsY = [];
%     for iX = 2:length(x)-1
%         for iY = 2:length(y)-1  
%             VtotSquare = Vtot(iX-1:iX+1, iY-1:iY+1);
%             [maxVal, maxInd] = max(VtotSquare(:));
%             if maxInd == 5 && Vtot(iX,iY)>Vmin
%                 minIndsX = [minIndsX, iX];
%                 minIndsY = [minIndsY, iY];
%             end
%         end
%     end
%     minIndsOldMethod = sort(sub2ind(size(Vtot), minIndsX, minIndsY));
    
    VTotSquare = zeros([length(x)-2, length(y)-2, 3, 3]);
    for ii = 1:3
        for jj = 1:3            
            VTotSquare(:,:,ii,jj) = Vtot(ii:(end-3+ii), jj:(end-3+jj));
        end
    end
    [maxVal, maxIndDim1] = max(VTotSquare(:,:,:), [], 3);
    [maxIndsX1, maxIndsY1] = find(maxIndDim1==5);
    maxInds = sub2ind(size(Vtot), maxIndsX1+1, maxIndsY1+1);
    indsOverVMin = find(Vtot>Vmin);
    maxIndsOverVMin = intersect(maxInds, indsOverVMin);
    [maxIndsX, maxIndsY] = ind2sub(size(Vtot), maxIndsOverVMin);
    xLattice = x(maxIndsX);
    yLattice = y(maxIndsY);

    if strcmp(boundary, 'circ')
        rPosMag = sqrt(xLattice.^2 + yLattice.^2);
        xLattice = xLattice(find(rPosMag<=radius));
        yLattice = yLattice(find(rPosMag<=radius));
    end
    
    rDips(:,1) = xLattice(:);
    rDips(:,2) = yLattice(:);
    rDips(:,3) = 0;
    nDips = size(rDips,1);

    
    
elseif strfind(shape,'rect')==1
    %% rectangular array, 'shape' = 'rect_x_y_z'
    %  where x,y,z are the relative dimension sizes
    
        % total number of atoms
    nDips = nX*nY*nZ;

    rDips = zeros(nDips,3); % positions matrix
    nXs = -nX/2+0.5:1:nX/2-0.5; % indices of positions in x,y,z
    nYs = -nY/2+0.5:1:nY/2-0.5;
    nZs = -nZ/2+0.5:1:nZ/2-0.5;

    % shift at an angle
%     shift_atom = -(nDips-1)/2:1:(nDips-1)/2;
    shift_atom = zeros(nDips);
    
    Nth = 1; % Nth atom
    for ll = 1:nZ % loop over atoms
        for jj = 1:nY 
            for kk = 1:nX
                rDips(Nth,1) = nXs(kk)+shift_atom(Nth); % x coordinate 
                rDips(Nth,2) = nYs(jj); % y coordinate
                rDips(Nth,3) = nZs(ll); % z coordinate
                Nth = Nth+1; % goes to next atom
            end
        end
    end
    
    % squash/extend out different dimensions
    dims = find(shape=='_');
    dim1 = str2num(shape(dims(1)+1:dims(2)-1));
    dim2 = str2num(shape(dims(2)+1:dims(3)-1));
    dim3 = str2num(shape(dims(3)+1:end));
    rDips(:,1) = rDips(:,1).*dim1;
    rDips(:,2) = rDips(:,2).*dim2;
    rDips(:,3) = rDips(:,3).*dim3;
    
elseif strcmp(shape,'slide')==1;
    %% slide one atom wrt to the other, keeping the first fixed in place at (x,y,z) = (0,0,0)
    rDips(1,:) = [0,0,0];
    Nth = 1;
    for ii = 1:nX
        for jj = 1:nY
            rDips(Nth,1) = (ii-1);
            rDips(Nth,2) = (jj-1);
            rDips(Nth,3) = 0;
            Nth = Nth+1;
        end
    end
    nDips = size(rDips,1);
    
    
elseif strcmp(shape,'shiftedplane')==1;
    %% multiple square planes each shifted w.r.t. the one in front
    % total number of atoms
    nDips = nX*nY*nZ;

    rDips = zeros(nDips,3); % positions matrix
    nXs = -nX/2+0.5:1:nX/2-0.5; % indices of positions in x,y,z
    nYs = -nY/2+0.5:1:nY/2-0.5;
    nZs = -nZ/2+0.5:1:nZ/2-0.5;

    % shift at an angle
%     shift_atom = -(nDips-1)/2:1:(nDips-1)/2;
    shift_atom = zeros(nDips);
    
    Nth = 1; % Nth atom
    for ll = 1:nZ % loop over atoms
        for jj = 1:nY 
            for kk = 1:nX
                rDips(Nth,1) = nXs(kk)+shift_atom(Nth); % x coordinate 
                rDips(Nth,2) = nYs(jj); % y coordinate
                rDips(Nth,3) = nZs(ll); % z coordinate                
                % shift each plane with respect to those in front and behind
                if nZ>1
                    oddeven = mod(ll,2);
                    if oddeven == 1
                        rDips(Nth,1) = rDips(Nth,1) + 0.25;
                        rDips(Nth,2) = rDips(Nth,2) + 0.25;
                    else
                        rDips(Nth,1) = rDips(Nth,1) - 0.25;
                        rDips(Nth,2) = rDips(Nth,2) - 0.25;
                    end
                end
                
                Nth = Nth+1; % goes to next atom

            end
        end
    end    
    

    
    
    
elseif strcmp(shape,'t')==1;   
    %% triagonal lattice (overall shape = hexagon)
    % nY controls height of array - has to be even, otherwise rounded down
    % nX controls width
    
%     nDips = nX*nY-nY^2/4+0.25;
%     rDips = zeros(nDips,3);
    
    Nth = 1;
    % loop over y axis
    for y_ind = 1:nY/2+1/2;
        x_lim = nX/2-y_ind/2; % limit for x range
        % loop over x axis
        for x_ind = -x_lim:1:x_lim;
            rDips(Nth,1) = x_ind;
            rDips(Nth,2) = (y_ind-1)*sqrt(3)/2; 
            Nth = Nth + 1;
            if y_ind ~= 1
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = -(y_ind-1)*sqrt(3)/2;
                Nth = Nth + 1;
            end
        end
    end  
    nDips = size(rDips,1);
    rDips(:,3) = 0;  % add z-dimension
    
    
elseif strcmp(shape,'t2')==1;   
    %% triagonal lattice (overall shape = triangle)
    % nY controls height of array
        
    Nth = 1;
    % loop over y axis
    for y_ind = 1:nY
        for x_ind = 1:(nX+1-y_ind)
            rDips(Nth,1) = x_ind - 1 +(y_ind-1)/2;
            rDips(Nth,2) = sqrt(3)/2*(y_ind-1);
            Nth = Nth + 1;
        end
    end
    
    % recentralise
    rDips(:,1) = rDips(:,1) - (nX-1)/2;
    rDips(:,2) = rDips(:,2) - (nY-1)*sqrt(3)/6;
    
    % add z component
    nDips = size(rDips,1);
    rDips(:,3) = zeros(1,length(nDips));
        
    
elseif strcmp(shape,'h')==1;
    %% if hexagonal lattice
    % doesn't matter what nY is

    nDips = 6*nX^2;
    rDips = zeros(nDips,3);
    
    Nth = 1;
    
    % position 1
    ii = nX;
    for jj = 1:ii
        % jj = row index
        x_lim = ii-jj/2;
        for x_val = -x_lim:1:x_lim
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = (1.5*jj-1);
            Nth = Nth + 1;
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = -(1.5*jj-1);
            Nth = Nth + 1;
        end
    end

    % position 2
    for jj = 1:ii
        x_lim = (ii-jj/2-1/2);
        for x_val = -x_lim:1:x_lim
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = (1.5*jj-1/2);
            Nth = Nth + 1;
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = -(1.5*jj-1/2);
            Nth = Nth + 1;
        end
    end
  
elseif strcmp(shape,'h2')==1;
    %% if hexagonal lattice
    % same as h but with additional row added to break degenerate states

    nDips = 6*nX^2+nX;
    rDips = zeros(nDips,3);
    
    Nth = 1;
    
    % position 1
    ii = nX;
    for jj = 1:ii
        % jj = row index
        x_lim = ii-jj/2;
        for x_val = -x_lim:1:x_lim
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = (1.5*jj-1);
            Nth = Nth + 1;
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = -(1.5*jj-1);
            Nth = Nth + 1;
        end
    end

    % position 2
    for jj = 1:ii
        x_lim = (ii-jj/2-1/2);
        for x_val = -x_lim:1:x_lim
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = (1.5*jj-1/2);
            Nth = Nth + 1;
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = -(1.5*jj-1/2);
            Nth = Nth + 1;
        end
    end
    [NN,yesno] = find(rDips(:,2)==max(rDips(:,2)));
    for Nind = 1:length(NN)
        rDips(Nth,1) = rDips(NN(Nind),1);
        rDips(Nth,2) = rDips(NN(Nind),2)+1;
        Nth = Nth + 1;
    end
    
    
elseif strcmp(shape,'h30')==1;
    %% if hexagonal lattice
    % doesn't matter what nY is

    nDips = 6*nX^2;
    rDips = zeros(nDips,3);
    
    Nth = 1;
    
    % position 1
    ii = nX;
    for jj = 1:ii
        % jj = row index
        x_lim = ii-jj/2;
        for x_val = -x_lim:1:x_lim
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = (1.5*jj-1);
            Nth = Nth + 1;
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = -(1.5*jj-1);
            Nth = Nth + 1;
        end
    end

    % position 2
    for jj = 1:ii
        x_lim = (ii-jj/2-1/2);
        for x_val = -x_lim:1:x_lim
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = (1.5*jj-1/2);
            Nth = Nth + 1;
            rDips(Nth,1) = x_val.*sqrt(3);
            rDips(Nth,2) = -(1.5*jj-1/2);
            Nth = Nth + 1;
        end
    end
    
    % rotate through 30 degrees (added 8.4.16)
    tt = pi/6;
    for ii = 1:size(rDips,1)
        rDips(ii,:) = ([[cos(tt),-sin(tt),0];[sin(tt),cos(tt),0];[0,0,1]]*rDips(ii,:)')';
    end
    
       
elseif strcmp(shape,'k')==1;
    %% kagome lattice (hexagonal in shape)
    % nX must be a multiple of 4n-1 where n=1,2,3,... i.e. % nX=3,7,11,...
    % nY makes no difference
    
    if mod(nX,2) == 1 % odd
        nDips = 9*nX^2/4 - 2*nX + 3/4;
    else % even
        nDips = 3*nX*(3*nX-2)/4;
    end
    rDips = zeros(nDips,3);
    
    Nth = 1;
    

    % even nDips
    if mod(nX,2) == 0

        % position 1
        for y_ind = 1:nX/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 

    % odd nDips
    else
        % position 1
        for y_ind = 1:nX/2+1/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2-1/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 
    end
    
    
elseif strcmp(shape,'k3')==1;
    %% kagome lattice built up from unit cell
    % nY makes no difference
    
    Nth = 1;
    
    x0 = [0];
    y0 = [0];
    for ii = 1:length(x0)
        for jj = 1:length(y0)
            rDips(Nth,1) = x0(ii); % atom 1
            rDips(Nth,2) = y0(jj);
            rDips(Nth+1,1) = x0(ii) + 1; % atom 2
            rDips(Nth+1,2) = y0(jj);
            rDips(Nth+2,1) = x0(ii) + 1/2; % atom 3
            rDips(Nth+2,2) = y0(jj) - sqrt(3)/2;
            rDips(Nth+3,1) = x0(ii) - 1; % atom 4
            rDips(Nth+3,2) = y0(jj);
            rDips(Nth+4,1) = x0(ii) - 1/2; % atom 5
            rDips(Nth+4,2) = y0(jj) + sqrt(3)/2;
            Nth = Nth+5;
        end
    end
    
%     x0 = [0:2:2*(nX-1)]
%     y0 = [0:2:2*(nY-1)]
%     for ii = 1:length(x0)
%         for jj = 1:length(y0)
%             rDips(Nth,1) = x0(ii);
%             rDips(Nth,2) = y0(jj);
%             rDips(Nth+1,1) = x0(ii);
%             rDips(Nth+1,2) = y0(jj)+1;
%             rDips(Nth+2,1) = x0(ii)+1;
%             rDips(Nth+2,2) = y0(jj)+1;
%             Nth = Nth+3;
%         end
%     end
    
    nDips = size(rDips,1);
    rDips(:,3) = zeros(1,length(nDips));
    
    
elseif strcmp(shape,'k30')==1;
    %% kagome lattice built up from unit cell - rotate by 30deg
    % nY makes no difference
    
    Nth = 1;
    
    x0 = [0];
    y0 = [0];
    for ii = 1:length(x0)
        for jj = 1:length(y0)
            rDips(Nth,1) = x0(ii); % atom 1
            rDips(Nth,2) = y0(jj);
            rDips(Nth+1,1) = x0(ii) + 1; % atom 2
            rDips(Nth+1,2) = y0(jj);
            rDips(Nth+2,1) = x0(ii) + 1/2; % atom 3
            rDips(Nth+2,2) = y0(jj) - sqrt(3)/2;
            rDips(Nth+3,1) = x0(ii) - 1; % atom 4
            rDips(Nth+3,2) = y0(jj);
            rDips(Nth+4,1) = x0(ii) - 1/2; % atom 5
            rDips(Nth+4,2) = y0(jj) + sqrt(3)/2;
            Nth = Nth+5;
        end
    end
    
%     x0 = [0:2:2*(nX-1)]
%     y0 = [0:2:2*(nY-1)]
%     for ii = 1:length(x0)
%         for jj = 1:length(y0)
%             rDips(Nth,1) = x0(ii);
%             rDips(Nth,2) = y0(jj);
%             rDips(Nth+1,1) = x0(ii);
%             rDips(Nth+1,2) = y0(jj)+1;
%             rDips(Nth+2,1) = x0(ii)+1;
%             rDips(Nth+2,2) = y0(jj)+1;
%             Nth = Nth+3;
%         end
%     end
    
    nDips = size(rDips,1);
    rDips(:,3) = zeros(1,length(nDips));
    
    % rotate through 30 degrees (added 8.4.16)
    tt = pi/6;
    for ii = 1:size(rDips,1)
        rDips(ii,:) = ([[cos(tt),-sin(tt),0];[sin(tt),cos(tt),0];[0,0,1]]*rDips(ii,:)')';
    end
    
    
elseif strcmp(shape,'k4')==1;
    %% kagome lattice built up from unit cell
    % nY makes no difference
    
    Nth = 1;
    
    rDips(1:2,1) = [-1,1];
    rDips(1:2,2) = [0,0];
    rDips(3:6,1) = [-1.5:1.5];
    rDips(3:6,2) = sqrt(3)/2.*ones(1,4);
    rDips(7:10,1) = rDips(3:6,1);
    rDips(7:10,2) = -rDips(3:6,2);
    rDips(11,1:2) = [0,sqrt(3)];
    rDips(12,1:2) = [0,-sqrt(3)];
    
    nDips = size(rDips,1);
    rDips(:,3) = zeros(1,length(nDips));

    
elseif strcmp(shape,'k5')==1;
    %% squashed kagome - angle between atoms now 54deg ('magic angle')
    % nX must be a multiple of 4n-1 where n=1,2,3,... i.e. % nX=3,7,11,...
    % nY makes no difference
    
    if mod(nX,2) == 1 % odd
        nDips = 9*nX^2/4 - 2*nX + 3/4;
    else % even
        nDips = 3*nX*(3*nX-2)/4;
    end
    rDips = zeros(nDips,3);
    
    Nth = 1;
    

    % even nDips
    if mod(nX,2) == 0

        % position 1
        for y_ind = 1:nX/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 

    % odd nDips
    else
        % position 1
        for y_ind = 1:nX/2+1/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2-1/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 
    end
    
    
elseif strcmp(shape,'k6')==1;
    %% kagome restricted to square lattice (nX,nY) in size
   
    % first make standard kagome lattice then remove atoms outside this

    Nth = 1;
    nX_square = nX;
    nY_square = nY;
    nX = ceil(nX + 10); % make much bigger than square will be

    % even nDips
    if mod(nX,2) == 0

        % position 1
        for y_ind = 1:nX/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 

    % odd nDips
    else
        % position 1
        for y_ind = 1:nX/2+1/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2-1/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 
    end
    
    
    % now remove atoms outside square
%     r = sqrt(abs(rDips(:,1)).^2 + abs(rDips(:,2)).^2 + abs(rDips(:,3)).^2);
%     rDips = rDips(rDips(:,1)<nX)
    rDips = rDips(find(abs(rDips(:,1))<nX_square),:); % fit to nX
    rDips = rDips(find(abs(rDips(:,2))<nY_square),:); % fit to nY

    nDips = size(rDips,1);
    rDips(:,3) = zeros(nDips,1);

    
elseif strcmp(shape,'k60')==1;
    %% kagome lattice - same as k5 but rotated 60 degrees
    % nX must be a multiple of 4n-1 where n=1,2,3,... i.e. % nX=3,7,11,...
    % nY makes no difference
    
    if mod(nX,2) == 1 % odd
        nDips = 9*nX^2/4 - 2*nX + 3/4;
    else % even
        nDips = 3*nX*(3*nX-2)/4;
    end
    rDips = zeros(nDips,3);
    
    Nth = 1;
    

    % even nDips
    if mod(nX,2) == 0

        % position 1
        for y_ind = 1:nX/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 

    % odd nDips
    else
        % position 1
        for y_ind = 1:nX/2+1/2;
            x_lim = nX-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:2:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1)*sqrt(3); 
                Nth = Nth + 1;
                if y_ind ~= 1
                    rDips(Nth,1) = x_ind;
                    rDips(Nth,2) = -(y_ind-1)*sqrt(3);
                    Nth = Nth + 1;
                end
            end
        end 

        % position 2
        for y_ind = 1:nX/2-1/2;
            x_lim = nX-1/2-y_ind; % limit for x range
            % loop over x axis
            for x_ind = -x_lim:1:x_lim;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (y_ind-1/2)*sqrt(3); 
                Nth = Nth + 1;
                rDips(Nth,1) = x_ind;
                rDips(Nth,2) = (1/2-y_ind)*sqrt(3);
                Nth = Nth + 1;
            end
        end 
    end
    
    theta = 60*pi/180; % 60 degrees
    rsep = sqrt(rDips(:,1).^2 + rDips(:,2).^2);
    theta1 = atan2(rDips(:,2),rDips(:,1));
    rDips_rot(:,1) = cos(theta+theta1).*rsep;
    rDips_rot(:,2) = sin(theta+theta1).*rsep;
    rDips_rot(:,3) = rDips(:,3);
    rDips = rDips_rot;
    
    
elseif strcmp(shape,'r')==1
    %% random array inside cubic volume of nX*nY*nZ
    nDips = nX*nY*nZ;
    rDips = zeros(nDips,3);
    rDips(:,1) = (rand(nDips,1)-0.5).*nX;
    rDips(:,2) = (rand(nDips,1)-0.5).*nY;
    rDips(:,3) = (rand(nDips,1)-0.5).*nZ;
    
    
elseif shape(1) == 'G' 
    %% Gaussian cloud
    nDips = str2num(shape(2:end)); % G followed by nDips

    rDips = zeros(nDips,3);
    rDips(:,1) = normrnd(zeros(nDips,1),nX);
    rDips(:,2) = normrnd(zeros(nDips,1),nY);
    rDips(:,3) = normrnd(zeros(nDips,1),nZ);
 

elseif shape(1:4) == 'Rcub' 
    %% random array inside cubic volume of nDips taken from shape = 'RnDips'
    nDips = str2num(shape(6:end));
    rDips = zeros(nDips,3);
    rDips(:,1) = (rand(nDips,1)-0.5).*nX;
    rDips(:,2) = (rand(nDips,1)-0.5).*nY;
    rDips(:,3) = (rand(nDips,1)-0.5).*nZ;
    
    
elseif shape(1:4) == 'Rcyl' 
    %% random array inside cylindrical volume of nDips = nX
    % Rcyl_dz
    % where dz = thickness w.r.t. dR, radius
    % dR should be set to 1 to enforce N2D = 1 for a = 1
    nDips = nX; % nX,nY,nZ = nDips,1,1
    NN = 1;
    rDips = zeros(nDips,3);
    
    % extract dR and dz from shape
    dz = str2num(shape(5:end));
    dR = sqrt(nDips/pi); % normalise s.t. N2D = 1
    dz = dz.*dR;
    
    while NN <= nDips
        insideCyl = 0; % change truth value to >1 when value is inside cylinder 
        while insideCyl == 0
            rDips(NN,:) = (rand(1,3,1)-0.5).*2.*[dR,dR,dz];
            if sqrt(sum(rDips(NN,1:2).^2))<=dR
                insideCyl = 1;
            end
        end
        NN = NN + 1;
    end

    
    
elseif any(strfind(shape,'RZCyl'))
    %% random array inside cylindrical volume
    % 'shape' has the form RZCyl_dR_dz
    % cylindrical volume has dimensions length dz and radius dR
    % atom number, nDips = nX

    nDips = nX; 
    NN = 1; % counter
    rDips = zeros(nDips,3);
    
    % extract dz and dR from 'shape'
    underscoreLocs = strfind(shape,'_'); % find indices where '_' appears
    dR = str2num(shape(underscoreLocs(1)+1:underscoreLocs(2)-1));
    dz = str2num(shape(underscoreLocs(2)+1:end));
    
    % loop over atoms, creating atom locations
    while NN <= nDips
        insideCyl = 0; % change truth value to >1 when value is inside cylinder 
        while insideCyl == 0
            rDips(NN,:) = (rand(1,3,1)-0.5).*[2.*dR,2.*dR,dz];
            if sqrt(sum(rDips(NN,1:2).^2))<=dR % check whether position is within cylinder
                insideCyl = 1;
            end
        end
        NN = NN + 1;
    end

elseif strcmp(shape,'Lieb')==1
    %% old full 2D Lieb lattice
    
    % first make full lattice:
    xdips = -(nX):(nX);
    ydips = -(nY):(nY);
    [X,Y] = meshgrid(xdips,ydips);
    rDips(:,1) = X(:);
    rDips(:,2) = Y(:);
    
    % then make secondary lattice to remove from full lattice
    xdips2 = -(nX-1):2:(nX-1);
    ydips2 = -(nY-1):2:(nY-1);
    [X2,Y2] = meshgrid(xdips2,ydips2);
    rDips2(:,1) = X2(:);
    rDips2(:,2) = Y2(:);
    
    % find common elements and remove them
    [commonVals,ind1,ind2] = intersect(rDips,rDips2,'rows');
    rDips(ind1,:) = [];
    
    % add z dimension
    rDips(:,3) = [0];
    
    nDips = size(rDips,1);
    
    
elseif strcmp(shape,'Lieb1D')==1
    %% 1D Lieb lattice (alternating nearest neighbour spacing of a, 2a, a, 2a, etc.)
    
    if length(find([nX,nY,nZ]>1)) > 1
        error('Only one of [nX,nY,nZ] can be >1 for Lieb1D')
    end
    
    [nDips, dimInd] = max([nX,nY,nZ]);    
        
    leftSites = 1:ceil(nDips/2); % create dipole positions
    rightSites = 1:floor(nDips/2);
    dipPos = sort([3.*leftSites, 3.*rightSites + 1]);
    
    meanDipPos = mean(dipPos);
    dipPos = dipPos - meanDipPos; % centre dipole positions
    
    rDips = zeros([nDips, 3]);
    rDips(:, dimInd) = dipPos;

    
    
elseif strcmp(shape,'Lieb2D')==1
    %% 2D Lieb lattice
    
    % basis
    basisAtoms = [[0,0,0]; [1,0,0]; [0,1,0]];
    nBasisAtoms = size(basisAtoms, 1);
    
    % create unit lattice
    xUnitLattice = 2.*(1:nX);
    yUnitLattice = 2.*(1:nY);
    zUnitLattice = 2.*(1:nZ);    
    [X,Y,Z] = meshgrid(xUnitLattice, yUnitLattice, zUnitLattice); 
    nUnitLattice = length(X(:)); % number of unit lattice sites        
    
    iDip = 0;
    rDips = zeros([nBasisAtoms.*nUnitLattice, 3]);
    for iUnitLattice = 1:nUnitLattice
        for iBasisAtom = 1:nBasisAtoms;
            iDip = iDip + 1;
            rDips(iDip, 1) = X(iUnitLattice) + basisAtoms(iBasisAtom, 1);
            rDips(iDip, 2) = Y(iUnitLattice) + basisAtoms(iBasisAtom, 2);
            rDips(iDip, 3) = Z(iUnitLattice) + basisAtoms(iBasisAtom, 3);
        end
    end
    
    nDips = iDip;
    
    % reposition to be centered on (0,0,0)
    meanPos = mean(rDips, 1);
    meanPosRep = repmat(meanPos, [nDips, 1]);
    rDips = rDips - meanPosRep;

    
    
elseif strcmp(shape,'Lieb3D')==1
    %% 3D Lieb lattice
    
    % basis
    basisAtoms = [[0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]];
    nBasisAtoms = size(basisAtoms, 1);
    
    % create unit lattice
    xUnitLattice = 2.*(1:nX);
    yUnitLattice = 2.*(1:nY);
    zUnitLattice = 2.*(1:nZ);    
    [X,Y,Z] = meshgrid(xUnitLattice, yUnitLattice, zUnitLattice); 
    nUnitLattice = length(X(:)); % number of unit lattice sites        
    
    iDip = 0;
    rDips = zeros([nBasisAtoms.*nUnitLattice, 3]);
    for iUnitLattice = 1:nUnitLattice
        for iBasisAtom = 1:nBasisAtoms;
            iDip = iDip + 1;
            rDips(iDip, 1) = X(iUnitLattice) + basisAtoms(iBasisAtom, 1);
            rDips(iDip, 2) = Y(iUnitLattice) + basisAtoms(iBasisAtom, 2);
            rDips(iDip, 3) = Z(iUnitLattice) + basisAtoms(iBasisAtom, 3);
        end
    end
    
    nDips = iDip;
    
    % reposition to be centered on (0,0,0)
    meanPos = mean(rDips, 1);
    meanPosRep = repmat(meanPos, [nDips, 1]);
    rDips = rDips - meanPosRep;
    
    
elseif strcmp(shape,'Lieb3DFull')==1
    %% 3D Lieb lattice
    
    % basis
    basisAtoms = [[0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]];
    nBasisAtoms = size(basisAtoms, 1);
    
    % create unit lattice
    xUnitLattice = 2.*(1:nX+1);
    yUnitLattice = 2.*(1:nY+1);
    zUnitLattice = 2.*(1:nZ+1);    
    [X,Y,Z] = meshgrid(xUnitLattice, yUnitLattice, zUnitLattice); 
    nUnitLattice = length(X(:)); % number of unit lattice sites        
    
    iDip = 0;
%     rDips = zeros([nBasisAtoms.*nUnitLattice, 3]);
    for iUnitLattice = 1:nUnitLattice
        for iBasisAtom = 1:nBasisAtoms;
            iDip = iDip + 1;
            rDips(iDip, 1) = X(iUnitLattice) + basisAtoms(iBasisAtom, 1);
            rDips(iDip, 2) = Y(iUnitLattice) + basisAtoms(iBasisAtom, 2);
            rDips(iDip, 3) = Z(iUnitLattice) + basisAtoms(iBasisAtom, 3);
        end
    end
    
    % remove outer edges to smooth off boundaries
    cutoffX = find(rDips(:,1) > xUnitLattice(end) + 0.1);
    rDips(cutoffX,:) = [];
    cutoffY = find(rDips(:,2) > yUnitLattice(end) + 0.1);
    rDips(cutoffY,:) = [];
    cutoffZ = find(rDips(:,3) > zUnitLattice(end) + 0.1);
    rDips(cutoffZ,:) = [];
    
    nDips = size(rDips, 1);
    
    % reposition to be centered on (0,0,0)
    meanPos = mean(rDips, 1);
    meanPosRep = repmat(meanPos, [nDips, 1]);
    rDips = rDips - meanPosRep;
    
    
elseif strcmp(shape,'SqOct')==1
    %% square-octagon semiregular lattice
    
    % first start by making base lattice
    dim = 1+1/sqrt(2);
    baseX = [-(nX-1):2:(nX-1)].*dim;
    baseY = [-(nY-1):2:(nY-1)].*dim;
    [X,Y] = meshgrid(baseX,baseY);
    baseXY(:,1) = X(:);
    baseXY(:,2) = Y(:);
    
    % at each lattice point, build 5 squares of dipoles
    square_points = [[1,1];[1,-1];[-1,1];[-1,-1]]; % corners on a square
    centres(1,:) = [0,0];
    centres(2:5,:) = square_points;
    centres = centres.*dim;
    for ii = 1:5 % loop over centres
        squares((ii-1)*4+1:ii*4,:) = repmat(centres(ii,:),4,1) + 0.5.*square_points;
    end
    
    % apply squares to each base lattice point
    for jj = 1:size(baseXY,1)
        rDips((jj-1)*20+1:jj*20,:) = repmat(baseXY(jj,:),20,1) + squares;
    end
    
    % remove double values
    rDips = unique(rDips,'rows');
        
    
    rDips(:,3) = [0];
    nDips = size(rDips,1)
    
    
elseif shape(1:5) == 'sRoty' 
    %% square lattice rotated in y axis 
    theta = str2num(shape(6:end)); % angle, in degrees
    theta = theta.*pi./180; % degrees to radians

    % total number of atoms
    nDips = nX*nY*nZ;

    rDips = zeros(nDips,3); % positions matrix
    nXs = -nX/2+0.5:1:nX/2-0.5; % indices of positions in x,y,z
    nYs = -nY/2+0.5:1:nY/2-0.5;
    nZs = -nZ/2+0.5:1:nZ/2-0.5;

    % shift at an angle
%     shift_atom = -(nDips-1)/2:1:(nDips-1)/2;
    shift_atom = zeros(nDips);
    
    Nth = 1; % Nth atom
    for ll = 1:nZ % loop over atoms
        for jj = 1:nY 
            for kk = 1:nX
                rDips(Nth,1) = nXs(kk)+shift_atom(Nth); % x coordinate 
                rDips(Nth,2) = nYs(jj); % y coordinate
                rDips(Nth,3) = nZs(ll); % z coordinate
                Nth = Nth+1; % goes to next atom
            end
        end
    end
    
    % rotate around y axis
    rDips_rot(:,1) = rDips(:,1);
    rDips_rot(:,2) = rDips(:,2).*cos(theta) - rDips(:,3).*sin(theta);
    rDips_rot(:,3) = rDips(:,2).*sin(theta) + rDips(:,3).*cos(theta);
    rDips = rDips_rot;
    

elseif shape(1:5) == 'sRotz'     
    %% square lattice rotated in z axis 
    theta = str2num(shape(6:end)); % angle, in degrees
    theta = theta.*pi./180; % degrees to radians

    % total number of atoms
    nDips = nX*nY*nZ;

    rDips = zeros(nDips,3); % positions matrix
    nXs = -nX/2+0.5:1:nX/2-0.5; % indices of positions in x,y,z
    nYs = -nY/2+0.5:1:nY/2-0.5;
    nZs = -nZ/2+0.5:1:nZ/2-0.5;

    % shift at an angle
%     shift_atom = -(nDips-1)/2:1:(nDips-1)/2;
    shift_atom = zeros(nDips);
    
    Nth = 1; % Nth atom
    for ll = 1:nZ % loop over atoms
        for jj = 1:nY 
            for kk = 1:nX
                rDips(Nth,1) = nXs(kk)+shift_atom(Nth); % x coordinate 
                rDips(Nth,2) = nYs(jj); % y coordinate
                rDips(Nth,3) = nZs(ll); % z coordinate
                Nth = Nth+1; % goes to next atom
            end
        end
    end
    
    % rotate around z axis
    rDips_rot(:,1) = rDips(:,1).*cos(theta) - rDips(:,2).*sin(theta);
    rDips_rot(:,2) = rDips(:,1).*sin(theta) + rDips(:,2).*cos(theta);
    rDips_rot(:,3) = rDips(:,3);
    rDips = rDips_rot;

    
end


% plot custom lattice
if nargin == 0
    
    figure(1000)
    clf
    hold on, box on
    
    scatter3(rDips(:,1),rDips(:,2),rDips(:,3), 100, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [180,199,231]./255, 'LineWidth', 2)  
    daspect([1,1,1])
    XLIM = get(gca,'XLim');
%     xlim(1.1.*XLIM)
%     view([30,20])
    
    % calculate density
    calculateDensity = 0;
    if calculateDensity == 1
        max(rDips(:,1))
        radiusDensity = 0.8*max(rDips(:,1));
        rDipsMag = sqrt(sum(abs(rDips).^2, 2));
        indsInsideRadius = find(rDipsMag <= radiusDensity);
        nInsideRadius = length(indsInsideRadius);
        density = nInsideRadius./(pi*radiusDensity^2)
        scatter3(rDips(indsInsideRadius,1),rDips(indsInsideRadius,2),rDips(indsInsideRadius,3), 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [255,130,130]./255, 'LineWidth', 2)  
    end
    
end

nDips = size(rDips, 1);

end
