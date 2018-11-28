%% function for calculating the dipole-dipole interaction matrix
function Gij = CalculateInteractionMatrixGij(rDips,k,dims,e0,polBasis,justReIm)
    nDips = size(rDips,1); % number of atoms
    nDims = length(dims);
    
    % create u=kr as 'cell' vector
    x = rDips(:,1);
    y = rDips(:,2);
    z = rDips(:,3);
    
    nDips = size(rDips,1);

    % array method (fastest)
    X = repmat(x,[1,nDips]);
    Y = repmat(y,[1,nDips]);
    Z = repmat(z,[1,nDips]);
    ux = k.*(X-transpose(X));
    uy = k.*(Y-transpose(Y));
    uz = k.*(Z-transpose(Z));
    nu = sqrt(ux.^2 + uy.^2 + uz.^2);
    uu = cat(3,ux,uy,uz); % cat links ux uy uz together into one cell uu

    delta_xy = kron(ones(nDips,nDips)-eye(nDips),eye(3));
    delta_ij = kron(ones(nDips,nDips)-eye(nDips),ones(3));
    g1 = kron( (3/2).*exp(1j.*nu)./(nu.^3).*(nu.^2+1j.*nu-1) , ones(3,3));
    g2 = kron( (3/2).*exp(1j.*nu)./(nu.^3).*(-nu.^2-3j.*nu+3) , ones(3,3));
    xy_r2 = zeros(3*nDips);
    xy_r2(1:3:end,1:3:end) = ux.*ux./nu.^2;
    xy_r2(1:3:end,2:3:end) = ux.*uy./nu.^2;
    xy_r2(1:3:end,3:3:end) = ux.*uz./nu.^2;
    xy_r2(2:3:end,1:3:end) = uy.*ux./nu.^2;
    xy_r2(2:3:end,2:3:end) = uy.*uy./nu.^2;
    xy_r2(2:3:end,3:3:end) = uy.*uz./nu.^2;
    xy_r2(3:3:end,1:3:end) = uz.*ux./nu.^2;
    xy_r2(3:3:end,2:3:end) = uz.*uy./nu.^2;
    xy_r2(3:3:end,3:3:end) = uz.*uz./nu.^2;

    G = g1.*delta_xy + g2.*delta_ij.*xy_r2;
    G(find(isnan(G))) = 0;
    G = G.*k^3/6/pi;  
    
    if strcmp(justReIm, 'Re')
        disp('Just taking Re(Gij)')
        G = real(G);
    elseif strcmp(justReIm, 'Im')
        disp('Just taking Im(Gij)')
        G = imag(G);
    end
    
    % rotate from xyz to +-z basis
    if strcmp(polBasis, '+-z')
        nDips = size(rDips,1);
        Gij_pmz = zeros([3*nDips,3*nDips]);

        Gij_pmz(1:3:end,1:3:end) = (G(1:3:end,1:3:end) + G(2:3:end,2:3:end))./2;
        Gij_pmz(2:3:end,2:3:end) = (G(1:3:end,1:3:end) + G(2:3:end,2:3:end))./2;
        Gij_pmz(3:3:end,3:3:end) = G(3:3:end,3:3:end);
%         Gij_pmz(1:3:end,2:3:end) = 1j.*G(1:3:end,2:3:end) + (G(1:3:end,1:3:end) - G(2:3:end,2:3:end))./2; % OLD - incorrect
%         Gij_pmz(2:3:end,1:3:end) = -1j.*G(1:3:end,2:3:end) + (G(1:3:end,1:3:end) - G(2:3:end,2:3:end))./2;
%         Gij_pmz(1:3:end,3:3:end) = (G(1:3:end,3:3:end) + 1j.*G(2:3:end,3:3:end))./sqrt(2);
%         Gij_pmz(2:3:end,3:3:end) = (G(1:3:end,3:3:end) - 1j.*G(2:3:end,3:3:end))./sqrt(2);
%         Gij_pmz(3:3:end,1:3:end) = (G(1:3:end,3:3:end) - 1j.*G(2:3:end,3:3:end))./sqrt(2);
%         Gij_pmz(3:3:end,2:3:end) = (G(1:3:end,3:3:end) + 1j.*G(2:3:end,3:3:end))./sqrt(2);        
        Gij_pmz(1:3:end,2:3:end) = -1j.*G(1:3:end,2:3:end) + (G(1:3:end,1:3:end) - G(2:3:end,2:3:end))./2; % NEW - correct
        Gij_pmz(2:3:end,1:3:end) = 1j.*G(1:3:end,2:3:end) + (G(1:3:end,1:3:end) - G(2:3:end,2:3:end))./2;
        Gij_pmz(1:3:end,3:3:end) = (G(1:3:end,3:3:end) - 1j.*G(2:3:end,3:3:end))./sqrt(2);
        Gij_pmz(2:3:end,3:3:end) = (G(1:3:end,3:3:end) + 1j.*G(2:3:end,3:3:end))./sqrt(2);
        Gij_pmz(3:3:end,1:3:end) = (G(1:3:end,3:3:end) + 1j.*G(2:3:end,3:3:end))./sqrt(2);
        Gij_pmz(3:3:end,2:3:end) = (G(1:3:end,3:3:end) - 1j.*G(2:3:end,3:3:end))./sqrt(2);
        G = Gij_pmz;  

    end

    dimInds = dims(1):3:nDips*3; % indices for extracting correct dimensions from Gij
    for iDim = 2:nDims
        dimInds = sort([dimInds,dims(iDim):3:nDips*3]);
    end
    Gij = G(dimInds,dimInds);
    
    % Gij = Gij; % Ruostekoski model
    Gij = Gij./e0; % Bettles model

end