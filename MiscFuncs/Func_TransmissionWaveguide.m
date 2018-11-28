%% function for calculating electric field collected by a waveguide
function [PL, PdTot, PTot] = TransmissionWaveguide(EFieldModel, polVec, polBasis, propagationDirection, c, lambda0, w0, EF, kappa, zLensIn, radiusLensIn, nRhoLensIn, zLensOut, radiusLensOut, nRhoLensOut, nkt, dataFileName, overwriteSavedData, rDips, dips, dims, e0);

    k0 = 2*pi/lambda0;

    % create lens positions
    xLensOut = linspace(-radiusLensOut, radiusLensOut, nRhoLensOut);
    yLensOut = xLensOut;
    dX = abs(xLensOut(2) - xLensOut(1));
    dY = abs(yLensOut(2) - yLensOut(1));
    dA = dX*dY;
    [YLensOut, XLensOut] = meshgrid(yLensOut, xLensOut);
    RhoLensOut = sqrt(XLensOut.^2 + YLensOut.^2);
    withinLens = find(RhoLensOut <= 100*radiusLensOut);
    rLensOut(:,1) = XLensOut(withinLens);
    rLensOut(:,2) = YLensOut(withinLens);
    rLensOut(:,3) = zLensOut;       

    % create laser field across waveguide/lens
    typeOfEField = ['EFieldAtOutput_z=', num2str(zLensOut./lambda0),...
                    '_R=', num2str(radiusLensOut/lambda0), '_nR=', num2str(nRhoLensOut)];
    [EL, kappa] = Func_EFieldLaser(rLensOut, EFieldModel, polVec, 'xyz', propagationDirection, c, lambda0, w0, EF, kappa, zLensIn, radiusLensIn, nRhoLensIn, nkt, dataFileName, overwriteSavedData, typeOfEField);
    
    % calculate scattered field across waveguide/lens    
    dipsOut = zeros([size(dips,1), 3, size(dips,3), size(dips,4)]);
    for iDim = 1:length(dims)
        dipsOut(:,dims(iDim),:,:) = dips(:,iDim,:,:); % [nDips, 3, nDetuning]
    end   
    if strcmp(polBasis, '+-z') % convert to xyz basis
        dipsxyz = zeros(size(dipsOut));
        dipsxyz(:, 1, :, :) = (dipsOut(:, 1, :, :) + dipsOut(:, 2, :, :))./sqrt(2);
        dipsxyz(:, 2, :, :) = 1j.*(dipsOut(:, 1, :, :) - dipsOut(:, 2, :, :))./sqrt(2);
        dipsxyz(:, 3, :, :) = dipsOut(:, 3, :, :);
        dipsOut = dipsxyz;
    end
                
    nLens = size(rLensOut, 1);
    nDips = size(rDips, 1);
    nDim1 = size(dips, 3);
    nDim2 = size(dips, 4);
    % [nDips, nLens, 3, nDim1, nDim2]
    dipsRep = repmat(permute(dipsOut, [1, 5, 2, 3, 4]), [1, nLens, 1, 1, 1]);
    rDipsRep = repmat(permute(rDips, [1, 3, 2]), [1, nLens, 1, nDim1, nDim2]);
    rLensRep = repmat(permute(rLensOut, [3, 1, 2]), [nDips, 1, 1, nDim1, nDim2]);
    rLensMag = repmat(sqrt(sum(abs(rLensRep).^2, 3)), [1,1,3,1,1]);
    rLensUnit = rLensRep./rLensMag;    
    
    rSepVec = rLensRep - rDipsRep; % vector from dipoles to positions on lens
    rSepMag = sqrt(sum(abs(rSepVec).^2, 3)); % separation magnitude
    XVec = rSepVec(:, :, 1, :, :, :, :, :);
    YVec = rSepVec(:, :, 2, :, :, :, :, :);
    ZVec = rSepVec(:, :, 3, :, :, :, :, :);
    XSep = XVec./rSepMag; % x y and z components of separations
    YSep = YVec./rSepMag;
    ZSep = ZVec./rSepMag;
    XVec = repmat(XVec, [1, 1, 3]);
    YVec = repmat(YVec, [1, 1, 3]);
    ZVec = repmat(ZVec, [1, 1, 3]);
    
    Goff = 3./(k0.*rSepMag).^3 - 3j./(k0.*rSepMag).^2 - 1./(k0.*rSepMag); % p ~= q
    Gdir = -1./(k0.*rSepMag).^3 + 1j./(k0.*rSepMag).^2 + 1./(k0.*rSepMag); % p = q
    % Goff =  - 1./(k0.*rsep); disp('Long range field'); % p ~= q
    % Gdir = 1./(k0.*rsep); % p = q

    % field from interacting dipoles
    Ed = zeros(size(rLensRep));  
    Ed(:,:,1,:,:,:,:,:) = (XSep.*XSep.*Goff + Gdir).*dipsRep(:,:,1,:,:,:,:,:) + ...
                           XSep.*YSep.*Goff.*dipsRep(:,:,2,:,:,:,:,:) + ...
                           XSep.*ZSep.*Goff.*dipsRep(:,:,3,:,:,:,:,:);
    Ed(:,:,2,:,:,:,:,:) = (YSep.*YSep.*Goff + Gdir).*dipsRep(:,:,2,:,:,:,:,:) + ...
                           YSep.*XSep.*Goff.*dipsRep(:,:,1,:,:,:,:,:) + ...
                           YSep.*ZSep.*Goff.*dipsRep(:,:,3,:,:,:,:,:);
    Ed(:,:,3,:,:,:,:,:) = (ZSep.*ZSep.*Goff + Gdir).*dipsRep(:,:,3,:,:,:,:,:) + ...
                           ZSep.*XSep.*Goff.*dipsRep(:,:,1,:,:,:,:,:) + ...
                           ZSep.*YSep.*Goff.*dipsRep(:,:,2,:,:,:,:,:);
    Ed = Ed .* k0.^3./e0./(4*pi).*exp(1j.*k0.*repmat(rSepMag,[1,1,3]));
    EdTot = permute(sum(Ed,1),[2,3,4,5,6,7,8,1]); % sum over dipoles       
    
    % calculate signal
    ELRep = repmat(EL, [1, 1, nDim1, nDim2]);
    ELSteadyState = ELRep;
    kLUnit = permute(rLensUnit(1,:,:,:,:), [2,3,4,5,1]); % laser wavevector, [nLens, 3, dim1, dim2]
     
    RWhere = 1; 
    dA = repmat([0,0,1].*dA, [nLens, 1, nDim1, nDim2]);
    couplingAmplitudeFunc = @(g, E, kg) e0*c/2.*abs(sum(dot(g, E, 2) .* dot(conj(kg), dA, 2) .* RWhere, 1)).^2;                                        
    PInAna = 1/4*e0*pi*c*EF^2*w0^2; % analytic input beam power
    g = ELSteadyState./sqrt(PInAna);   
    PL = couplingAmplitudeFunc(g, ELRep, kLUnit);
    PdTot = couplingAmplitudeFunc(g, EdTot, kLUnit);
    PTot = couplingAmplitudeFunc(g, ELRep + EdTot, kLUnit); 
    
end