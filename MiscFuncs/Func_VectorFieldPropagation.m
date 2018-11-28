% Calculate full vector propagation of an electric field 
% Uses a basis of angular momentum cylindrically symmetric field modes
%
% Theory based on Van Enk & Kimble PRA 63 023809 (2001) and Tey et al. NJP 11 043011 (2009)
% Steps:
% (1) Decompose given electric field in some plane z into the cylindrical modes
% (2) The cylindrical modes satisfy Maxwell's equations so are valid at all r.
%     Therefore can recompose them into a total field at any plane z
%
% Author: Rob Bettles
% Date: 1.10.18

function [EOut, kappa, rho_kt_params] = VectorFieldPropagation(rOut, zLens, polVec, polBasis, k, w0, kappa, EAtom, rho_kt_params); 
    % Inputs:
    %   rOut - positions to calculate output field
    %   zIn - z value of input plane
    %   polVec - polarization vector
    %   polBasis - the basis of polarization (either +-z or xyz)
    %   k - light wavenumber (=2*pi/lambda)
    %   w0 - beam width at z=0 (or at least what it would be assuming paraxial propagation)
    %   kappa - expansion coefficients at input plane in +-z basis ( = [] if not yet calculated)
    %   EAtom - E-field amplitude at atom location
    %   rho_kt_params [optional parameter] - parameters used for rho and kt (limits, number of points)
    %
    % Outputs:
    %   EOut - output electric field over rOut  
    %   kappa - calculated expansion coefficients
    %
    % Assumptions in this code:
    %   1.  The input field polarization components have the following dependences on phi:
    %       F+ ~ 1
    %       Fz ~ exp(1j*phi)
    %       F- ~ exp(2j*phi)
    %       This is important because it means that they cancel the exp(1j*m*phi) term in GFunc
    %       Otherwise, would need to consider additionally the integral over phi
    %   2.  The input polarization must be |+> = (|x> + 1j*|y>)/sqrt(2)
    %       To consider different polarizations, we would need to include explicitly exp(1j*m*phi) in the GFunc
    %       and this would not necessarily cancel with the phase terms in the polarization of the driving field (therefore we would again need to consider the phi dependence explicitly)
    
    newMethod = 1; % new method = 1, old method = 0;

    % allow us to run script stand-alone
    if nargin == 0
        lambda = 780e-9;
%         lambda = 1;
        
%         rhoMax = (30e-6)./(780e-9).*lambda;
%         rhoMax = (10e-6)./(780e-9).*lambda;
%         rhoMax = (2.5e-6)./(780e-9).*lambda;
        rhoMax = (1e-6)./(780e-9).*lambda;
        rhoOut = linspace(0, rhoMax, 50);
        zOut = 1.2.*lambda;
        rOut(:, 1) = rhoOut;
        rOut(:, 2) = 0;
        rOut(:, 3) = zOut;        
        
        EInModel = 1;
        polVec = [1, 1j, 0]./sqrt(2);
        polBasis = '+-z';
%         polBasis = 'xyz';
        angularMomentum = 1;
        k = 2*pi/lambda;
        
        fLens = 100.5*lambda;
%         wL = (0.1e-3)./(780e-9).*lambda;
%         wL = (0.3e-3)./(780e-9).*lambda;
%         wL = (1.1e-3)./(780e-9).*lambda;
%         fLens = (4.5e-3)./(780e-9).*lambda;
%         w0 = 1.2*lambda;
        u = 2.22;
        wL = u*fLens;
        w0 = sqrt((wL^2 - sqrt(wL^4 - 4*fLens^2*lambda^2/pi^2))/2);
        
        zLens = -fLens; 
        EAtom = 1;
        
        kappa = [];         
        plotFigs = 1;
    else
        plotFigs = 0;
        fLens = abs(zLens);
        
    end
        
    if exist('rho_kt_params') == 1
        nkt = rho_kt_params.nkt;
        nRhoLens = rho_kt_params.nRhoLens;
        rhoLensMax = rho_kt_params.rhoLensMax;
    else
        nkt = 500;
        nRhoLens = 500;
        rhoLensMax = 2; % multiples of wL
    end            
        
    % determine angular momentum (the only reason we consider additional polarization possibilities is for if we want to later change this code to account for other polarizations)
    if all(polVec == [1, 1j, 0]./sqrt(2))
        angularMomentum = 1;
    elseif all(polVec == [1, -1j, 0]./sqrt(2))
        angularMomentum = -1;
    elseif all(polVec==[1,0,0]) || all(polVec==[0,1,0])
        angularMomentum = 0;
    end
        
    % check whether polarization is |+>
    if abs(dot([1, 1j, 0]./sqrt(2), polVec) - 1) > 1e-8 
        error('### Polarization must be |+> = (|x> + 1j|y>)/sqrt(2) for vector focussing code to work ###')
    end
   
    
    % calculate beam radius at lens (assuming paraxial propagation - this is just to give us an estimate of the input beam size)
    lambda = 2*pi/k;
    zR = pi.*w0.^2./lambda; % Rayleigh range
    wz = w0.*sqrt(1+(zLens./zR).^2); % beam radius at z
    wL = wz;

    % calculate amplitude of field at lens
    ELens = EAtom.*w0./wz;

    % kt = transverse component of wavevector, kz = longitudinal component of wavevector, k = sqrt(kz^2+kt^2) = magnitude of wavevector
    ktMax = k; % upper integral limit for kt points
%     nkt = 2001; % number of points in kt integral
    kt = linspace(0, ktMax, nkt); % kt points
    dkt = kt(2) - kt(1);
    
    % mode function
    GFunc = @(m, kt, kz, rho, z, phi) besselj(m, kt.*rho).*exp(1j.*kz.*z).*exp(1j.*m.*phi);

    % calculate overlap of EIn with field modes, kappa
    if length(kappa) == 0

        u = wL/fLens;
%         rhoLensMax = 2.5*wz; % limit of rho integral (should be large enough that integral converges)
        rhoLensMax = rhoLensMax.*wz;
%         nRhoLens = 2001; % number of points in rho intengral
        rhoLens = linspace(0, rhoLensMax, nRhoLens); % rho points
        dRhoLens = rhoLens(2) - rhoLens(1); % difference in rho points

        % add dimensions [nkt, nrho0]
        rhoLensRep = repmat(rhoLens, [nkt, 1]); % rho at lens
        ktRep = repmat(permute(kt, [2, 1]), [1, nRhoLens]);

        kzRep = sqrt(k.^2 - ktRep.^2); % longitudinal component

        % calculate bessel function terms
        % GFunc = @(m, kt, kz, rho, z, phi)
        if newMethod == 1
            GPlus = GFunc(angularMomentum - 1, ktRep, kzRep, rhoLensRep, zLens, 0); % NEW METHOD
            GMinus = GFunc(angularMomentum + 1, ktRep, kzRep, rhoLensRep, zLens, 0);
            GZ = GFunc(angularMomentum, ktRep, kzRep, rhoLensRep, zLens, 0);        
        else
            GPlus = GFunc(angularMomentum - 1, ktRep, kzRep, rhoLensRep, 0, 0); % OLD METHOD
            GMinus = GFunc(angularMomentum + 1, ktRep, kzRep, rhoLensRep, 0, 0);
            GZ = GFunc(angularMomentum, ktRep, kzRep, rhoLensRep, 0, 0);
        end
        
        % create EField over input lens - at the moment this just does a naiive lens focusing
        simpleFocussing = 0;
        if simpleFocussing == 1
            disp('### change this part to include proper lens focusing ###')
            fLens = abs(zLens); % lens focal length
            lensPhase = exp(-1j.*k.*rhoLensRep.^2./2./fLens);
            FInPlus = exp(-rhoLensRep.^2./wz.^2).*lensPhase;
            FInMinus = zeros(size(FInPlus));
            FInZ = zeros(size(FInPlus));
        else
            thetaLens = atan2(rhoLensRep,-zLens);
            phiLens = 0;
%             beamProfile = exp(-rhoLensRep.^2./wL.^2).*exp(1j.*sign(zLens).*(k.*sqrt(rhoLensRep.^2+fLens.^2)-pi/2));
            beamProfile = exp(-rhoLensRep.^2./wL.^2).*exp(1j.*sign(zLens).*(k.*sqrt(rhoLensRep.^2+fLens.^2)));
            FInPlus = (1-sign(zLens).*cos(thetaLens))./2./sqrt(cos(thetaLens)).*beamProfile;
            FInMinus = (-sign(zLens).*cos(thetaLens)-1)./2.*exp(1j.*2.*phiLens)./sqrt(cos(thetaLens)).*beamProfile;
            FInZ = -sign(zLens).*sin(thetaLens).*exp(1j.*phiLens)./sqrt(2)./sqrt(cos(thetaLens)).*beamProfile;
        end

        % kappa = [nkt, s=+-1, +-z, :, :, :, :, ...]        
        F0_polP_spinP = 2*pi/4/pi.*(+1*k + kzRep)./k.*GPlus; % pol=+, s=+1 (the 2*pi takes into account the integral over phi)
        F0_polP_spinM = 2*pi/4/pi.*(-1*k + kzRep)./k.*GPlus; % pol=+, s=-1
        F0_polM_spinP = 2*pi/4/pi.*(+1*k - kzRep)./k.*GMinus; % pol=-, s=+1
        F0_polM_spinM = 2*pi/4/pi.*(-1*k - kzRep)./k.*GMinus; % pol=-, s=-1
        F0_polZ_spinP = -1j.*2*pi.*sqrt(2)./4./pi.*ktRep./k.*GZ; % pol=z, s=+1
        F0_polZ_spinM = -1j.*2*pi.*sqrt(2)./4./pi.*ktRep./k.*GZ; % pol=z, s=-1
        kappa = zeros([nkt, 2, 3]);
        % to speed up: could bring the ktRep out of the sum (in fact probably don't need it in the sum in the first place)
        kappa(:,1,1) = sum(2*pi*ktRep.*rhoLensRep.*FInPlus.*conj(F0_polP_spinP), 2);
        kappa(:,2,1) = sum(2*pi*ktRep.*rhoLensRep.*FInPlus.*conj(F0_polP_spinM), 2);
        kappa(:,1,2) = sum(2*pi*ktRep.*rhoLensRep.*FInMinus.*conj(F0_polM_spinP), 2);
        kappa(:,2,2) = sum(2*pi*ktRep.*rhoLensRep.*FInMinus.*conj(F0_polM_spinM), 2);
        kappa(:,1,3) = sum(2*pi*ktRep.*rhoLensRep.*FInZ.*conj(F0_polZ_spinP), 2);
        kappa(:,2,3) = sum(2*pi*ktRep.*rhoLensRep.*FInZ.*conj(F0_polZ_spinM), 2);        
        kappa = sum(kappa,3).*dRhoLens;        
        
    end
   
    % plot convergence of rho integral
    plotRhoIntegralConvergence = 1;
    if plotRhoIntegralConvergence == 1 && plotFigs == 1
        figure(678)
        clf
        hold on, box on
        kappaCum(:,:,1,1) = cumsum(2*pi*ktRep.*rhoLensRep.*FInPlus.*conj(F0_polP_spinP), 2);
        kappaCum(:,:,2,1) = cumsum(2*pi*ktRep.*rhoLensRep.*FInPlus.*conj(F0_polP_spinM), 2);
        kappaCum(:,:,1,2) = cumsum(2*pi*ktRep.*rhoLensRep.*FInMinus.*conj(F0_polM_spinP), 2);
        kappaCum(:,:,2,2) = cumsum(2*pi*ktRep.*rhoLensRep.*FInMinus.*conj(F0_polM_spinM), 2);
        kappaCum(:,:,1,3) = cumsum(2*pi*ktRep.*rhoLensRep.*FInZ.*conj(F0_polZ_spinP), 2);
        kappaCum(:,:,2,3) = cumsum(2*pi*ktRep.*rhoLensRep.*FInZ.*conj(F0_polZ_spinM), 2);  
        
        ktPlot = [100:100:nkt];
        for ikt = 1:length(ktPlot)            
            for iSpin = 1:2
                for iPol = 1:3
                    plot(rhoLens, abs(kappaCum(ikt,:,iSpin,iPol)))
                end
            end
        end
    end
    
    % add dimensions of positions calculating field
    nOut = size(rOut, 1); % number of locations in rOut
    rhoOut = repmat(sqrt(rOut(:,1).^2 + rOut(:,2).^2), [1, nkt]);
    phiOut = repmat(atan2(rOut(:,2), rOut(:,1)), [1, nkt]);
    zOut = repmat(rOut(:,3), [1, nkt]);
        
    ktRep2 = repmat(kt, [nOut, 1]);
    kzRep2 = sqrt(k.^2 - ktRep2.^2);
    
    if newMethod == 1
        GPlus2 = GFunc(angularMomentum-1, ktRep2, kzRep2, rhoOut, zOut, phiOut);
        GMinus2 = GFunc(angularMomentum+1, ktRep2, kzRep2, rhoOut, zOut, phiOut); 
        GZ2 = GFunc(angularMomentum, ktRep2, kzRep2, rhoOut, zOut, phiOut);
    else
        GPlus2 = GFunc(angularMomentum-1, ktRep2, kzRep2, rhoOut, zOut + fLens, phiOut); % OLD METHOD
        GMinus2 = GFunc(angularMomentum+1, ktRep2, kzRep2, rhoOut, zOut + fLens, phiOut); 
        GZ2 = GFunc(angularMomentum, ktRep2, kzRep2, rhoOut, zOut + fLens, phiOut);
    end
    kappaRep = repmat(permute(kappa, [4, 1, 2, 3]), [nOut, 1, 1, 1]);
        
    EOut(:, 1) = sum(dkt./4./pi.*(+k+kzRep2)./k.*GPlus2.*kappaRep(:, :, 1), 2) + ...
                 sum(dkt./4./pi.*(-k+kzRep2)./k.*GPlus2.*kappaRep(:, :, 2), 2);
    EOut(:, 2) = sum(dkt./4./pi.*(+k-kzRep2)./k.*GMinus2.*kappaRep(:, :, 1), 2) + ...
                 sum(dkt./4./pi.*(-k-kzRep2)./k.*GMinus2.*kappaRep(:, :, 2), 2);
    EOut(:, 3) = sum(dkt./4./pi.*(-1j*sqrt(2)).*ktRep2./k.*GZ2.*kappaRep(:, :, 1), 2) + ...
                 sum(dkt./4./pi.*(-1j*sqrt(2)).*ktRep2./k.*GZ2.*kappaRep(:, :, 2), 2);                

    if newMethod == 0 % old method
        EOut = EOut.*1j;
    end        
          
    % plot field in focal plane
    plotFieldInFocalPlane = 1;
    if plotFieldInFocalPlane == 1 && plotFigs == 1
        figure(256)
        clf
        hold on, box on
        
        rhoOut = sqrt(rOut(:,1).^2 + rOut(:,2).^2);
%         scaleFactor = lambda;
        scaleFactor = 1e-6;
        plot(rhoOut./scaleFactor, abs(EOut(:,1)), 'r')
        plot(rhoOut./scaleFactor, abs(EOut(:,2)), 'g')
        plot(rhoOut./scaleFactor, abs(EOut(:,3)), 'b')
        plot(rhoOut./scaleFactor, wL/w0.*exp(-(rhoOut./w0).^2), 'k--')
    end    
    
    % convert EIn from xyz to +-z
    if strcmp(polBasis, 'xyz')
        EOutxyz = zeros(size(EOut));
        EOutxyz(:, 1) = (EOut(:, 1) + EOut(:, 2))./sqrt(2);
        EOutxyz(:, 2) = 1j.*(EOut(:, 1) - EOut(:, 2))./sqrt(2);
        EOutxyz(:, 3) = EOut(:, 3);
        EOut = EOutxyz;
    end  
    
    EOut = EOut.*ELens;
        
    if nargin == 0
        EOut = [];
    end

end