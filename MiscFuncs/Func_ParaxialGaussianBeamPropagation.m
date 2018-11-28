%% calculate electric and magnetic field of paraxial laser beam
function [EL, BL] = ParaxialGaussianBeamPropagation(rVecs, polVec, c, lambda0, w0, EF, propagationDirection, polDim, polBasis);
    % polDim = which dimension corresponds to the polarization

    X = repmat(rVecs(:,1), [1, 3, 1, 1, 1, 1, 1]);
    Y = repmat(rVecs(:,2), [1, 3, 1, 1, 1, 1, 1]);
    Z = repmat(rVecs(:,3), [1, 3, 1, 1, 1, 1, 1]);
    k0 = 2*pi/lambda0;
    zR = pi.*w0.^2./lambda0;
    wz = w0.*(1+(Z./zR).^2).^0.5; % beam width at z
    RCurv = (Z.^2+zR.^2)./Z; % beam radius of curvature
    rho = sqrt(X.^2+Y.^2); % radius
    Gouy = atan(Z./zR); % Gouy phase
    EL = polVec.*EF.*exp(-rho.^2./wz.^2).*w0./wz.*...
         exp(1j.*k0.*Z).*...
         exp(1j.*k0.*rho.^2./2./RCurv).*...
         exp(-1j.*Gouy);      
      
    % B-field
%     BL = cross(propagationDirection, EL, polDim)./c;
    BL = 0;
    
    % rotate from xyz to +-z basis
    if strcmp(polBasis, '+-z')
        EL_pmz = zeros(size(EL));
        EL_pmz(:,1) = (EL(:,1) - 1j.*EL(:,2))./sqrt(2); 
        EL_pmz(:,2) = (EL(:,1) + 1j.*EL(:,2))./sqrt(2);
        EL_pmz(:,3) = EL(:,3);
        EL = EL_pmz;
    end
    
end