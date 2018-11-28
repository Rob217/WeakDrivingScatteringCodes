function [EL, kappa] = EFieldLaser(rPos, EFieldModel, polVec, polBasis, propagationDirection, c, lambda0, w0, EF, kappa, zLensIn, radiusLensIn, nRhoLensIn, nkt, dataFileName, overwriteSavedData, typeOfEField);
    % function for calculating the laser field at positions rPos
    
    nPos = size(rPos, 1);
    k0 = 2*pi/lambda0;
    
    polVecRep = repmat(polVec, [nPos, 1]);
    if strcmp(EFieldModel, 'UniformBeam')
        propagationDirectionRep = repmat(propagationDirection, [nPos, 1]);                
        eikr = exp(1j.*dot(propagationDirectionRep, rPos, 2));
        EL = EF.*polVecRep.*repmat(eikr, [1, 3]);
        if strcmp(polBasis, '+-z')
            EL_pmz = zeros(size(EL));
            EL_pmz(:,1) = (EL(:,1) - 1j.*EL(:,2))./sqrt(2); 
            EL_pmz(:,2) = (EL(:,1) + 1j.*EL(:,2))./sqrt(2);
            EL_pmz(:,3) = EL(:,3);
            EL = EL_pmz;
        end                

    elseif strcmp(EFieldModel, 'ParaxialGaussianBeam')
        [EL, BL] = Func_ParaxialGaussianBeamPropagation(rPos, polVecRep, c, lambda0, w0, EF, propagationDirection, 2, polBasis);

    elseif strcmp(EFieldModel, 'VectorGaussianBeam')
        rho_kt_params.nkt = nkt;
        rho_kt_params.nRhoLens = nRhoLensIn;
        rho_kt_params.rhoLensMax = radiusLensIn; % multiples of wL
        EFieldFileName = [dataFileName,'/',typeOfEField,'.mat'];
        if exist(EFieldFileName)==0 || overwriteSavedData==1 % check whether already exists
            [EL, kappa] = Func_VectorFieldPropagation(rPos, zLensIn, polVec, polBasis, k0, w0, kappa, EF, rho_kt_params);                            
            save(EFieldFileName, 'EL')
        else
            load(EFieldFileName)            
        end

    elseif strcmp(EFieldModel, 'Localized')
        EL = zeros([nPos, 3]);
        [maxX, maxInd] = max(rPos(:,1));                
        EL(maxInd,:) = polVec.*EF;
        if strcmp(polBasis, '+-z')
            EL_pmz = zeros(size(EL));
            EL_pmz(:,1) = (EL(:,1) - 1j.*EL(:,2))./sqrt(2); 
            EL_pmz(:,2) = (EL(:,1) + 1j.*EL(:,2))./sqrt(2);
            EL_pmz(:,3) = EL(:,3);
            EL = EL_pmz;
        end  
    end
    
    
    
end