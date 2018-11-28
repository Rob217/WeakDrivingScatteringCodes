% Code for calculating the scattering from a cloud of cold atoms
%
% Author: Robert Bettles
% Date: September 2018

    
%% CHOOSE INITIAL PARAMETERS %%  

addpath('./InputVariables/')
addpath('./PlottingFuncs/')
addpath('./MiscFuncs/')

% load variables 
if exist('runFromInputVariablesScript')==0
    clear all
    runFromInputVariablesScript = 0;
    Script_InputVariables_1
end

clear 'runFromInputVariablesScript'
cols = Func_CustomCols(); % load custom colors

tic;
t1 = toc;


%% PERFORM CALCULATIONS %%

disp('-- running calculations')
for iRep = 1:nReps % loop over repetitions
    ithRep = repsList(iRep);

    % use repetition index as seed for random number generator - this way the random numbers generated for any given repetition index will always be the same
    rng(ithRep) 

    % recalculate atom positions
    [nDips, rDips0] = Func_LatticePositions(nX, nY, nZ, shape);
    nModes = nDims*nDips;
    
    rDipsWithoutFluctuations = rDips0;
    [rDips0] = Func_PositionFluctuations(rDips0, positionFluctuations);

    % create Zeeman shift matrix
    if iRep == 1
        ZeemanShiftMat = zeros([3,3]);
        if strcmp(polBasis, '+-z')
            ZeemanShiftMat(1,1) = muB;
            ZeemanShiftMat(2,2) = -muB;
        elseif strcmp(polBasis, 'xyz')
            ZeemanShiftMat(1,2) = -1j*muB;
            ZeemanShiftMat(2,1) = muB;
        end
        ZeemanShiftMat = ZeemanShiftMat(dims, dims);
        ZeemanShiftMat = kron(eye(nDips), ZeemanShiftMat);
        ZeemanShiftMat = ZeemanShiftMat./alpha0./gamma0; % convert to these units
    end

    for iSep = 1:nSeps % loop over nearest-neighbor separations
        rSep = rSepRange(iSep);           

        % scale atom positions
        rDips = rDips0.*rSep;   

        % create filename
        EFieldInfo = 0; % for calculating the dipoles we don't need to know information about collection the lens properties
        fileName = Func_CreateFilename(calcSteadyState,calcDynamics,nX,nY,nZ,nDips,shape,detuningRange,rSep,ithRep,polVec,polBasis,dims,w0,EF,tSteps,tPulseEdge,tOff,edgeShape,tau0,Gamma0,lambda0,EFieldModel,radiusLensIn,fLensIn,zLensIn,nRhoLensIn,nkt,EFieldInfo,muB,justReIm,positionFluctuations);
        dataFileName = ['./Data/',fileName];
        if exist([dataFileName])==0
            mkdir([dataFileName])
        end      

        % calculate electric field at atom positions
        propagationDirection = [0, 0, 1]; % direction of beam propagation                        

        typeOfEField = 'EFieldAtDipoles';
        [EL, kappa] = Func_EFieldLaser(rDips, EFieldModel, polVec, polBasis, propagationDirection, c, lambda0, w0, EF, kappa, zLensIn, radiusLensIn, nRhoLensIn, nkt, dataFileName, overwriteSavedData, typeOfEField);            

        for iDim = 1:3
            EFVec(iDim:3:nDips*3,1) = EL(:,iDim);
        end           
        redInds = []; % reduced dimension indices, e.g. [1,2,4,5,7,8,...]
        for iDim = 1:nDims
            redInds = sort([redInds,dims(iDim):3:3*nDips]);
        end
        EFVec = EFVec(redInds);          

        clear 'dipSS_allDetuning'
        for iDetuning = 1:nDetunings % loop over detunings                
            % print progress
            progressText = [];
            if nReps>1
                progressText = [progressText, 'iRep=',num2str(iRep),'/',num2str(nReps),'  '];
            end
            if nSeps > 1
                progressText = [progressText, 'iSep=',num2str(iSep),'/',num2str(nSeps),'  '];
            end
            if nDetunings > 1
                progressText = [progressText, 'iDetuning=',num2str(iDetuning),'/',num2str(nDetunings),'  '];
            end
            t2 = toc;
            dTime = t2 - t1;
            progressText = [progressText, 'tStep=', num2str(dTime), '  '];
            if displayProgress == 1 && dTime>1       
                disp(progressText)
                t1 = t2;
            end
            Delta = detuningRange(iDetuning);                

            calcDips = 0; % set to 0 - this is changed if we do actually calculate dipoles
            
            %% calculate steady state using matrix inversion method %%
            if calcSteadyState == 1                         

                % calculate/load dipole data
                if exist([dataFileName,'/dipSS.mat'])==0 || overwriteSavedData==1 % check whether already exists
                    calcDips = 1;

                    % calculate interaction matrix
                    if iDetuning==1
                        Gij = Func_CalculateInteractionMatrixGij(rDips, k0, dims, e0, polBasis, justReIm);
                    end

                    % diagonal elements in coupling matrix M - driving and singl-atom decay and frequencies
                    alphaInv = -((Delta+1j.*gamma0).*eye(nDims*nDips))./alpha0./gamma0  + ZeemanShiftMat;
                    if strcmp(justReIm, 'Re')
                        disp('Just taking real part of alphaInv')
                        alphaInv = real(alphaInv);
                    elseif strcmp(justReIm, 'Im')
                        disp('Just taking imag part of alphaInv')
                        alphaInv = imag(alphaInv);
                    end

                    % coupling matrix - diagonal elements include driving and bare atom frequency; off-diagonal elements include interactions
                    M = alphaInv - Gij;   

                    % calculate dipole moments 
                    dipSol = M\EFVec; % dipole moments, dimension is [d1x, d1y, d1z, d2x, d2y, d2z, ...] where, e.g., d1x is the x component of atom 1
                    clear 'dip'
                    for iDim = 1:nDims
                        dip(:,iDim) = dipSol(iDim:nDims:end); % convert to [nDips, nDims]
                    end
                    dipSS_allDetuning(:,:,iDetuning) = dip; % save steady state dipoles

                    % save dipole data
                    if iDetuning==nDetunings
                        save([dataFileName,'/dipSS.mat'], 'dipSS_allDetuning')
                    end

                % load dipole data
                elseif iDetuning==1                        
                    load([dataFileName,'/dipSS.mat'])                        
                end
                if iDetuning == nDetunings
                    dipSS(:,:,:,iRep,iSep) = dipSS_allDetuning;
                end

                % calculate/load extinction cross-section
                if calcCrossSection == 1
                    if exist([dataFileName,'/crossSectionSS.mat'])==0 || overwriteSavedData==1

                        % calculate cross-section
                        dip = dipSS_allDetuning(:,:,iDetuning);
                        for iDim = 1:nDims
                            dipSol(iDim:nDims:nDims*nDips,1) = dip(:,iDim);
                        end
                        crossSectionSS_allDetuning(iDetuning) = imag(EFVec'*dipSol); % steady state scattering cross-section

                        % save cross-section data
                        if iDetuning==nDetunings
                            save([dataFileName,'/crossSectionSS.mat'], 'crossSectionSS_allDetuning')                                
                        end       

                    % load cross-section data
                    elseif iDetuning==1 
                        load([dataFileName,'/crossSectionSS.mat'])
                    end
                    if iDetuning == nDetunings
                        crossSectionSS(:,iRep,iSep) = crossSectionSS_allDetuning;
                    end
                end

                % calculate eigenvalues and eigenvectors - only calculate for a single detuning since same for all detunings     
                if calcEigs == 1 && iDetuning == 1

                    if exist([dataFileName,'/eigShifts.mat'])==0 || overwriteSavedData==1

                        % check whether coupling matrix has already been defined
                        if calcDips == 0;
                            % calculate interaction matrix
                            if iSep==1 && iDetuning==1
                                Gij = Func_CalculateInteractionMatrixGij(rDips, k0, dims, e0, polBasis, justReIm);
                            end

                            % diagonal elements in coupling matrix M - driving and singl-atom decay and frequencies
                            alphaInv = -((Delta+1j.*gamma0).*eye(nDims*nDips))./alpha0./gamma0  + ZeemanShiftMat;

                            if strcmp(justReIm, 'Re')
                                disp('Just taking real part of alphaInv')
                                alphaInv = real(alphaInv);
                            elseif strcmp(justReIm, 'Im')
                                disp('Just taking imag part of alphaInv')
                                alphaInv = imag(alphaInv);
                            end

                            % coupling matrix - diagonal elements include driving and bare atom frequency; off-diagonal elements include interactions
                            M = alphaInv - Gij;
                        end                            

                        if strcmp(justReIm, 'Re')
                            disp('Just taking real part of alphaInv')
                            alphaInv = real(alphaInv);
                        elseif strcmp(justReIm, 'Im')
                            disp('Just taking imag part of alphaInv')
                            alphaInv = imag(alphaInv);
                        end

                        M = alphaInv - Gij;                            

                        % interaction matrix Gij is complex symmetric (NOT hermitian), so eigenvectors are nonorthogonal and eigenvalues are complex                            
%                             M = M.*alpha0./2; % rescale to make these match with other code
                        [eigVecsSingleRun,eigValsSingleRun] = eig(M);

%                             eigValsSingleRun = eigValsSingleRun./alpha0.*2; % rescale
                        eigValsSingleRun = diag(eigValsSingleRun); % convert from diagonal matrix to list
                        eigShiftsSingleRun = alpha0.*gamma0.*real(eigValsSingleRun) + detuningRange(iDetuning); % eigenmode energy shifts
                        eigGammasSingleRun = -alpha0.*Gamma0.*imag(eigValsSingleRun); % eigenmode decay rates
                        eigVecOverlap = transpose(eigVecsSingleRun)*eigVecsSingleRun; % overlap of eigenvectors with each other (left eigenvector for complex symmetric matrix = transpose of right eigenvector)
                        bCoefficient = transpose(eigVecsSingleRun)*EFVec./diag(eigVecOverlap); % overlap of eigenvectors with driving field - EFVec = sum_l b_l eigVec_l, i.e. b_l is the expansion coefficient 
%                         [dummy,sortInds] = sort(eigGammasSingleRun); % sort modes in ascending order of decay rates
                        [dummy,sortInds] = sort(bCoefficient, 'ascend'); % sort modes in ascending order of mode overlap with driving field
                        eigShifts_allDetuning(:,iDetuning) = eigShiftsSingleRun(sortInds); % sort shifts 
                        eigGammas_allDetuning(:,iDetuning) = eigGammasSingleRun(sortInds); % sort decay rates
                        eigVecs_allDetuning(:,:,iDetuning) = eigVecsSingleRun(:,sortInds);
                        modePops_allDetuning(:,iDetuning) = abs(bCoefficient(sortInds)).^2./EF^2; % relative overlap of each mode with driving field

                        % save cross-section data
                        save([dataFileName,'/eigShifts.mat'], 'eigShifts_allDetuning')
                        save([dataFileName,'/eigGammas.mat'], 'eigGammas_allDetuning')
                        save([dataFileName,'/eigVecs.mat'], 'eigVecs_allDetuning')
                        save([dataFileName,'/modePops.mat'], 'modePops_allDetuning', 'bCoefficient')

                    % load cross-section data
                    elseif iDetuning==1 
                        load([dataFileName,'/eigShifts.mat'])
                        load([dataFileName,'/eigGammas.mat'])
                        load([dataFileName,'/eigVecs.mat'])
                        load([dataFileName,'/modePops.mat'])
                    end       
%                         if iDetuning == nDetuning
                    eigShifts(:,:,iRep,iSep) = eigShifts_allDetuning;
                    eigGammas(:,:,iRep,iSep) = eigGammas_allDetuning;
                    eigVecs(:,:,:,iRep,iSep) = eigVecs_allDetuning;
                    modePops(:,:,iRep,iSep) = modePops_allDetuning;
%                         end

                    % calculate cross-section of individual modes
                    if calcCrossSection == 1 && iDetuning == 1

                        if exist([dataFileName,'/crossSectionEigsSS.mat'])==0 || overwriteSavedData==1      

                            % OLD 
%                                 bModeRep = repmat(transpose(bCoefficient), [nModes, 1]);
%                                 eigValsRep = repmat(transpose(eigValsSingleRun), [nModes, 1]);
%                                 CSWithoutEigVals = (eigVecsSingleRun'*eigVecsSingleRun).*bModeRep.*bModeRep';
%                                 eigValsRepWithDetuning = repmat(transpose(eigValsSingleRun), [nModes, 1, nDetuning]);
%                                 detuningRangeRep = repmat(permute(detuningRange, [1, 3, 2]), [nModes, nModes, 1]) - detuningRange(1);                                                                        
%                                 eigValsRepWithDetuning = eigValsRepWithDetuning - detuningRangeRep./alpha0./gamma0;
%                                 CSWithoutEigValsRep = repmat(CSWithoutEigVals, [1, 1, nDetuning]);
%                                 crossSectionEigsSS_allDetuning = imag(CSWithoutEigValsRep./eigValsRepWithDetuning);
                            bModeRep = repmat(transpose(bCoefficient), [nModes, 1]);
                            eigValsRep = repmat(transpose(eigValsSingleRun), [nModes, 1]);                                
                            eigValsRepWithDetuning = repmat(transpose(eigValsSingleRun), [nModes, 1, nDetunings]);
                            detuningRangeRep = repmat(permute(detuningRange, [1, 3, 2]), [nModes, nModes, 1]) - detuningRange(1);                                                                        
                            eigValsRepWithDetuning = eigValsRepWithDetuning - detuningRangeRep./alpha0./gamma0;
                            CSWithoutEigVals = (eigVecsSingleRun'*eigVecsSingleRun).*bModeRep.*bModeRep';
%                                 CSWithoutEigValsRep = repmat(CSWithoutEigVals, [1, 1, nDetunings]);
                            if calcEigInterferences == 1 % eigenmode cross-section interferences as well
                                CSWithoutEigValsRep = repmat(CSWithoutEigVals, [1, 1, nDetunings]);
                                crossSectionEigsSS_allDetuning = imag(CSWithoutEigValsRep./eigValsRepWithDetuning);
                            else % just single mode terms
                                CSWithoutEigValsRep = repmat(diag(CSWithoutEigVals), [1, 1, nDetunings]);
                                crossSectionEigsSS_allDetuning = imag(CSWithoutEigValsRep./permute(eigValsRepWithDetuning(1,:,:), [2,1,3])); % problem is with this line
                            end
                            save([dataFileName,'/crossSectionEigsSS.mat'], 'crossSectionEigsSS_allDetuning')                                    
                        else
                            load([dataFileName,'/crossSectionEigsSS.mat'])
                        end
                        if calcEigInterferences == 1
                            crossSectionEigsSS_inclInterferences(:,:,:,iRep,iSep) = crossSectionEigsSS_allDetuning;
                            for iMode = 1:nModes
                                crossSectionEigsSS(iMode,1,:,iRep,iSep) = crossSectionEigsSS_allDetuning(iMode,iMode,:);                                   
                            end
                        else
                            crossSectionEigsSS(:,:,:,iRep,iSep) = crossSectionEigsSS_allDetuning;
                        end
                    end

                end                                                            
            end                    


            %% calculate dynamics %%
            if calcDynamics == 1

                % calculate dipole dynamics
                if exist([dataFileName,'/dipDynamics.mat'])==0 || overwriteSavedData==1

                    if calcDips==0 && iSep==1 && iDetuning==1
                        Gij = Func_CalculateInteractionMatrixGij(rDips, k0, dims, e0, polBasis, justReIm);
                    end

                    % diagonal elements in coupling matrix M - driving and singl-atom decay and frequencies
                    alphaInv = -((Delta+1j.*gamma0).*eye(nDims*nDips))./alpha0./gamma0 + ZeemanShiftMat;

                    if strcmp(justReIm, 'Re')
                        disp('Just taking real part of alphaInv')
                        alphaInv = real(alphaInv);
                    elseif strcmp(justReIm, 'Im')
                        disp('Just taking imag part of alphaInv')
                        alphaInv = imag(alphaInv);
                    end

                    dIn = zeros([nDips*nDims, 1]); % input dipole vector
                    tol = 1e-8;
                    optionsODE45 = odeset('RelTol',tol,'AbsTol',tol); % options for ODE45 solver
                    pgeIn = dIn./D0;
                    [tOut, pgeOut] = ode45(@(t, pge) Func_WeakDrivingEOM(t, pge, gamma0, alpha0, alphaInv, D0, EFVec, Gij, tPulseEdge, tOff, edgeShape), tSteps, pgeIn, optionsODE45);
                    dipDynamics_allDetuning(:,:,iDetuning,iRep,iSep) = pgeOut.*D0;

                    % save dipolar dynamics data
                    if iDetuning == nDetunings
                        save([dataFileName,'/dipDynamics.mat'], 'dipDynamics_allDetuning')                                   
                    end

                % load dipolar dynamics data
                elseif iDetuning==1 
                    load([dataFileName,'/dipDynamics.mat'])
                end       

                % assign to object
                if iDetuning == nDetunings
                    dipDynamics(:,:,:,iRep,iSep) = dipDynamics_allDetuning;
                end


                % calculate cross-section dynamics
                if exist([dataFileName,'/crossSectionDynamics.mat'])==0 || overwriteSavedData==1

                    if iDetuning == nDetunings
                        EFVecRep = repmat(permute(EFVec, [2,1]), [nTSteps, 1, nDetunings]);
                        crossSectionDynamics_allDetuning = permute(imag(dot(EFVecRep,dipDynamics_allDetuning,2)), [1,3,2]);

                        % save dipolar dynamics data                        
                        save([dataFileName,'/crossSectionDynamics.mat'], 'crossSectionDynamics_allDetuning')                                   
                    end

                % load dipolar dynamics data
                elseif iDetuning==1 
                    load([dataFileName,'/crossSectionDynamics.mat'])
                end       

                % assign to object
                if iDetuning == nDetunings
                    crossSectionDynamics(:,:,iRep,iSep) = crossSectionDynamics_allDetuning;
                end

            end

        end  

        %% calculate scattered electric fields
        if calcEFields == 1            
            disp('-- calculate E-fields')

            % check for whether data already exists                                
            lensOutFileName = [dataFileName,'/LensOut_R=', num2str(radiusLensOut./lambda0), '_z=',num2str(zLensOut./lambda0), '_nR=', num2str(nRhoLensOut)];
            if exist(lensOutFileName)==0
                mkdir(lensOutFileName)
            end                
            waveguideTransmissionFileName = ['TransmissionWaveguide'];

            if exist([lensOutFileName,'/',waveguideTransmissionFileName,'.mat']) == 0 || overwriteSavedData == 1
                if calcSteadyState == 1
                    dip = dipSS_allDetuning;
                elseif calcDynamics == 1
                    dip = permute(dipDynamics_allDetuning, [2, 4, 3, 1]);
                    for iDim = 1:nDims
                        dip2D(:,iDim,:,:) = dip(iDim:nDims:end,:,:,:);
                    end
                    dip = dip2D;
                end
                [PL_allDetuning, PdTot_allDetuning, PTot_allDetuning] = Func_TransmissionWaveguide(EFieldModel, polVec, polBasis, propagationDirection, c, lambda0, w0, EF, kappa, zLensIn, radiusLensIn, nRhoLensIn, zLensOut, radiusLensOut, nRhoLensOut, nkt, dataFileName, overwriteSavedData, rDips, dip, dims, e0);
                save([lensOutFileName,'/',waveguideTransmissionFileName,'.mat'], 'PL_allDetuning', 'PdTot_allDetuning', 'PTot_allDetuning')
            else
                load([lensOutFileName,'/',waveguideTransmissionFileName,'.mat']);
            end
            PL(:, iRep, iSep, :) = PL_allDetuning;
            PdTot(:, iRep, iSep, :) = PdTot_allDetuning;
            PTot(:, iRep, iSep, :) = PTot_allDetuning;

        end
    end        
end


%% PLOTTING %%
disp('-- plotting figures')

iFig = 0; % figure counter
extraPlotInfo = ''; % extra information for plot filenames
EFieldPlotNames(1) = 0; % list of figures requiring EField information in figurename
textHandles(1) = 0; % list of handles of text objects that I don't want the fontsize to be changed
textSizes(1) = 10; % list of fontsizes of above text objects
if exist('plotFigs')==0
    plotFigs = 1;
end

% plot beam waist over atoms
plotBeamWaistAtomPlane = 0;
if plotBeamWaistAtomPlane == 1 && plotFigs==1
    figNum = 1101;
    Script_PlotBeamWaist;
end    

if calcSteadyState == 1 && plotFigs==1
    % plot cross-section
    plotCrossSectionLineshape = 0;
    if plotCrossSectionLineshape == 1 && calcCrossSection == 1
        figNum = 101;
        Script_PlotCrossSectionLineshape;
    end

    % plot transmission
    plotTransmissionLineshape = 0;
    if plotTransmissionLineshape == 1 && calcEFields == 1
        figNum = 253;
        Script_PlotTransmissionLineshape;
    end        

    % plot eigenvalues
    plotEigenvalueSpectrum = 0;
    if plotEigenvalueSpectrum == 1 && calcEigs == 1
        figNum = 102;
        Script_PlotEigenvalueSpectrum;
    end

    % plot eigenvector size
    plotEigenVectorSize = 0;
    if plotEigenVectorSize == 1 && calcEigs == 1
        figNum = 105;
        Script_PlotEigenvectorSize;
    end

    % plot band structure (quasimomentum vs E) and density of states
    plotBandStructure = 1;
    if plotBandStructure == 1
        Script_PlotBandStructure;
    end
    
    % plot band structure with Fourier transform of eigenmodes wrt kx and ky
    plotBandStructureFourier = 1;
    if plotBandStructureFourier == 1
        figNum = 1100;
        Script_PlotBandStructureFourier;
    end

end

if calcDynamics == 1 && plotFigs==1

    % plot cross-section dynamics
    plotCrossSectionDynamics = 0;
    if plotCrossSectionDynamics == 1
        figNum = 201;
        Script_PlotCrossSectionDynamics;
    end

    % plot transmission 
    plotTransmissionDynamics = 1;
    if plotTransmissionDynamics == 1
        figNum = 205;
        Script_PlotTransmissionDynamics;        
    end

    % plot dipoles on lattice
    plotDipoleManitudesOnLattice = 0;
    if plotDipoleManitudesOnLattice == 1
        figNum = 202;
        Script_PlotDipoleMagnitudesOnLattice;                
    end
end

% save figures
if plotFigs == 1
    Script_FigureDetails
end

clear runFromInputVariablesScript
endedScriptSuccessfully = 1;