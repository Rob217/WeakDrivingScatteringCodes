% create filename

function [fileName] = CreateFilename(calcSteadyState,calcDynamics,nX,nY,nZ,nDips,shape,detuningRange,rSep,ithRep,polVec,polBasis,dims,w0,EF,tSteps,tPulseEdge,tOff,edgeShape,tau0,Gamma0,lambda0,EFieldModel,radiusLensIn,fLensIn,zLensIn,nRhoLensIn,nkt,EFieldInfo,muB,justReIm,positionFluctuations)

    fileName = '';
    fileName = [fileName,shape,'_N=',num2str(nX),'x',num2str(nY),'x',num2str(nZ),'=',num2str(nDips),...
                '_r=',num2str(rSep)];
            
    if any(positionFluctuations.parameters ~= 0)
        fileName = [fileName, '_posFluc=', num2str(positionFluctuations.parameters(1)), '_', num2str(positionFluctuations.parameters(2)), '_', num2str(positionFluctuations.parameters(3)), '_', positionFluctuations.type];
    end
            
    if EF~=1
        fileName = [fileName,'_EF=',num2str(EF)];
    end
    
    fileName = [fileName, fileNameVaryingDim('det', detuningRange./Gamma0)];
    
    if all(polVec == [1,0,0]); polString = ['polX'];
    elseif all(polVec == [0,1,0]); polString = ['polY'];
    elseif all(polVec == [0,0,1]); polString = ['polZ'];
    elseif all(polVec == [1,1j,0]./sqrt(2)); polString = ['pol+'];
    elseif all(polVec == [1,-1j,0]./sqrt(2)); polString = ['pol-'];
    else polString = ['pol=',numstr(polVec(1)),'_',num2str(polVec(2)),'_',num2str(polVec(3))];
    end
    fileName = [fileName, '_', polString];
    
    dimsString = '';
    for iDim = 1:length(dims)
        dimsString = [dimsString, polBasis(dims(iDim))];
    end
    fileName = [fileName, '_dims=', dimsString];
    
    if muB ~= 0
        fileName = [fileName, '_B=', num2str(muB./Gamma0)];
    end
    
    fileName = [fileName, '_EModel=', EFieldModel(1:3)];
    
    if w0./lambda0<1000 && strcmp(EFieldModel,'UniformBeam')~=1
        fileName = [fileName, '_w0=', num2str(w0./lambda0)];
    end
    
    if strcmp(EFieldModel, 'VectorGaussianBeam') == 1         
        fileName = [fileName, '_LensIn_f=',num2str(fLensIn./lambda0)];        
        if zLensIn ~= - fLensIn
            fileName = [fileName, '_z=', num2str(zLensIn./lambda0)];
        end
        fileName = [fileName, '_R=', num2str(radiusLensIn./lambda0),...
                    '_nR=', num2str(nRhoLensIn), '_nkt=', num2str(nkt)];
    end
    
    if strcmp(justReIm, 'Re') || strcmp(justReIm, 'Im')
        fileName = [fileName, '_', justReIm];
    end
    
    if calcDynamics == 1           
        fileName = [fileName, '_t=', num2str(tSteps(1)./tau0), '-', num2str(tSteps(end)./tau0),...
                    'x', num2str(length(tSteps)), '_tEdge=', num2str(tPulseEdge(1)./tau0), '_', num2str(tPulseEdge(2)./tau0)];
        if tOff < tSteps(end)
            fileName = [fileName, '_tOff=', num2str(tOff./tau0)];
        end
        fileName = [fileName,'_edge=', num2str(edgeShape)];                
    end
    
    if EFieldInfo == 1
        fileName = [fileName, '_R=',num2str(radiusLens/lambda0), '_n=', num2str(nLensDim)];
    end                        
    
    fileName = [fileName, '_rep=', num2str(ithRep)];
            
end





%% filename varying dimension
function [fileNameString] = fileNameVaryingDim(dataType,data)
    % creates a string with values of data;
    % - for section of values that are evenly spaced, string = 'data(1)_data(2)-data(1)_data(end)'
    % - otherwise string = 'data(1)_data(2)_data(3)_...'

    fileNameString = ['_',dataType,'='];

    diffData = diff(data);
    diffData = round(diffData,8,'significant'); % round to nearest 7 significant figures (this avoids floating point errors)

    if length(data) == 1 % just a single value
        fileNameString = [fileNameString,num2str(data)];
    else
        if length(unique(diffData)) == 1 % constant difference between values
            fileNameString = [fileNameString,num2str(data(1)),'_',num2str(data(2)-data(1)),'_',num2str(data(end))];
        else % varying 
            iData = 1; % data index
            while iData <= length(data)
                if iData == length(data) - 1
                    fileNameString = [fileNameString,num2str(data(iData)),',',num2str(data(iData+1))];
                    iData = iData + 1;
                elseif iData == length(data)
                    fileNameString = [fileNameString,num2str(data(iData))];
                elseif abs(diffData(iData) - diffData(iData+1)) > 1e-13
                    fileNameString = [fileNameString,num2str(data(iData))];
                else
                    fileNameString = [fileNameString,num2str(data(iData)),'_',num2str(data(iData+1)-data(iData)),'_'];
                    continueIterations = 1;
                    while continueIterations == 1
                        if iData <= (length(data) - 2)
                            if (abs(diffData(iData) - diffData(iData + 1)) < 1e-13) 
                                iData = iData + 1;
                            else
                                continueIterations = 0;
                            end
                        else
                            continueIterations = 0;
                        end
                    end
                    iData = iData + 1;
                    fileNameString = [fileNameString,num2str(data(iData))];
                end
                if iData ~= length(data)
                    fileNameString = [fileNameString,','];
                end
                iData = iData + 1;
            end
        end
    end

end