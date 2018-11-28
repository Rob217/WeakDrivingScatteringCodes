%% driving field function - pulse shape
function EFt = PulseShape(t,EF,tStep,tEnd,edgeShape)
% edgeShape = 2; % 1 = linear, 2 = sin^2, 3 = tanh

% linear increase and decrease
if edgeShape == 1 
    if t<tStep % linear increase
        EFt = EF.*t./tStep;
    elseif t<tEnd % continuously on
        EFt = EF;
    elseif t<tEnd+tStep % linear decrease
        EFt = EF.*(1-(t-tEnd)./tStep);
    else % off
        EFt = 0; 
    end
end

% sin^2 increase and decrease
if edgeShape == 2
    if t<tStep % linear increase
        EFt = EF.*sin(t./tStep.*pi/2).^2;
    elseif t<tEnd % continuously on
        EFt = EF;
    elseif t<tEnd+tStep % linear decrease
        EFt = EF.*(1-sin((t-tEnd)./tStep.*pi/2).^2);
    else % off
        EFt = 0; 
    end
end


% tanh increase and decrease fitted to data
if edgeShape == 3
    
    tStepOn = 1.97;
    if t < tStep(1) % switch on 
        t1 = 0; % or -0.53 if want pulse to switch on at t=0
        t2 = tStep(1); % 2.5, or 1.97 if want pulse to switch on at t=0
        risingEdgeWidth = 8; % tanh goes from -risingEdgeWidth/2 to + risingEdgeWidth/2
        tp = risingEdgeWidth*(t-t2)/(t2-t1) + risingEdgeWidth/2;
        EFtMax = (tanh(risingEdgeWidth/2)./2 + 0.5).*1.1 - 0.1;
%         tp = 4*(t-t2)/(t2-t1) + 2;
%         EFtMax = (tanh(2)./2 + 0.5).*1.1 - 0.1;
        EFAmplitude = ((tanh(tp)./2 + 0.5).*1.1 - 0.1)./EFtMax;
%         EFt = max(EF.*((tanh(tp)./2 + 0.5).*1.1 - 0.1)./EFtMax,0);
%         EFt = EF.*max(EFAmplitude,0);
        
    elseif t < tEnd % continuously on
%         EFt = EF;
        EFAmplitude = 1;
        
    elseif t < tEnd + tStep(2) % switch off
%     elseif t < 17.9 %%% I was getting problems here
        t1 = tEnd; % 16.8
        t2 = tEnd + tStep(2); % 18.5 = 16.8 + 1.7
        fallingEdgeWidth = 8; % tanh goes from -risingEdgeWidth/2 to + risingEdgeWidth/2
        tp = fallingEdgeWidth*(t-t2)/(t2-t1) + fallingEdgeWidth/2;
        EFtMax = (tanh(fallingEdgeWidth/2)./2 + 0.5).*1.1 - 0.1;   
%         tp = 4*(t-t2)/(t2-t1) + 2;
%         EFtMax = (tanh(2)./2 + 0.5).*1.1 - 0.1;
        EFAmplitude = ((tanh(-tp)./2 + 0.5).*1.1 - 0.1)./EFtMax;
%         EFt = max(EF.*((tanh(-tp)./2 + 0.5).*1.1 - 0.1)./EFtMax,0);
%         EFt = EF.*max(EFAmplitude,0);
        
%         disp(['t=',num2str(t),', E=',num2str(EFAmplitude)])

    else % off
%         EFt = 0;        
        EFAmplitude = 0;
    end   
    
    EFt = EF.*max(EFAmplitude,0);
    
%     if EFt < 0
%         EFt = 0;
%     end
        
end


end