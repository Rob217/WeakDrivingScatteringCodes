% add random kick to position of atoms
function [rDipsOut] = Func_PositionFluctuations(rDipsIn, positionFluctuations)

if strcmp(positionFluctuations.type, 'stdDev')
    stdDev = positionFluctuations.parameters;
    
elseif strcmp(positionFluctuations.type, 'trapDepth')
    % trap depth: V = s*Erec 
    % Wernier wavefucntion: |psi|^2 = exp(-rho^2/l^2)
    % rms width: l = a/(pi*s^0.25) 
    % kick=1/(pi*s^0.25) 
    % l = a*kick  OR  l = s
    % e.g. see Jenkins PRA 86 031602 (2012)
    % typical values for s = 35 -> 5500 (Bakr Nat 462 74 (2009))
    % -> kick = 0.13, 0.037

    % MODIFIED...
    % the relation between kick and V0 is:
    % kick = (1/pi)*1/(4*s)^0.25
    % or
    % (V0/Er) = 1/(4*pi^4*kick^4)
    
    trapDepth = positionFluctuations.parameters;
    a = 1; % lattice spacing - since we are using this for rDips0

    kick = (1/pi)*1./(4.*trapDepth).^0.25;
    kick(find(isinf(kick))) = 0; % account for infs
    stdDev = a.*kick;
    
%     if length(kick) == 1 % same kick in all directions
%         l = a*kick;
%     else % different kick in x,y,z
%         l = zeros(size(a));
%         l(:,1,:,:) = a(:,1,:,:)*kick(1); 
%         l(:,2,:,:) = a(:,2,:,:)*kick(2);
%         l(:,3,:,:) = a(:,3,:,:)*kick(3);
%     end
end


rDipsOut(:,1) = normrnd(rDipsIn(:,1), stdDev(1));
rDipsOut(:,2) = normrnd(rDipsIn(:,2), stdDev(2));
rDipsOut(:,3) = normrnd(rDipsIn(:,3), stdDev(3));

end