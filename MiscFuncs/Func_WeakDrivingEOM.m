%% weak driving equations of motion
function dpge_dt = WeakDrivingEOM(t,pgeIn,gamma0,alpha0,alphaInv,D0,EFVec,Gij,tPulseEdge,tOff,edgeShape);
    % this function calculates the time derivative for the coherences of each atom in an ensemble of driven interating dipoles
    
    EFPulse = Func_PulseShape(t,EFVec,tPulseEdge,tOff,edgeShape);
    
    % calculate time derivative
    dpge_dt = 1j.*(abs(alpha0)).*gamma0.*(-alphaInv*pgeIn + EFPulse./D0 + Gij*pgeIn); 

end