% if not working properly, just type 'clear all' into command window

if exist('runFromInputVariablesScript')==0 && exist('EndedScriptSuccessfully')==0
    clear all
    runFromInputVariablesScript = 1;
    PWD = pwd;
    if strcmp(PWD(end-13:end), 'InputVariables')
        cd ..
    end
else
    runFromInputVariablesScript = 0;
end

disp('-- loading input variables')

% natural constants
lambda0 = 1; % wavelength
k0 = 2*pi/lambda0; % wavenumber
c = 3e8; % speed of light
Gamma0 = 1; % spontaneous decay rate
gamma0 = Gamma0/2; % half decay rate
tau0 = 1/Gamma0; % natural lifetime
e0 = 8.85418782e-12; % permittivity of free space
alpha0 = 6*pi*e0/k0^3; % polarizability coefficient
hbar = 1.0545718e-34; % Planck constant
D0 = sqrt(6*pi*e0*hbar*gamma0/k0^3); % dipole matrix element


%% choose what to calculate
calcSteadyState = 1; % calculate steady state using matrix inversion method
calcDynamics = 0; % calculate dynamics using ODE45 solver 
calcEigs = 1; % calculate eigenmodes
calcEigInterferences = 0; % when calculating eigenmode cross-section, also save/load interference terms (take up large space)
calcEFields = 0; % calculate scattered and total electric fields   
calcCrossSection = 0; % calculate scattering cross-section
overwriteSavedData = 0; % perform calculation even if data file already exists
saveFigs = 1; % save figures (1) or not (0)

% number of repetitions
repsList = [1]; % list of repetition indices - each repetition is a seed for the random number generator, so the ith repetition will always give the same result
nReps = length(repsList); % total number of repetitions


%% atom positions
% see descriptions of lattices in Func_LatticePositions for details

% nX = 9;
% nY = 10;
% nZ = 1;
% shape = 'h';

% nX = 300;
% nY = 1;
% nZ = 1;
% dR = 3.9088;
% dZ = 0;
% shape = ['RZCyl_',num2str(dR),'_',num2str(dZ)];

% nX = 17;
% nY = 17;
% nZ = 1;
% shape = 's';

% nX = 1;
% nY = 1;
% nZ = 1;
% %     shape = 'q1_45_25_4.5_circ'; % different lattice geometries are defined in LatticePositions()
% shape = 'q1_45_5_4.5_circ'; % different lattice geometries are defined in LatticePositions()
% % %     shape = 'q1_45_20_4.8_circ'; % different lattice geometries are defined in LatticePositions()

nX = 1;
nY = 1;
nZ = 1;
shape = 'Penrose1_5';
% shape = 'Penrose2_4';

% nX = 10;
% nY = 10;
% nZ = 1;
% shape = 'Lieb';

% nearest neighbour spacing
rSepRange = 0.05;
nSeps = length(rSepRange); % number of separations considered

% atomic position fluctuations
positionFluctuations.parameters = 0.*[1,1,1];
positionFluctuations.type = 'stdDev'; % 'stdDev' (relative to lattice spacing), or 'trapDepth' in units of Er (recoil energies)



%% laser beam parameters

% laser beam parameters
w0 = 3.*lambda0; % 1/e waist
% w0Vals = [1.5:0.5:8];
% w0Vals = [1.5,2,2.5,3:10];
EF = 1; % field amplitude (because this code is weak driving only, the magnitude of EF does not matter)      
sigma0 = EF^2*alpha0; % single atom cross-section coefficient
kappa = []; % initialise mode function calculation (for VectorGaussianBeam calculation)
EFieldModel = 'UniformBeam'; 
% EFieldModel = 'ParaxialGaussianBeam';
% EFieldModel = 'VectorGaussianBeam';
% EFieldModel = 'Localized';

% detuning, Delta = omegaL - omega0, where omegaL and omega0 are the laser and atomic frequencies
if calcSteadyState == 1
%     detuningRange = 1.*linspace(-1,1,51).*Gamma0;
    detuningRange = 0;
elseif calcDynamics == 1        
    detuningRange = [-2.5:0.5:2.5].*Gamma0;
end
nDetunings = length(detuningRange);

% polarization 
% polVec = [0, 1, 0]; % polarization vector in [x,y,z]
polVec = [1,1j,0]./sqrt(2);
% polBasis = 'xyz'; % choose which basis to use
polBasis = '+-z'; 
Upmz = [[1/sqrt(2),1j/sqrt(2),0];[1/sqrt(2),-1j/sqrt(2),0];[0,0,1]]; % |+-z> = Upmz * |xyz>, rotate from xyz basis to +-z basis

% choose which dimensions to include from polBasis, i.e. if polBasis = '+-z' and dims = [1,2], then will just use + and - dimensions
dims = [1,2]; 
nDims = length(dims); % number of dimensions

% time limits and steps (for calculating dynamics)
tStart = 0; % starting time, default = 0
tEnd = 30*tau0; % end time
nTSteps = 31; % number of steps in time - if this is just 2 then ODE45 automatically chooses values for t, varying resolution as needed
tSteps = linspace(tStart, tEnd, nTSteps); % time steps (ODE45 calculates dipoles at each of these times)

% pulse profile
tPulseEdge = [0,0].*tau0; % duration of rising edge and lowering edge
tOff = 15*tau0; % time at which pulse is switched off (switched on at t=0)
edgeShape = 1; % 1 = linear, 2 = sin^2, 3 = tanh, see PulseShape()

% input lens parameters (needed if EFieldModel=='VectorGaussianBeam'        
radiusLensIn = exp(1); % units of w0, beam waist in paraxial approximation
fLensIn = 250*lambda0; % lens focal length
zLensIn = -fLensIn; % z position of lens
nRhoLensIn = 501; % number of points in input lens
nkt = 500; % number of points in k-space (used for mode calculation)

% output lens parameters (needed if calculating transmission or scattering)
radiusLensOut = 125*lambda0; % lens radius
fLensOut = 250*lambda0; % lens focal length 
zLensOut = fLensOut; % z position of lens
nRhoLensOut = 21; % lens made up of grid of (nLensDim x nLensDim) points


%% extra variables

% external B-field - produces a Zeeman shift
muB = 50.*Gamma0;
% muB = 0;

% justReIm = 'Re'; % just calculate real part - this is faster
% justReIm = 'Im'; % just calculate imag part
justReIm = 'ReIm'; % calculate real and imaginary parts


%% run calculations
plotFigs = 1;
displayProgress = 1; % display progress through calculation

if runFromInputVariablesScript == 1
    CalculateScattering
end

% if runFromInputVariablesScript == 1
%     clear 'PTot_varying_w0'
%     clear 'PL_varying_w0'
%     for iw0 = 1:length(w0Vals)        
%         w0 = w0Vals(iw0)        
%         if w0 <= 3*lambda0            
%             EFieldModel = 'VectorGaussianBeam';
%         else
%             EFieldModel = 'ParaxialGaussianBeam';
%         end
%         runFromInputVariablesScript = 1;
%         CalculateScattering
%         PTot_varying_w0(:,iw0) = PTot;
%         PL_varying_w0(:,iw0) = PL;
%     end        
%     
%     figNum = 301;
%     Script_Plot_T_vs_w0;       
%     
% end
