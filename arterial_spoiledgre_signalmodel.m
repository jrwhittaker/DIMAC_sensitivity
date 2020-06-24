function [S] = arterial_spoiledgre_signalmodel(varargin)
% Calculate the steady-state spoiled GRE signal from a single voxel
% orientated perpendicularly to an artery that is contained with some
% partial volume. The model has two compartments, intravascular
% (blood volume) and extravascular (CSF).
%
% Usage:
%
%   [S] = arterial_spoiledgre_signalmodel(varargin)
%
% Compulsory arguments: PV     - Arterial blood partial volume (0 - 1).
%                       CBFV   - Arterial blood flow velocity (cm/s).
%                       FA     - Flip angle (degrees).
%                       TR     - Repetition time (s).
%                       TE     - Echo time (s)
%
% Optional parameters:  ST     - Slice thickness (cm).
%                       M0     - Equilibrium magnetisation (a.u).
%                       T1i    - Longitudinal relaxation time of
%                                arterial blood (s).
%                       T1e    - Longitudinal relaxation time of
%                                brain csf (s).
%                       RHOi   - Relative proton density of blood.
%                       RHOe   - Relative proton density of CSF.
%                       T2Si   - Apparent relaxation time of arterial
%                                blood (s).
%                       T2Se   - Apparent relxation time of CSF (s).
%                       SNR    - Signal to noise ratio.
%
% Output:    A single value of the MR magnitude signal.
%
% Author:    JR Whittaker, 4th June 2020.

%% Model parameters
% Read in mandatory input parameters
p = inputParser;
addRequired(p,      'PV'); % Artery partial volume
addRequired(p,   'CBFV');  % Artery blood flow velocity
addRequired(p,     'FA');  % Flip angle
addRequired(p,     'TR');  % Repetition time
addRequired(p,      'TE'); % Echo time

% Read in optional user input parameters
addOptional(p, 'ST'    , 1 ,       @isnumeric); % Slice thickness
addOptional(p, 'M0'    , 1 ,       @isnumeric); % Equilibrium magnetisation
addOptional(p, 'T1i'   , 1.660 ,   @isnumeric); % Arterial blood T1
addOptional(p, 'T1e'   , 3.820 ,   @isnumeric); % CSF T1
addOptional(p, 'RHOi'  , 0.850 ,    @isnumeric); % Arterial blood proton density
addOptional(p, 'RHOe'  , 1 ,       @isnumeric); % CSF proton density
addOptional(p, 'T2Si'  , 0.0480 ,  @isnumeric); % Arterial blood T2*
addOptional(p, 'T2Se'  , 0.4000 ,  @isnumeric); % CSF T2*
addOptional(p, 'SNR'   , Inf    ,  @isnumeric); % SNR

% N.B. Relaxation times and proton densities
%
% T1 of arterial blood at 3T is 1664 ms at Hct=0.42 - Lu et al, 2004. 
% Determining the longitudinal relxation time (T1) of blood at 3.0 Tesla. MRM.
%
% T1 of CSF at 3T is 3817 ms - Lu et al, 2005. Routine clinical brain MRI
% sequences for use at 3.0 Tesla. JMRI. 
%
% T2* of arterial blood is 48.4 ms for 98% oxygenation and Hct=0.44 - Zhao
% et al, 2007. Oxygenation and hematocrit dependence of transverse
% relaxation rates of blood at 3T. MRM.
%
% T2* of CSF is 400 ms
%
% Relative proton density of arterial blood according to equation:
% Cbl = 0.95 g/ml - (1 - Hct)*Hct. So for Hct=0.44, Cbl=0.85

% Parse
parse(p,varargin{:});
args = p.Results;

% Model parameters
PV      =   args.PV;
CBFV    =   args.CBFV;
FA      =   args.FA;
TR      =   args.TR;
TE      =   args.TE;
ST      =   args.ST;
M0      =   args.M0;
T1i     =   args.T1i;
T1e     =   args.T1e;
RHOi    =   args.RHOi;
RHOe    =   args.RHOe;
T2Si    =   args.T2Si;
T2Se    =   args.T2Se;
SNR     =   args.SNR;

%% Checks
% Convert flip angle to radians
FA = deg2rad(FA);
% Make sure 0 <= PV <= 1
if (PV < 0), PV=0; end
if (PV > 1), PV=1; end

%% Functions
% Function handles
longmag = @inflow_longitudinal_magnetisation;
transmag = @transverse_magnetisation;

%% Evaluation model
% Intravascular signal
long_params_i = [ST TR M0 T1i FA CBFV];
Mzi = longmag(long_params_i);
trans_params_i = [Mzi PV RHOi TE T2Si FA];
Mxyi = transmag(trans_params_i);

% Extravascular signal
long_params_e = [ST TR M0 T1e FA 0];
Mze = longmag(long_params_e);
trans_params_e = [Mze (1-PV) RHOe TE T2Se FA];
Mxye = transmag(trans_params_e);

% Total signal = intravascular signal + extravascular signal
S = Mxyi + Mxye;
% Sampled signal = total signal + noise
S = abs(S + (M0/SNR)*randn);


end % main function

%% Local function definitions
% Function to calculate inflow weighted longitudinal magnetisation
function Mz = inflow_longitudinal_magnetisation(params)

ST = params(1); % Slice thickness
TR = params(2); % Repetition time
M0 = params(3); % Equilibrium magnetisation
T1 = params(4); % Longitudinal relaxation
FA = params(5); % Flip angle
V  = params(6); % Spin velocity

Vc = ST/TR; % critical velocity

q = exp(-TR/T1)*cos(FA);

if (V == 0)
    
    Mz = (M0*(1-exp(-TR/T1))/(1-q));
    
elseif (V < Vc)
    
    Mz0 = (M0*(1-exp(-TR/T1))/(1-q));
    Mz = Mz0+(M0 - Mz0)*(1-q^(Vc/V))/((Vc/V)*(1-q));
    
else
    
    Mz = M0;
    
end

end % inflow_longitudinal_magnetisation function

% Function to calculate the measured trasnverse magnetisation
function Mxy = transverse_magnetisation(params)

MZ  = params(1); % Longitudinal magnetisation
PV  = params(2); % Partial volume
RHO = params(3); % Spin density
TE  = params(4); % Echo time
T2S = params(5); % Apparent relaxation time;
FA  = params(6); % Flip angle

Mxy = MZ*PV*RHO*exp(-TE/T2S)*sin(FA);

end % transverse_magnetisation function
