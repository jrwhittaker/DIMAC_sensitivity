function [Y] = pulsatile_scaling_function(varargin)
% Simulate a periodic physiological signal given some baseline value and
% maximume percentage signal change.
%
% Usage:
%
%   [Y] = pulsatile_scaling_function(varargin)
%
% Compulsory arguments: S0  - Baseline signal value.
%                       DS  - Maximum signal value as a percentage change
%                             from baseline (e.g. 1 = 100% change, i.e. max
%                             signal is double of baseline = S0 + DS*S0)
%
% Optional parameters:  N   - Number of samples (default = 64).
%                       K   - Number of Fourier terms
%
% Output:   A vector of length N representing a single cycle of a periodic
%           pulsatile signal. 
%
% Author: JR Whittaker 4th June 2020. 

%% Model parameters
% Read in mandatory input parameters
p = inputParser;
addRequired(p,  'S0'); % Signal baseline
addRequired(p,  'DS'); % Max signal value

% Read in optional user input parameters
addOptional(p,  'N'     , 64 ,  @isnumeric); % Number of samples
addOptional(p,  'K'     , 3 ,   @isnumeric); % number of Fouerier terms

% Parse
parse(p,varargin{:});
args = p.Results;

% Model parameters
S0  =   args.S0;
DS  =   args.DS;
N   =   args.N;
K   =   args.K;

% Max 10 terms
if (K > 10)
    K=10;
end

% Fourier coefficients
% ... taken from Yang et al. BioMed Eng OnLine (2019) 
% https://doi.org/10.1186/S12938-018-0620-3
dc=0.9490;
a=[0.5150; 0.2893; 0.2838; 0.1498; 0.1227; ...
    0.1177; 0.0569; 0.0105; 0.0008; 0.0056];
b=[0.5873; 0.1871; 0.0459; -0.0032; -0.1109; ...
    -0.0743; -0.0196; -0.0086; -0.0043; -0.0007];

%% Evaluate model

for k=1:K
    
    freq=k/N;
    pulseComps(:,k)=a(k)*cos(2*pi*freq*[1:N])'+b(k)*sin(2*pi*freq*[1:N])';
    
end

pulseS=dc+sum(pulseComps,2);

% normalise at start at wave foot
[~,minIdx]=min(pulseS);
pulseS=pulseS([minIdx:end 1:(minIdx-1)],:);
pulseS=pulseS-min(pulseS);
pulseS=pulseS./max(pulseS);

Y = S0 + DS*pulseS*S0;






