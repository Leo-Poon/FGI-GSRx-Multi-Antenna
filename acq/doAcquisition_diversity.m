%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright 2015-2021 Finnish Geospatial Research Institute FGI, National
%% Land Survey of Finland. This file is part of FGI-GSRx software-defined
%% receiver. FGI-GSRx is a free software: you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as published
%% by the Free Software Foundation, either version 3 of the License, or any
%% later version. FGI-GSRx software receiver is distributed in the hope
%% that it will be useful, but WITHOUT ANY WARRANTY, without even the
%% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%% See the GNU General Public License for more details. You should have
%% received a copy of the GNU General Public License along with FGI-GSRx
%% software-defined receiver. If not, please visit the following website 
%% for further information: https://www.gnu.org/licenses/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acqResults = doAcquisition(allSettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function acquires all enabled signals
%
% Inputs: 
%   allSettings         - Receiver settings

% Outputs:
%   acqResults          - Acquisition results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------- DIVERSITY MODE --------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise acquisition structure
acqResults = initAcquisition(allSettings);

% Loop over all signals
for i = 1:allSettings.sys.nrOfSignals
    
    % Temporary variables
    signal = allSettings.sys.enabledSignals{i};
    msToSkip = allSettings.sys.msToSkip;
    
    % Extract block of parameters for one signal from settings
    signalSettings = allSettings.(signal);
    
    % Read RF Data
    % Added new functions for channels
    [pRfData,sampleCount] = getDataForAcquisition(signalSettings,100);
    [pRfData2,sampleCount] = getDataForAcquisition2(signalSettings,100);  
    [pRfData3,sampleCount] = getDataForAcquisition3(signalSettings,100);  
    [pRfData4,sampleCount] = getDataForAcquisition4(signalSettings,100);  
    
    % Execute acquisition for one signal
    acqResults.(signal) = acquireSignal(pRfData,signalSettings);
    acqResults(2).(signal) = acquireSignal(pRfData2,signalSettings);
    acqResults(3).(signal) = acquireSignal(pRfData3,signalSettings);
    acqResults(4).(signal) = acquireSignal(pRfData4,signalSettings);
    
end





