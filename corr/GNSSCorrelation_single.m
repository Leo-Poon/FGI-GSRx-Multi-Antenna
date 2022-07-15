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
function [tR]= GNSSCorrelation(tR, ch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs code and carrier correlation for GNSS data
%
% Inputs:
%   tR             - Results from signal tracking for one signals
%   ch             - Channel index
%
%   after acq, ant contains the strongest signal, this variable is updates
%   as it tracks
%
% Outputs:
%   tR             - Results from signal tracking for one signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------- DIVERSITY MODE SINGLE -------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set local variables
fid = tR.fid;
%add new channels
fid2 = tR.fid2;
fid3 = tR.fid3;
fid4 = tR.fid4;

loopCnt = tR.channel(ch).loopCnt;


% Read RF data from file
% I dont think anything in tR changes in the different signal acquisition
% this is for consistency
[tR, rawSignal] = getDataForCorrelation(fid,tR,ch);
[tR, rawSignal2] = getDataForCorrelation(fid2,tR,ch);
[tR, rawSignal3] = getDataForCorrelation(fid3,tR,ch);
[tR, rawSignal4] = getDataForCorrelation(fid4,tR,ch);

% Generate finger data
[fingers,tR] = corrFingerGeneration(tR,ch);

%creates temp tR
tR1 = tR;
tR2 = tR;
tR3 = tR;
tR4 = tR;


% Carrier generation + correlation and mixing with code signal
tR1 = carrierMixing(tR1,ch, rawSignal);
tR2 = carrierMixing(tR2,ch, rawSignal2);
tR3 = carrierMixing(tR3,ch, rawSignal3);
tR4 = carrierMixing(tR4,ch, rawSignal4);

%loading tracking results
% in phase
I_P_res(1) = tR1.channel(ch).I_P(end);
I_P_res(2) = tR2.channel(ch).I_P(end);
I_P_res(3) = tR3.channel(ch).I_P(end);
I_P_res(4) = tR4.channel(ch).I_P(end);

Q_P_res(1) = tR1.channel(ch).Q_P(end);
Q_P_res(2) = tR2.channel(ch).Q_P(end);
Q_P_res(3) = tR3.channel(ch).Q_P(end);
Q_P_res(4) = tR4.channel(ch).Q_P(end);

% early results

I_E_res(1) = tR1.channel(ch).I_E(end);
I_E_res(2) = tR2.channel(ch).I_E(end);
I_E_res(3) = tR3.channel(ch).I_E(end);
I_E_res(4) = tR4.channel(ch).I_E(end);

Q_E_res(1) = tR1.channel(ch).Q_E(end);
Q_E_res(2) = tR2.channel(ch).Q_E(end);
Q_E_res(3) = tR3.channel(ch).Q_E(end);
Q_E_res(4) = tR4.channel(ch).Q_E(end);

%late results

I_L_res(1) = tR1.channel(ch).I_L(end);
I_L_res(2) = tR2.channel(ch).I_L(end);
I_L_res(3) = tR3.channel(ch).I_L(end);
I_L_res(4) = tR4.channel(ch).I_L(end);

Q_L_res(1) = tR1.channel(ch).Q_L(end);
Q_L_res(2) = tR2.channel(ch).Q_L(end);
Q_L_res(3) = tR3.channel(ch).Q_L(end);
Q_L_res(4) = tR4.channel(ch).Q_L(end);

% time to combine the results into complex numbers

earlyRes = I_E_res + Q_E_res.*j;
promptRes = I_P_res + Q_P_res.*j;
lateRes = I_L_res + Q_L_res.*j;

%returns max value and index representing antenna

[maxEarly, i_early] = max(earlyRes);
[maxPrompt, i_prompt] = max(promptRes);
[maxLate, i_late] = max(lateRes);

[maxVals, index] = max([maxEarly maxPrompt maxLate]);

switch index
    case 1
        antSel = i_early;
    case 2
        antSel = i_prompt;
    case 3
        antSel = i_late;
    otherwise
        antSel = i_prompt;
end

switch antSel
    case 1
        tR = tR1;
    case 2
        tR = tR2;
    case 3
        tR = tR3;
    case 4
        tR = tR4;
    otherwise
        tR = tR1;
end

tR.channel(ch).antNum(end+1) = antSel;

% Check if user have requested multi correlator tracking
% This should be disabled
if(tR.enableMultiCorrelatorTracking)    
    tR = multiFingerTracking(tR,ch,fingers); % Generate fingers for multi finger tracking
    if(mod(tR.channel(ch).loopCnt,tR.multiCorrelatorTrackingRate) == 0)
        % Plot output
        tR.channel(ch) = plotMultiFingerTracking(tR.channel(ch));
    end
end

loopCnt = tR.channel(ch).loopCnt;

% Check where we are in data file
tR.channel(ch).absoluteSample(loopCnt) =(ftell(fid))/(tR.sampleSize/8);
tR.channel(ch).prevAbsoluteSample =tR.channel(ch).absoluteSample(loopCnt);




