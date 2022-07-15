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

% ----------------------- DIVERSITY MODE HYBRID -------------------------%

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

%returns max value and index
%finds the largest early, prompt and late
[maxEarly, i_early] = maxk(earlyRes,2);
[maxPrompt, i_prompt] = maxk(promptRes,2);
[maxLate, i_late] = maxk(lateRes,2);

%finds largest 2 values
[maxVals, index] = maxk([maxEarly maxPrompt maxLate],2);

% [maxEarly maxPrompt maxLate] is defined as Early[1,2], Prompt[1,2], 
% late [1,2] where [1,2] is [largest, 2nd largest]. Index from 1-6

% we do not care if it is early, prompt or late, just which ant has the
% strongest signal

% largest and magnitude
%it should not go to 2,4,6 as those are 2nd largest but in here incase
switch index(1) 
    case 1
        strongest_ant = i_early(1);
        strongest_mag = maxEarly(1);
    case 2
        strongest_ant = i_early(2);
        strongest_mag = maxEarly(2);
    case 3
        strongest_ant = i_prompt(1);
        strongest_mag = maxPrompt(1);
    case 4
        strongest_ant = i_prompt(2);
        strongest_mag = maxPrompt(2);
    case 5
        strongest_ant = i_late(1);
        strongest_mag = maxLate(1);
    case 6
        strongest_ant = i_late(2);
        strongest_mag = maxLate(2);
    otherwise
        strongest_ant = i_prompt(1);
        strongest_mag = maxPrompt(1);
end

%the 2nd strongest, this can be any value except above case
switch index(2)
    case 1
        strongest_ant2 = i_early(1);
        strongest_mag2 = maxEarly(1);
    case 2
        strongest_ant2 = i_early(2);
        strongest_mag2 = maxEarly(2);
    case 3
        strongest_ant2 = i_prompt(1);
        strongest_mag2 = maxPrompt(1);
    case 4
        strongest_ant2 = i_prompt(2);
        strongest_mag2 = maxPrompt(2);
    case 5
        strongest_ant2 = i_late(1);
        strongest_mag2 = maxLate(1);
    case 6
        strongest_ant2 = i_late(2);
        strongest_mag2 = maxLate(2);
    otherwise
        strongest_ant2 = i_prompt(1);
        strongest_mag2 = maxPrompt(1);
end

%if the 2 strongest antenna are the same or the 2nd largest power is less than
%largest/sqrt(2) use 1 antenna 

%sqrt(2) can be changed to another number

if (strongest_ant == strongest_ant2) || (abs(strongest_mag2) <= abs(strongest_mag)/(sqrt(2)))
    %use 1 active antenna configuration
    switch strongest_ant
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
    antSel = strongest_ant;
elseif (strongest_ant == 4 && strongest_ant2 == 1) || (strongest_ant == 1 && strongest_ant2 == 4)
    %calculates by adding the antenna number and dividing by 2 except in
    %1,4 or 4,1 where it overflows/underflow, called antNum 4.5
    
    %consistent between runs so we just pick 1

    tR.channel(ch).earlyCode = tR4.channel.earlyCode;
    tR.channel(ch).lateCode = tR4.channel.lateCode;
    tR.channel(ch).promptCode = tR4.channel.promptCode;
    tR.channel(ch).twoChipEarlyCode = tR4.channel(ch).twoChipEarlyCode;

    %adding stuff together now

    tR.channel(ch).qBasebandSignal = (tR4.channel(ch).qBasebandSignal + tR1.channel(ch).qBasebandSignal)/2;
    tR.channel(ch).iBasebandSignal = (tR4.channel(ch).iBasebandSignal + tR1.channel(ch).iBasebandSignal)/2;
    tR.channel(ch).I_E = (tR4.channel(ch).I_E + tR1.channel(ch).I_E)/2;
    tR.channel(ch).I_P = (tR4.channel(ch).I_P + tR1.channel(ch).I_P)/2;
    tR.channel(ch).I_L = (tR4.channel(ch).I_L + tR1.channel(ch).I_L)/2;
    tR.channel(ch).Q_E = (tR4.channel(ch).Q_E + tR1.channel(ch).Q_E)/2;
    tR.channel(ch).Q_P = (tR4.channel(ch).Q_P + tR1.channel(ch).Q_P)/2;
    tR.channel(ch).Q_L = (tR4.channel(ch).Q_L + tR1.channel(ch).Q_L)/2;
    tR.channel(ch).I_E_E = (tR4.channel(ch).I_E_E + tR1.channel(ch).I_E_E)/2;
    tR.channel(ch).Q_E_E = (tR4.channel(ch).Q_E_E + tR1.channel(ch).Q_E_E)/2;
    antSel = 4.5;
else
    %other cases should end in .5
    
    antSel = (strongest_ant + strongest_ant2)/2;
    switch antSel
        case 1.5
        %consistent between runs
            tR.channel(ch).earlyCode = tR1.channel.earlyCode;
            tR.channel(ch).lateCode = tR1.channel.lateCode;
            tR.channel(ch).promptCode = tR1.channel.promptCode;
            tR.channel(ch).twoChipEarlyCode = tR1.channel(ch).twoChipEarlyCode;

        %adding stuff together now
            tR.channel(ch).qBasebandSignal = (tR1.channel(ch).qBasebandSignal + tR2.channel(ch).qBasebandSignal)/2;
            tR.channel(ch).iBasebandSignal = (tR1.channel(ch).iBasebandSignal + tR2.channel(ch).iBasebandSignal)/2;
            tR.channel(ch).I_E = (tR1.channel(ch).I_E + tR2.channel(ch).I_E)/2;
            tR.channel(ch).I_P = (tR1.channel(ch).I_P + tR2.channel(ch).I_P)/2;
            tR.channel(ch).I_L = (tR1.channel(ch).I_L + tR2.channel(ch).I_L)/2;
            tR.channel(ch).Q_E = (tR1.channel(ch).Q_E + tR2.channel(ch).Q_E)/2;
            tR.channel(ch).Q_P = (tR1.channel(ch).Q_P + tR2.channel(ch).Q_P)/2;
            tR.channel(ch).Q_L = (tR1.channel(ch).Q_L + tR2.channel(ch).Q_L)/2;
            tR.channel(ch).I_E_E = (tR1.channel(ch).I_E_E + tR2.channel(ch).I_E_E)/2;
            tR.channel(ch).Q_E_E = (tR1.channel(ch).Q_E_E + tR2.channel(ch).Q_E_E)/2;
        case 2.5
        %consistent between runs
            tR.channel(ch).earlyCode = tR2.channel.earlyCode;
            tR.channel(ch).lateCode = tR2.channel.lateCode;
            tR.channel(ch).promptCode = tR2.channel.promptCode;
            tR.channel(ch).twoChipEarlyCode = tR2.channel(ch).twoChipEarlyCode;

        %adding stuff together now
            tR.channel(ch).qBasebandSignal = (tR2.channel(ch).qBasebandSignal + tR3.channel(ch).qBasebandSignal)/2;
            tR.channel(ch).iBasebandSignal = (tR2.channel(ch).iBasebandSignal + tR3.channel(ch).iBasebandSignal)/2;
            tR.channel(ch).I_E = (tR2.channel(ch).I_E + tR3.channel(ch).I_E)/2;
            tR.channel(ch).I_P = (tR2.channel(ch).I_P + tR3.channel(ch).I_P)/2;
            tR.channel(ch).I_L = (tR2.channel(ch).I_L + tR3.channel(ch).I_L)/2;
            tR.channel(ch).Q_E = (tR2.channel(ch).Q_E + tR3.channel(ch).Q_E)/2;
            tR.channel(ch).Q_P = (tR2.channel(ch).Q_P + tR3.channel(ch).Q_P)/2;
            tR.channel(ch).Q_L = (tR2.channel(ch).Q_L + tR3.channel(ch).Q_L)/2;
            tR.channel(ch).I_E_E = (tR2.channel(ch).I_E_E + tR3.channel(ch).I_E_E)/2;
            tR.channel(ch).Q_E_E = (tR2.channel(ch).Q_E_E + tR3.channel(ch).Q_E_E)/2;
        case 3.5
        %consistent between runs
            tR.channel(ch).earlyCode = tR3.channel.earlyCode;
            tR.channel(ch).lateCode = tR3.channel.lateCode;
            tR.channel(ch).promptCode = tR3.channel.promptCode;
            tR.channel(ch).twoChipEarlyCode = tR3.channel(ch).twoChipEarlyCode;

        %adding stuff together now
            tR.channel(ch).qBasebandSignal = (tR3.channel(ch).qBasebandSignal + tR4.channel(ch).qBasebandSignal)/2;
            tR.channel(ch).iBasebandSignal = (tR3.channel(ch).iBasebandSignal + tR4.channel(ch).iBasebandSignal)/2;
            tR.channel(ch).I_E = (tR3.channel(ch).I_E + tR4.channel(ch).I_E)/2;
            tR.channel(ch).I_P = (tR3.channel(ch).I_P + tR4.channel(ch).I_P)/2;
            tR.channel(ch).I_L = (tR3.channel(ch).I_L + tR4.channel(ch).I_L)/2;
            tR.channel(ch).Q_E = (tR3.channel(ch).Q_E + tR4.channel(ch).Q_E)/2;
            tR.channel(ch).Q_P = (tR3.channel(ch).Q_P + tR4.channel(ch).Q_P)/2;
            tR.channel(ch).Q_L = (tR3.channel(ch).Q_L + tR4.channel(ch).Q_L)/2;
            tR.channel(ch).I_E_E = (tR3.channel(ch).I_E_E + tR4.channel(ch).I_E_E)/2;
            tR.channel(ch).Q_E_E = (tR3.channel(ch).Q_E_E + tR4.channel(ch).Q_E_E)/2;
        otherwise
            tR = tR1;
    end
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




