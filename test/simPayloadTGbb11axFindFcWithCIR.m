function [packetErrorRate,snr,tgbbChannel] = simPayloadTGbb11axFindFcWithCIR(MCS,centerFreq,snr,maxNumErrors,maxNumPackets,varargin)

restoredefaultpath
addpath('./../lib/filter')
addpath('./../lib/misc')
addpath('./../lib/simulator')
addpath('./../lib/tgbb')
addpath('./../data/freqresponse')
addpath(fullfile(matlabroot,'examples','wlan','main'))

% setup configuration for HESU waveform 
cfgHE = wlanHESUConfig;
cfgHE.ChannelBandwidth = 'CBW20';  % Channel bandwidth
cfgHE.APEPLength = 1e3;            % Payload length in bytes
% cfgHE.APEPLength = 4095;            % Payload length in bytes
cfgHE.MCS = MCS;

% centerFreq = 10.5e6; % in Hz

% snr = 1:0.5:5; % in dB
% maxNumErrors = 1e3;   % The maximum number of packet errors at an SNR point
% maxNumPackets = 1e4; % The Maximum number of packets at an SNR point

disp('> Running payload11axSim...(wait)')
if ~isempty(varargin)
    str = varargin{1};
else
    str = '';
end

whichCase = './../data/tgbbcirs/simulation scenario hospital ward/individual cirs/optical cirs/s12/D5/Run1.mat';

[packetErrorRate,tgbbChannel] = hesuSimFindFcWithCIR(snr, cfgHE, centerFreq, whichCase, maxNumErrors, maxNumPackets,str);

