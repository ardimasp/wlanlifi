function [packetErrorRate,snr] = simPayloadTGbb11axAWGN(MCS,centerFreq,snr,maxNumErrors,maxNumPackets,varargin)

% restoredefaultpath
% addpath('./../lib/filter')
% addpath('./../lib/misc')
% addpath('./../lib/simulator')
% addpath('./../lib/tgbb')
% addpath('./../data/freqresponse')
% addpath(fullfile(matlabroot,'examples','wlan','main'))

addpath('./../wlan');

% setup configuration for HESU waveform 
cfgHE = wlanHESUConfig;
cfgHE.ChannelBandwidth = 'CBW20';  % Channel bandwidth
% cfgHE.APEPLength = 1e3;            % Payload length in bytes
cfgHE.APEPLength = 4095;            % Payload length in bytes
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

[packetErrorRate] = hesuSimFindFc(snr, cfgHE, centerFreq, maxNumErrors, maxNumPackets,str);

% clear
% close all

% restoredefaultpath
% addpath('./../lib/filter')
% addpath('./../lib/misc')
% addpath('./../lib/simulator')
% addpath('./../lib/tgbb')
% addpath('./../data/freqresponse')
% addpath(fullfile(matlabroot,'examples','wlan','main'))


% % setup configuration for HESU waveform 
% cfgHE = wlanHESUConfig;
% cfgHE.ChannelBandwidth = 'CBW20';  % Channel bandwidth
% % cfgHE.APEPLength = 1e3;            % Payload length in bytes
% cfgHE.APEPLength = 4095;            % Payload length in bytes
% cfgHE.MCS = 0;

% centerFreq = 10.5e6; % in Hz

% snr = 1:0.5:5; % in dB
% maxNumErrors = 1e3;   % The maximum number of packet errors at an SNR point
% maxNumPackets = 1e4; % The Maximum number of packets at an SNR point

% disp('> Running payload11axSim...(wait)')
% [packetErrorRate,tgbbChannel] = hesuSimFindFc(...
%     snr, cfgHE, centerFreq, maxNumErrors, maxNumPackets);

% if numel(snr)~= 1
%     figure;
%     semilogy(snr,packetErrorRate(1,:).','-*');
%     hold on;
%     grid on;
%     xlabel('SNR (dB)');
%     ylabel('PER');
%     title(sprintf('PER (payload) for TGbb with 802.11ax, fc:%2.2f MHz',centerFreq/1e6));
% end
