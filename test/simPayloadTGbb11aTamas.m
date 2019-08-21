%simPayloadTGbb11aTamas
% A test file to run the payload simulation for TGbb that is based on 802.11a.
% The simulation setup mostly follow Tamas', especially the AFEs and the implementations
% of DC bias and AGC.
%
%   See also payload11aSimTamas.

% MCS table of 802.11a
% MCS_idx 	MCS 		dataRate (Mbs)
% 0			BPSK 1/2		6
% 1			BPSK 3/4		9
% 2			QPSK 1/2	   12
% 3			QPSK 3/4	   18
% 4			16QAM 1/2	   24
% 5			16QAM 3/4	   36
% 6			64QAM 2/3	   48
% 7			64QAM 3/4      54

%   Logs:
%   20 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Created

clear
close all

addpath('./../lib/filter')
addpath('./../lib/misc')
addpath('./../lib/simulator')
addpath('./../lib/tgbb')
addpath('./../data/freqresponse')

cfgNHT = wlanNonHTConfig;				% 802.11a is NHT waveform
cfgNHT.ChannelBandwidth = 'CBW20';    
cfgNHT.PSDULength = 1000;
cfgNHT.MCS = 0;                       
cfgNHT.Modulation = 'OFDM';             % OFDM modulation

snr = 0:0.5:4; % in dB
whichCase = 'S1_D1'; % {'S1_D1','S1_D2','S3_D1','S3_D2','industrial_d7'}
SamplesPerSymbol = 50; % 20M x 50 = 1 Gsamples/sec
maxNumErrors = 1e0;   % The maximum number of packet errors at an SNR point
maxNumPackets = 1e1; % The Maximum number of packets at an SNR point

disp('> Running payload11aSimTamas...(wait)')
[packetErrorRate,tgbbChannel] = payload11aSimTamas(...
	snr, cfgNHT, whichCase, SamplesPerSymbol, maxNumErrors, maxNumPackets);

if numel(snr)~= 1
	figure;
    semilogy(snr,packetErrorRate(1,:).','-*');
    hold on;
    grid on;
    xlabel('SNR (dB)');
    ylabel('PER');
	title(sprintf('PER (payload) for TGbb with 802.11a, Case: %s',whichCase));
end