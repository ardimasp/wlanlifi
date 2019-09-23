function [packetErrorRate,varargout] = hesuSimFindFc(...
    snr, cfgHE, centerFreq, maxNumErrors, maxNumPackets, varargin)
%hesuSim simulates the 802.11ax-based TGbb

%   Logs:
%   12 Sept 19, ardimas (ardimasandipurwita@outlook.com):
%       Created

% Specify waveform parameters
numTx = cfgHE.NumTransmitAntennas;
numRx = 1; % Number of receive antennas

% setup waveform recovery parameters
chanBW = cfgHE.ChannelBandwidth; % Assume channel bandwidth is known
sr = wlanSampleRate(cfgHE); % Sample rate
indHESU = wlanFieldIndices(cfgHE);

% specify propagation channel
tgbbChannel                     = wlanTGbbChannel;
tgbbChannel.SymbolRate          = sr;
tgbbChannel.SamplesPerSymbol    = floor(1e9/sr);
tgbbChannel.beta                = 0.1;
tgbbChannel.awgnonlymode        = true;
tgbbChannel.simMethod           = 'snrOptical';

tgbbChannel.isApplyNoise        = true;
tgbbChannel.isShotNoise         = true;
% tgbbChannel.isApplyNoise        = false;
% tgbbChannel.isShotNoise         = false;

tgbbChannel.Offset              = centerFreq - sr/2;

% tgbbChannel.methodDCBias        = 'max';
% tgbbChannel.methodAGC           = 'max';
% tgbbChannel.whichFE             = 'tamas';

if ~isempty(varargin)
    filename = ['hesuFindFc_' 'MCS' num2str(cfgHE.MCS) '_fc' num2str(centerFreq/6) 'MHz' '_' varargin{1} '.mat'];
else
    filename = ['hesuFindFc_' 'MCS' num2str(cfgHE.MCS) '_fc' num2str(centerFreq/6) 'MHz' '.mat'];
end

numSNR = numel(snr); % Number of SNR points
packetErrorRate = zeros(1,numSNR);

%parfor isnr = 1:numSNR % Use 'parfor' to speed up the simulation
for isnr = 1:numSNR
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',1412);
    stream.Substream = isnr;
    RandStream.setGlobalStream(stream);

    % Loop to simulate multiple packets
    numPacketErrors = 0;
    numPkt = 1; % Index of packet transmitted
    while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
        % Generate a packet with random PSDU
        psduLength = getPSDULength(cfgHE); % PSDU length in bytes
        txPSDU = randi([0 1],psduLength*8,1);
        txWaveform = wlanWaveformGenerator(txPSDU,cfgHE);

        % padding symbol to take into account the delay due to the filters
        % expecting that the delay will be less than 10 baseband symbols
        txWaveformPadded = [txWaveform; zeros(tgbbChannel.SamplesPerSymbol*10,1)];

        tgbbChannel.SNRdB = snr(isnr);
        rx = tgbbChannel(txWaveformPadded);

        % The recovery configuration object is used to get the start and end
        % indices of the pre-HE-SIG-B field.
        ind = wlanFieldIndices(cfgHE);

        % Minimum packet length is 10 OFDM symbols
        lstfLength = double(ind.LSTF(2));
        minPktLen = lstfLength*5; % Number of samples in L-STF

        rxWaveLen = size(rx,1);

        searchOffset = 0; % Offset from start of waveform in samples
        while (searchOffset + minPktLen) <= rxWaveLen

            % Packet detection
            pktOffset = wlanPacketDetect(rx,chanBW,searchOffset);

            % Adjust packet offset
            pktOffset = searchOffset + pktOffset;
            % if isempty(pktOffset) || (pktOffset + ind.LSIG(2) > rxWaveLen)
            if isempty(pktOffset) || (pktOffset<0) ||  (pktOffset + ind.LSIG(2) > rxWaveLen)

                isBreak = 1;

                break; 
            else
                isBreak = 0;
            end

            % Coarse frequency offset estimation and correction using L-STF
            LSTF = rx(pktOffset+(ind.LSTF(1):ind.LSTF(2)), :);
            coarseFreqOffset = wlanCoarseCFOEstimate(LSTF,chanBW);
            rx = helperFrequencyOffset(rx,sr,-coarseFreqOffset);

            % Symbol timing synchronization
            searchBufferLLTF = rx(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            pktOffset = pktOffset+wlanSymbolTimingEstimate(searchBufferLLTF,chanBW);

            %% ADDING THIS
            % If offset is without bounds of waveform  skip samples and continue
            % searching within remainder of the waveform
            if (pktOffset<0) || ((pktOffset+indHESU.HEData(2))>rxWaveLen)
                searchOffset = pktOffset+double(ind.LSTF(2))+1;
                continue;
            end


            % Fine frequency offset estimation and correction using L-STF
            rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
            fineFreqOffset = wlanFineCFOEstimate(rxLLTF,chanBW);
            rx = helperFrequencyOffset(rx,sr,-fineFreqOffset);

            % Display estimated carrier frequency offset
            cfoCorrection = coarseFreqOffset + fineFreqOffset; % Total CFO

            break; % Front-end processing complete, stop searching for a packet
        end

        if isBreak % If error in packet detection
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end

        %% L-LTF Channel and Noise Power Estimation
        % Extract and demodulate the L-LTF and perform channel and noise power estimation. The L-LTF channel and noise power estimates are used to decode the pre HE-LTF.
        rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        lltfDemod = wlanHEDemodulate(rxLLTF,'L-LTF',chanBW);
        lltfChanEst = wlanLLTFChannelEstimate(lltfDemod,chanBW);
        noiseVar = helperNoiseEstimate(lltfDemod);

        % L-SIG and RL-SIG Decoding
        % L-SIG and RL-SIG are processed simultaneously. This results in some improvement in terms of SNR
        
        % Extract L-SIG and RL-SIG fields
        rxLSIG = rx(pktOffset+(ind.LSIG(1):ind.RLSIG(2)),:);

        % OFDM demodulate
        helsigDemod = wlanHEDemodulate(rxLSIG,'L-SIG',chanBW);

        % Estimate CPE and phase correct symbols
        helsigDemod = preHECommonPhaseErrorTracking(helsigDemod,lltfChanEst,'L-SIG',chanBW);

        % Estimate channel on extra 4 subcarriers per subchannel and create full
        % channel estimate
        preheInfo = wlanHEOFDMInfo('L-SIG',chanBW);
        preHEChanEst = preHEChannelEstimate(helsigDemod,lltfChanEst,preheInfo.NumSubchannels);

        % Average L-SIG and RL-SIG before equalization
        helsigDemod = mean(helsigDemod,2);

        % Equalize data carrying subcarriers, merging 20 MHz subchannels
        [eqLSIGSym,csi] = preHESymbolEqualize(helsigDemod(preheInfo.DataIndices,:,:), ...
            preHEChanEst(preheInfo.DataIndices,:,:),noiseVar,preheInfo.NumSubchannels);

        % Decode L-SIG field
        [rxLSIGBits,failCheck,lsigInfo] = wlanLSIGBitRecover(eqLSIGSym,noiseVar,csi);
        [~,txLSIGBits] = wlan.internal.heLSIG(cfgHE);

        if failCheck && ~any(biterr(txLSIGBits,rxLSIGBits))
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end

        % Get the length information from the recovered L-SIG bits and update the
        % L-SIG length property of the recovery configuration object
        lsigLength = lsigInfo.Length;

        % Packet format detection
        rxSIGA = rx(pktOffset+(ind.HESIGA(1):ind.HESIGA(2)),:);
        pktFormat = hePacketFormat(rxSIGA,preHEChanEst,lsigLength,noiseVar,chanBW);

        % Update the packet format in the recovery object and recompute the field
        % indices
        ind = wlanFieldIndices(cfgHE);

        % HE-SIG-A decoding
        rxSIGA = rx(pktOffset+(ind.HESIGA(1):ind.HESIGA(2)),:);
        sigaDemod = wlanHEDemodulate(rxSIGA,'HE-SIG-A',chanBW);
        hesigaDemod = preHECommonPhaseErrorTracking(sigaDemod,preHEChanEst,'HE-SIG-A',chanBW);

        % Equalize data carrying subcarriers, merging 20 MHz subchannels
        preheInfo = wlanHEOFDMInfo('HE-SIG-A',chanBW);
        [eqSIGASym,csi] = preHESymbolEqualize(hesigaDemod(preheInfo.DataIndices,:,:), ...
                                              preHEChanEst(preheInfo.DataIndices,:,:), ...
                                              noiseVar,preheInfo.NumSubchannels);
        % Recover HE-SIG-A bits
        [rxSIGABits,failCRC] = wlanHESIGABitRecover(eqSIGASym,noiseVar,csi);
        [~,txSIGABits] = wlan.internal.heSIGA(cfgHE);

        % Perform the CRC on HE-SIG-A bits
        if failCRC && ~any(biterr(txSIGABits,rxSIGABits))
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end

        %% Interpret Recovered HE-SIG-A bits
        ind = wlanFieldIndices(cfgHE); % Update field indices

        numUsers = 1;

        constellationDiagram = heSigEqualizeSetupPlots(numUsers);

        % HE-LTF demodulation and channel estimation
        rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
        heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfgHE);
        [chanEst,pilotEst] = heLTFChannelEstimate(heltfDemod,cfgHE);

        % HE-Data demodulate
        rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
        demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE);

        % Pilot phase tracking. Average single-stream pilot estimates over
        % symbols (2nd dimension)
        pilotEstTrack = mean(pilotEst,2);
        demodSym = heCommonPhaseErrorTracking(demodSym,pilotEstTrack,cfgHE);

        % Estimate noise power in HE fields
        preheInfo = wlanHEOFDMInfo('HE-Data',cfgHE);
        demodPilotSym = demodSym(preheInfo.PilotIndices,:,:);
        nVarEst = heNoiseEstimate(demodPilotSym,pilotEst,cfgHE);

        % Equalize
        [eqSym,csi] = heEqualizeCombine(demodSym,chanEst,nVarEst,cfgHE);

        % Discard pilot subcarriers
        eqSymUser = eqSym(preheInfo.DataIndices,:,:);
        csiData = csi(preheInfo.DataIndices,:);

        % Demap and decode bits
        rxPSDU = wlanHEDataBitRecover(eqSymUser,nVarEst,csiData,cfgHE);

        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr([txPSDU], [rxPSDU]));
        numPacketErrors = numPacketErrors+packetError;
        numPkt = numPkt+1;

    end

    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(isnr) = numPacketErrors/(numPkt-1);
    disp(['> MCS ' num2str(cfgHE.MCS) ','...
          ' SNR ' num2str(snr(isnr)) ','...
          ' payloadLen ' num2str(cfgHE.APEPLength) ','...
          ' fc ' num2str(tgbbChannel.Fc/1e6) ' MHz,'...
          ' completed after ' num2str(numPkt-1) ' packets,'...
          ' PER:' num2str(packetErrorRate(isnr))]);

    fc = tgbbChannel.Fc;
    save(filename,'cfgHE','snr','fc','packetErrorRate','isnr','numSNR');

end

varargout{1} = tgbbChannel;

end