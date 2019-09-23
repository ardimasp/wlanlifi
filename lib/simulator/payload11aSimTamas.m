function [packetErrorRate,varargout] = payload11aSimTamas(...
    snr, cfgNHT, whichCase, SamplesPerSymbol, maxNumErrors, maxNumPackets)
%payload11aSimTamas Simulates the 802.11a-based TGbb following Tamas' setup.
% That is, the FEs are based on Tamas'. The implementations of DC and AGC also follow
% Tamas.

%   Logs:
%   20 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Created

% Get symbol rate from cfgNHT
fs = wlanSampleRate(cfgNHT);

% Instantiate the tgbb object
tgbbChannel                     = wlanTGbbChannel;
tgbbChannel.SymbolRate          = fs;
tgbbChannel.useData             = true;
tgbbChannel.whichCase           = whichCase;
tgbbChannel.isIncludeDelay      = true;
tgbbChannel.isApplyNoise        = true;
tgbbChannel.SamplesPerSymbol    = SamplesPerSymbol;

tgbbChannel.isShotNoise         = true;
tgbbChannel.beta                = 0.1;
tgbbChannel.methodDCBias        = 'max';
tgbbChannel.methodAGC           = 'max';
tgbbChannel.whichFE             = 'tamas';

numSNR = numel(snr); % Number of SNR points
packetErrorRate = zeros(1,numSNR);

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgNHT);

% Get the number of occupied subcarriers and FFT length
ofdmInfo = wlan.internal.wlanGetOFDMConfig(cfgNHT.ChannelBandwidth, 'Long', 'Legacy');
Nst = numel(ofdmInfo.DataIndices)+numel(ofdmInfo.PilotIndices); % Number of occupied subcarriers

% Generate the legacy preamble fields
lstf = wlanLSTF_custom(cfgNHT);
lltf = wlanLLTF_custom(cfgNHT);
lsig = wlanLSIG_custom(cfgNHT);
preamble = [lstf; lltf; lsig];

% NHT params for OFDM format
giType = 'Long'; % Always
FFTLen = 64;
wlength = 2*ceil(1e-7*fs/2);
numSamplesBeforeShortGI = size(preamble,1);
bLen = wlength/2; % Number of samples overlap at the end of the packet
aLen = bLen-1; % Number of samples overlap at start of packet
s = validateConfig(cfgNHT);
numPktSamples = real(s.NumPPDUSamples);

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
        % Generate a packet waveform
        inpPSDU = randi([0 1], cfgNHT.PSDULength*8, 1); % PSDULength in bytes
        
        % tx = wlanWaveformGenerator(inpPSDU, cfgNHT);
        % Copy-paste from wlanWaveformGenerator
        data = wlanNonHTData_custom(inpPSDU,cfgNHT);
        % Construct packet from preamble and data
        packet = [preamble; data];
        windowedPacket = wlan.internal.wlanWindowing(packet,FFTLen,wlength,giType,numSamplesBeforeShortGI);
        tx = windowedPacket(aLen+(1:numPktSamples), :);
        % Overlap start of packet with end
        tx(1:bLen,:) = tx(1:bLen,:)+windowedPacket(end-bLen+1:end,:);
        % Overlap end of packet with start
        tx(end-aLen+1:end,:) = tx(end-aLen+1:end,:)+windowedPacket(1:aLen,:);

        % Add trailing zeros to allow for channel delay
        padTx = [tx; zeros(100, 1)];

        tgbbChannel.SNRdB = snr(isnr);
        rx = tgbbChannel(padTx);
        % reset(tgbbChannel); % Reset channel for different realization

        % Pass the waveform through AWGN channel
        % rx = awgnChannel(rx);

        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rx, cfgNHT.ChannelBandwidth);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end

        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)), :);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf, cfgNHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx, fs, -coarseFreqOff);

        % Extract the non-HT fields and determine fine packet offset
        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)), :);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields, cfgNHT.ChannelBandwidth);
        
        % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;

        % If packet detected outside the range of expected delays from the
        % channel modeling; packet error
        if pktOffset>100
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end

        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)), :);
        fineFreqOff = wlanFineCFOEstimate(lltf, cfgNHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx, fs, -fineFreqOff);

        % Extract L-LTF samples from the waveform, demodulate and perform
        % channel estimation
        lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)), :);
        lltfDemod = wlanLLTFDemodulate_custom(lltf, cfgNHT, 1);
        chanEst = wlanLLTFChannelEstimate_custom(lltfDemod, cfgNHT);
        
        % Get estimate of the noise power from L-LTF
        nVar = helperNoiseEstimate(lltfDemod);
        
        % Recover L-SIG field bits
        [rxLSIGBits, failCheck, eqLSIGSym] = wlanLSIGRecover_custom(rx(pktOffset + (ind.LSIG(1):ind.LSIG(2)), :), ...
            chanEst, nVar, cfgNHT.ChannelBandwidth);
        % if failCheck % Skip L-STF length of samples and continue searching
        %     disp('** L-SIG check fail **');
        % else
        %     disp('L-SIG check pass');
        % end
        % Obtain the transmitted LSIG bits
        txLSIGBits = wlanLSIGRecover_custom(lsig,ones(52,1),0,cfgNHT.ChannelBandwidth);

        % Extract non-HT Data samples from the waveform and recover the PSDU
        nhtdata = rx(pktOffset+(ind.NonHTData(1):ind.NonHTData(2)), :);
        
        % Recover data
        rxPSDU = wlanNonHTDataRecover(nhtdata, chanEst, nVar, cfgNHT);
        
        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr([txLSIGBits;inpPSDU], [rxLSIGBits;rxPSDU]));
        numPacketErrors = numPacketErrors+packetError;
        numPkt = numPkt+1;

    end

    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(isnr) = numPacketErrors/(numPkt-1);
    disp(['MCS ' num2str(cfgNHT.MCS) ','...
          ' SNR ' num2str(snr(isnr)) ','...
          ' Case ' num2str(whichCase) ','...
          ' Sps ' num2str(SamplesPerSymbol) ','...
          ' payloadLen ' num2str(cfgNHT.PSDULength) ','...
          ' completed after ' num2str(numPkt-1) ' packets,'...
          ' PER:' num2str(packetErrorRate(isnr))]);

end

varargout{1} = tgbbChannel;

end