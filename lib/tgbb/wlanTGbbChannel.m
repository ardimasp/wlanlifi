classdef wlanTGbbChannel < matlab.System
%wlanTGbbChannel filters input signal and adds noise based on TGbb specification
%   tgbb = wlanTGbbChannel creates an object that follows TGbb specifications [1]. 
%   Details of it are provided in the methodology document as found in [2]. Initial 
%   simulation results are also presented in [3].
%
%   Step method syntax:
%
%   Y = step(tgbb,x) filter input signal X according to the TGbb sepcification.s
%   
%   tgbb.visualize() visualizes the signals of interest in both time domain (amplitude)
%       and frequency domain (PSD).
%   
%   wlanTGbbChannel properties:
%
%   SymbolRate          - Input signal symbol rate (symbols/sec)
%   SamplesPerSymbol    - Samples per symbol for upconversion 
%   Offset              - Freq. offset for the upconversion
%   beta                - Rollof factor for the pulse shaping
%   span                - The number of spanned symbols
%   ClippingRatio       - The clipping ratio for the DC bias
%   useData             - Indicator whether to use a pre-calculated WOCIR
%   whichCase           - Index of predefined case
%   whichMethod         - Method that is used to compare FD, TD and SP methods
%   useButterworth      - Indicator whether to use a butterworth filter
%   fcNormButter        - The cutoff frequency for the butterworth filter
%   isIncludeDelay      - Include delay due to filters or not
%   methodDCandAGC      - A method used for DC bias and AGC
%   varNoiseFloor       - The variance of the noise floor (UNUSED)
%   isApplyNoise        - Injecting noise before the RX FE
%   isShotNoise         - Indicator to apply a signal-dependent shot noise
%   whichFE             - Type of AFEs used
%   
%   % References:
%   % [1] http://www.ieee802.org/11/Reports/tgbb_update.htm
%   % [2] https://mentor.ieee.org/802.11/dcn/19/11-19-0187-04-00bb-evaluation-methodology-for-phy-and-mac-proposals.docx
%   % [3] https://mentor.ieee.org/802.11/dcn/19/11-19-1224-01-00bb-simulation-results-for-802-11a-phy-in-lc.ppt
%
%   See also wlanTGnChannel, wlanTGacChannel, wlanTGahChannel.
%

%   TODO:
%       - The most important thing is to make sure that the steps follow the TGbb spec
%       - Adds an example in the documentation
%       - Define a local random generator seed
%       - Pulse shaping should be replaced by a DAC and an ADC models, especisally for high sample rate

%   Logs:
%   16 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Created



properties (Nontunable)
    %Symbol rate (symbols/sec)
    SymbolRate = 20e6;
    %Samples per symbol 
    %   This is used the continuous-time simulation
    %   '11-19-0187-04-00bb-evaluation-methodology-for-phy-and-mac-proposals.docx' 
    %   suggest that the sample rate is 1 Gsamp/s. However, using 
    %   the default value of symbol rate, it needs 50 samples per symbol.
    %   This will lengthen the simulation time. I tried 4, 6, 8, 10, 20 and 50.
    %   The optimal value (trade-off between performance and runtime) is 10.
    SamplesPerSymbol = 10;
    %Offset (Hz)
    %   This is the offset for the upconversion
    %   see 11-19-0187-04-00bb-evaluation-methodology-for-phy-and-mac-proposals
    Offset = 1.5e6;
    %Rolloff factor
    %   This is for the pulse shaping using the raised cosine
    beta = 0.5;
    %Span
    %   Again this is for the pulse shaping using the raised cosine
    span = 8;
    %Clipping ratio
    %   Clipping ratio is used to determine the DC bias, i.e., DC = clipping ratio*std(sTx)
    ClippingRatio = 3;
    %Flag whether using the pre-calculated WOCIR
    useData = false;
    %If yes, which case
    whichCase = 1;
    %which method: {'fd','td','sp'}
    whichMethod = 'fd';
    %useButterworth
    useButterworth = false;
    %fcNormButter
    fcNormButter = 0.95;
    %including delay from the filters or not
    isIncludeDelay = true;
    %methodDCandAGC: {'ave','max'}
    %   'ave'
    %       DC ->
    %           the DC signal is multiple of the average power of its input
    %       AGC -> 
    %           the average power of AGC output is the same as that of before VGA
    %    'max'
    %       DC ->
    %           the DC signal is first carried out such that it  fits within [-1,1] 
    %           depending the max of the absolute value of the input signal and then make it all positive.
    %       AGC -> 
    %           the scaling factor of the input signal is a constant divided by the max of the input signal
    methodDCandAGC = 'ave'
    % noise floor
    %   modeling the thermal and shot noise
    %   see https://mentor.ieee.org/802.11/dcn/18/11-18-1423-08-00bb-tgbb-simulation-scenarios.docx
    %   noise floor is -70dBm
    %   Notes: 
    %       Please double check. If the SNR is measured after the wocir, then the power signal can be less 
    %       than db2pow(-100) = 1e-10. Remember that the DC gain of the wocir can be -100 dB.
    % varNoiseFloor = db2pow(-70-30);
    varNoiseFloor = 0;
    % apply optical noise
    isApplyNoise = false;
    % Add shot noise
    isShotNoise = false;
    %which FE: {'tgbb','tamas'}
    whichFE = 'tgbb'
end

properties (Nontunable)
    %Packet length or the number of symbols in a packet
    lenPacket
    %Sampling freq. (samples/sec)
    SampleRate
    %CIR of the raised cosine
    cirRcos
    %Freq. center of the upconversion
    Fc
    %Time bins at downsampled
    tDown
    %Time bins at upsampled
    tUp
    %Freq bins at downsampled
    fDown
    %Freq bins at upsampled
    fUp
    %Tx AFE
    %   see 11-18-1574-04-00bb-lc-frontend-models
    fetx
    %Rx AFE
    %   see 11-18-1574-04-00bb-lc-frontend-models
    ferx
end

properties 
    %SNR
    SNRdB = 100
    %Signal input
    sInput
    %Signal after upsampling
    sUpsampled
    %Signal after upsampling+pulse shape
    sPulseShaped
    %Signal after upconversion
    sUpconverted
    %Signal after AFE Tx
    sTx
    %Signal after DC bias
    sPos
    %Signal after WOCIR
    sWOCIR
    %Signal after AWGN
    sNoise
    %Signal after AFE Rx
    sRx
    %Signal after AGC
    sAGC
    %Signal after downconversion
    sDownconverted
    %Signal after match filter
    sMatchFiltered
    %Signal after downsampling
    sDownsampled
    %Freq responses
    freqResp_cirRcos
    freqResp_TxAFE
    freqResp_RxAFE
    freqResp_WOCIR
    %Scaling factor
    sfPulseShaped = 1
    sfUpConverted = 1
    sfAGC = 1
    sfDownConverted = 1
    sfMatchedFiltered = 1
    sfDownSampled = 1
    % Others
    varArtificialNoise
end

methods
    function obj = wlanTGbbChannel(varargin) % Constructor
        setProperties(obj, nargin, varargin{:});
    end 

    function visualize(obj,varargin)
    %Visualize all time domains and PSDs of signals of interest.    

        % To do: modify this function such that we can opt for which signal to be plotted
        % Now, just plot all of them.

        setWhichPlot = {'td','psd'};
        setWhichSignal = {...
            'sInput',...
            'sUpsampled',...
            'sPulseShaped',...
            'sUpconverted',...
            'sTx',...
            'sPos',...
            'sWOCIR',...
            'sNoise',...
            'sRx',...
            'sAGC',...
            'sDownconverted',...
            'sMatchFiltered',...
            'sDownsampled',...
            'freqResp_cirRcos',...
            'freqResp_TxAFE',...
            'freqResp_RxAFE',...
            'freqResp_WOCIR'};

        figure
        title('Discrete-time Signals')
        subplot(2,1,1), hold on; box on;
        stem(obj.tDown,real(obj.sInput),'DisplayName','sInput'); 
        stem(obj.tDown,real(obj.sDownsampled),'DisplayName','sDownsampled'); 
        ylabel('Re[x]')
        xlabel('time [s]')
        xlim([0.7e-5,0.8e-5])
        legend('-DynamicLegend')
        subplot(2,1,2), hold on; box on;
        stem(obj.tDown,imag(obj.sInput),'DisplayName','sInput'); 
        stem(obj.tDown,imag(obj.sDownsampled),'DisplayName','sDownsampled'); 
        ylabel('Imag[x]')
        xlabel('time [s]')
        legend('-DynamicLegend')
        xlim([0.7e-5,0.8e-5])


        figure
        title('Continuous-time Signals')
        subplot(2,1,1), hold on; box on;
        stem(obj.tDown,real(obj.sInput)/(max(real(obj.sInput))),'DisplayName','sInput'); 
        plot(obj.tUp,real(obj.sPulseShaped)/(max(real(obj.sPulseShaped))),'DisplayName','sPulseShaped'); 
        plot(obj.tUp,real(obj.sUpconverted)/(max(real(obj.sUpconverted))),'DisplayName','sUpconverted'); 
        plot(obj.tUp,real(obj.sTx)/(max(real(obj.sTx))),'DisplayName','sTx'); 
        plot(obj.tUp,real(obj.sPos)/(max(real(obj.sPos))),'DisplayName','sPos'); 
        plot(obj.tUp,real(obj.sWOCIR)/(max(real(obj.sWOCIR))),'DisplayName','sWOCIR'); 
        plot(obj.tUp,real(obj.sNoise)/(max(real(obj.sNoise))),'DisplayName','sNoise'); 
        plot(obj.tUp,real(obj.sRx)/(max(real(obj.sRx))),'DisplayName','sRx'); 
        plot(obj.tUp,real(obj.sAGC)/(max(real(obj.sAGC))),'DisplayName','sAGC'); 
        plot(obj.tUp,real(obj.sDownconverted)/(max(real(obj.sDownconverted))),'DisplayName','sDownconverted'); 
        plot(obj.tUp,real(obj.sMatchFiltered)/(max(real(obj.sMatchFiltered))),'DisplayName','sMatchFiltered'); 
        stem(obj.tDown,real(obj.sDownsampled)/(max(real(obj.sDownsampled))),'DisplayName','sDownsampled'); 
        ylabel('Re[x]')
        xlabel('time [s]')
        legend('-DynamicLegend')
        xlim([0.7e-5,0.8e-5])
        subplot(2,1,2), hold on; box on;
        stem(obj.tDown,imag(obj.sInput)/(max(imag(obj.sInput))),'DisplayName','sInput'); 
        plot(obj.tUp,imag(obj.sPulseShaped)/(max(imag(obj.sPulseShaped))),'DisplayName','sPulseShaped'); 
        plot(obj.tUp,imag(obj.sUpconverted)/(max(imag(obj.sUpconverted))),'DisplayName','sUpconverted'); 
        plot(obj.tUp,imag(obj.sTx)/(max(imag(obj.sTx))),'DisplayName','sTx'); 
        plot(obj.tUp,imag(obj.sPos)/(max(imag(obj.sPos))),'DisplayName','sPos'); 
        plot(obj.tUp,imag(obj.sWOCIR)/(max(imag(obj.sWOCIR))),'DisplayName','sWOCIR'); 
        plot(obj.tUp,imag(obj.sNoise)/(max(imag(obj.sNoise))),'DisplayName','sNoise'); 
        plot(obj.tUp,imag(obj.sRx)/(max(imag(obj.sRx))),'DisplayName','sRx'); 
        plot(obj.tUp,imag(obj.sAGC)/(max(imag(obj.sAGC))),'DisplayName','sAGC'); 
        plot(obj.tUp,imag(obj.sDownconverted)/(max(imag(obj.sDownconverted))),'DisplayName','sDownconverted'); 
        plot(obj.tUp,imag(obj.sMatchFiltered)/(max(imag(obj.sMatchFiltered))),'DisplayName','sMatchFiltered'); 
        stem(obj.tDown,imag(obj.sDownsampled)/(max(imag(obj.sDownsampled))),'DisplayName','sDownsampled'); 
        ylabel('Imag[x]')
        xlabel('time [s]')
        legend('-DynamicLegend')
        xlim([0.7e-5,0.8e-5])

        figure
        title('PSDs')
        subplot(11,1,1);
        x = obj.sInput;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SymbolRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SymbolRate,'centered');
        plot(f/1e6,pow2db(abs(H).^2))
        xlabel('f [MHz]')
        ylabel('sInput')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,2);
        x = obj.sPulseShaped;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        magRespFilt = abs(obj.freqResp_cirRcos).^2;
        magRespFilt = magRespFilt*max(abs(H).^2)/max(magRespFilt);
        plot(f/1e6,pow2db(abs(H).^2),obj.fUp/1e6,pow2db(magRespFilt)); 
        xlabel('f [MHz]')
        ylabel('sPulseShaped')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,3);
        x = obj.sUpconverted;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        plot(f/1e6,pow2db(abs(H).^2))
        xlabel('f [MHz]')
        ylabel('sUpconverted')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,4);
        x = obj.sTx;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        magRespFilt = abs(obj.freqResp_TxAFE).^2;
        magRespFilt = magRespFilt*max(abs(H).^2)/max(magRespFilt);
        plot(f/1e6,pow2db(abs(H).^2),obj.fUp/1e6,pow2db(magRespFilt));
        xlabel('f [MHz]')
        ylabel('sTx')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,5);
        x = obj.sWOCIR;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        magRespFilt = abs(obj.freqResp_WOCIR).^2;
        magRespFilt = magRespFilt*max(abs(H).^2)/max(magRespFilt);
        plot(f/1e6,pow2db(abs(H).^2),obj.fUp/1e6,pow2db(magRespFilt));
        xlabel('f [MHz]')
        ylabel('sWOCIR')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,6);
        x = obj.sNoise;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SymbolRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SymbolRate,'centered');
        plot(f/1e6,pow2db(abs(H).^2))
        xlabel('f [MHz]')
        ylabel('sNoise')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,7);
        x = obj.sRx;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        magRespFilt = abs(obj.freqResp_RxAFE).^2;
        magRespFilt = magRespFilt*max(abs(H).^2)/max(magRespFilt);
        plot(f/1e6,pow2db(abs(H).^2),obj.fUp/1e6,pow2db(magRespFilt));
        xlabel('f [MHz]')
        ylabel('sRx')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,8);
        x = obj.sAGC;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        plot(f/1e6,pow2db(abs(H).^2))
        xlabel('f [MHz]')
        ylabel('sAGC')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,9);
        x = obj.sDownconverted;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        plot(f/1e6,pow2db(abs(H).^2))
        xlabel('f [MHz]')
        ylabel('sDownconverted')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,10);
        x = obj.sMatchFiltered;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SampleRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SampleRate,'centered');
        magRespFilt = abs(obj.freqResp_cirRcos).^2;
        magRespFilt = magRespFilt*max(abs(H).^2)/max(magRespFilt);
        plot(f/1e6,pow2db(abs(H).^2),obj.fUp/1e6,pow2db(magRespFilt));
        xlabel('f [MHz]')
        ylabel('sMatchFiltered')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);
        subplot(11,1,11);
        x = obj.sDownsampled;
        % [H,f] = periodogram(x,hamming(floor(length(x)/4)),length(x),obj.SymbolRate,'centered');
        [H,f] = pwelch(x,hamming(floor(length(x)/4)),[],[],obj.SymbolRate,'centered');
        plot(f/1e6,pow2db(abs(H).^2))
        xlabel('f [MHz]')
        ylabel('sDownsampled')
        xlim([-obj.SymbolRate,obj.SymbolRate]/1e6);

    end

end

methods(Access = protected)
    function setupImpl(obj,varargin)

        % Get the packet length
        obj.lenPacket = length(varargin{1});

        % Sampling rate
        obj.SampleRate = obj.SymbolRate*obj.SamplesPerSymbol;

        % CIR for the raised cosine filter
        obj.cirRcos = rcosdesign(obj.beta,obj.span,obj.SamplesPerSymbol);

        % Freq. center for the upconversion
        obj.Fc = obj.SymbolRate/2 + obj.Offset;
        
        % Get the time and freq. bins
        [obj.fDown,obj.tDown] = freq_time(obj.lenPacket,obj.SymbolRate);
        [obj.fUp,obj.tUp] = freq_time(obj.lenPacket*obj.SamplesPerSymbol,obj.SampleRate);

        if strcmp(obj.whichFE,'tgbb')
            % Tx AFE
            obj.fetx = choose_filter('fetx');

            % Rx AFE
            obj.ferx = choose_filter('ferx');
        elseif strcmp(obj.whichFE,'tamas')
            
            % Tx AFE
            obj.fetx = choose_filter('fetxtamas');

            % Rx AFE
            obj.ferx = choose_filter('ferxtamas');
        end

        %% Magnitude Responses
        %cirRcos
        obj.freqResp_cirRcos = getfreqResp(obj.cirRcos,1,obj.fUp,obj.SampleRate,'isIncludeDelay',obj.isIncludeDelay);
        
        if obj.isIncludeDelay
            %Tx AFE
            obj.freqResp_TxAFE = (freqz(obj.fetx.num,obj.fetx.den,obj.fUp,1e9));
            %Rx AFE
            obj.freqResp_RxAFE = (freqz(obj.ferx.num,obj.ferx.den,obj.fUp,1e9));
        else
            % without delay : ignore the phase response
            %Tx AFE
            obj.freqResp_TxAFE = abs(freqz(obj.fetx.num,obj.fetx.den,obj.fUp,1e9));
            %Rx AFE
            obj.freqResp_RxAFE = abs(freqz(obj.ferx.num,obj.ferx.den,obj.fUp,1e9));
        end

        %WOCIR
        if obj.useData
            if strcmp(num2str(obj.whichCase),'S1_D1')
                % This is how Tamas used the loaded data
                load('S1_D1.mat')
                fs = 1e9;
                [h,f_data] = freqz(averun2,1,2001,'whole',fs);
                freqResp_data = (h);
            elseif strcmp(num2str(obj.whichCase),'S1_D2')
                % This is how Tamas used the loaded data
                load('S1_D2.mat')
                fs = 1e9;
                [h,f_data] = freqz(averun2,1,2001,'whole',fs);
                freqResp_data = (h);
            elseif strcmp(num2str(obj.whichCase),'S3_D1')
                % This is how Tamas used the loaded data
                load('S3_D1.mat')
                fs = 1e9;
                [h,f_data] = freqz(averun2,1,2001,'whole',fs);
                freqResp_data = (h);
            elseif strcmp(num2str(obj.whichCase),'S3_D2')
                % This is how Tamas used the loaded data
                load('S3_D2.mat')
                fs = 1e9;
                [h,f_data] = freqz(averun2,1,2001,'whole',fs);
                freqResp_data = (h);
            elseif strcmp(num2str(obj.whichCase),'industrial_d7')
                % This is how Tamas used the loaded data
                load('industrial_d7.mat')
                fs = 1e9;
                [h,f_data] = freqz(averun2,1,2001,'whole',fs);
                freqResp_data = (h);
            end
            if obj.isIncludeDelay % include the phase response, hence complex signal
                fitObjR=fit(f_data(:),real(freqResp_data(:)),'smoothingspline'); % real
                fitObjI=fit(f_data(:),imag(freqResp_data(:)),'smoothingspline'); % imag
                obj.freqResp_WOCIR = arrayfun(@(f) feval(fitObjR,abs(f))+1i*feval(fitObjI,abs(f)),obj.fUp);
            else
                freqResp_data = abs(freqResp_data);
                fitObj=fit(f_data(:),freqResp_data(:),'smoothingspline');
                obj.freqResp_WOCIR = arrayfun(@(f) feval(fitObj,abs(f)),obj.fUp);
            end
        else
            if obj.useButterworth
                butterFilt = choose_filter('butter','fcnorm',obj.fcNormButter,'nthOrder',4);
                obj.freqResp_WOCIR = getfreqResp(butterFilt.num,butterFilt.den,obj.fUp,obj.SampleRate,'isIncludeDelay',obj.isIncludeDelay);
            else
                obj.freqResp_WOCIR = ones(size(obj.fUp));
            end
        end

    end

    function varargout = stepImpl(obj,input)

        % I will repeatedly use a temporary variable x
        x = input;
        obj.sInput = x;

        % Upsampling 
        x = upsample(x,obj.SamplesPerSymbol);
        obj.sUpsampled = x;

        % Pulse shaping
        x = filtering(x,obj.freqResp_cirRcos);
        [x,obj.sfPulseShaped] = normalize(x,input); 
        obj.sPulseShaped = x;

        % Upconversion
        x = upConversion(x,obj.Fc,obj.tUp);
        [x,obj.sfUpConverted] = normalize(x,obj.sPulseShaped); 
        obj.sUpconverted = x;

        % Tx AFE
        x = filteringReal(x,obj.freqResp_TxAFE);
        obj.sTx = x;

        % DC bias
        if strcmp(obj.methodDCandAGC,'ave')
            x = x+obj.ClippingRatio*sqrt(var(x));
        elseif strcmp(obj.methodDCandAGC,'max')
            x = x/max(abs(x));
            x = x+1;
        end

        % Clip
        x(x<0) = 0; 

        obj.sPos = x;

        % Wireless optical CIR
        x = filteringReal(x,obj.freqResp_WOCIR);
        obj.sWOCIR = x;

        % Inject noise
        if obj.isApplyNoise
            varNoise = mean(abs(x).^2)/db2pow(obj.SNRdB);
            % UNUSED
            % First, check the issue in the note about noise floor in the 
            % properties section
            % obj.varArtificialNoise = varNoise-obj.varNoiseFloor;
            % At the meantime, use the following.
            obj.varArtificialNoise = varNoise;
            x = x + sqrt(varNoise)*randn(size(x));
        end
        if obj.isShotNoise
            x(x<0) = 0;
            x = x+sqrt(2*1.6021766208e-19*x*obj.SampleRate/2).*randn(size(x));
        end
        obj.sNoise = x;

        % Rx AFE
        x = filteringReal(x,obj.freqResp_RxAFE);
        obj.sRx = x;

        % AGC
        if strcmp(obj.methodDCandAGC,'ave')
            [x,obj.sfAGC] = normalize(x,obj.sUpconverted);
        elseif strcmp(obj.methodDCandAGC,'max')
            obj.sfAGC = 3/max(x);
            x = x*obj.sfAGC;
        end
        obj.sAGC = x;

        % Downconversion
        x = downConversion(x,obj.Fc,obj.tUp);
        obj.sDownconverted = x;

        % Anti-aliasing
        % Optimal matched filter
        x = filtering(x,conj(obj.freqResp_cirRcos));
        obj.sMatchFiltered = x;
        
        % Downsampling
        x = downsample(x,obj.SamplesPerSymbol);
        [x,obj.sfDownSampled] = normalize(x,obj.sDownconverted);
        obj.sDownsampled = x;

        varargout{1} = obj.sDownsampled;


    end
end

end

function freqResp = getfreqResp(num,den,f,fs,varargin)

    p = inputParser;
    addOptional(p,'isIncludeDelay',false,@(x) validateattributes(x,{'logical'},{'scalar'}))
    parse(p,varargin{:});
    isIncludeDelay = p.Results.isIncludeDelay;

    if isIncludeDelay
        freqResp = (freqz(num,den,f,fs));
    else
        freqResp = abs(freqz(num,den,f,fs)); % ignore the phase response
    end
end

function y = filtering(x,freqResp)
    %Filters the input signal x according to a frequency response
    % Doing it this way is more flexible in terms of setting whether 
    % we want to include the phase distortion (delay) or not.
    y = ifft( fft(x).* ifftshift( freqResp ) );
end

function y = filteringReal(x,freqResp)
    %The output is real. Hence, the symmetric argument.
    y = ifft( fft(x).* ifftshift( freqResp ), 'symmetric' );
end

function [y,scalingFactor] = normalize(x,xRef)
    %Normalize so that the energy of x is equal the energy of xRef
    scalingFactor = sqrt(mean(abs(xRef).^2)/mean(abs(x).^2));
    y = x*scalingFactor;

end

function y = upConversion(x,fc,t)
    y = real(x).*cos(2*pi*fc*t)+imag(x).*sin(2*pi*fc*t);
    y = real(y); % make sure that the output is real
end

function y = downConversion(x,fc,t)
    y = x.*cos(2*pi*fc*t)+1i*x.*sin(2*pi*fc*t);
end



