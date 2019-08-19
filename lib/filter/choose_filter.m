function filterStruct = choose_filter(fName,varargin)
% choose_filter obtains the numerator and denumerator coefficients of filters
%   
%   FILTERSTRUCT = choose_filter(FILTERNAME,NTHORDER,FTYPE,FCNORM) is 
%   a container/struct that describe the NTHORDER FILTERNAME with the cut-
%   off frequency of FCNORM. FCNORM is defined as f3dB/(fs/2), where fs is the 
%   symbol rate (downsampled) or sample rate (upsampled).
%
%   FILTERNAME must be one of:
%       'butter', 'fetx', 'ferx', 'fetxtamas', 'ferxtamas'
%
%   If the butterowrth filter is used, FTYPE must be one of:
%        'low', 'bandpass', 'high', 'stop'
%

%   Logs:
%   16 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Created

% Parsing input and validate
setFilterName = {'butter','fetx','ferx','fetxtamas','ferxtamas'};
setFilterType = {'low','bandpass','high','stop'};
p = inputParser;

addRequired(p,'fName',@(x) any(validatestring(x,setFilterName)));
addOptional(p,'fcNorm',0.5, @(x) validateattributes(x,{'numeric'},{'>',0,'<=',1}))
addOptional(p,'nthOrder',uint8(1),@(x) x == round(x))
addOptional(p,'fType','low',@(x) any(validatestring(x,setFilterType)))
addOptional(p,'verbose',false,@(x) validateattributes(x,{'logical'},{'scalar'}))

parse(p,fName,varargin{:});

fName = p.Results.fName;
if strcmp(fName,'fetx')
    fcNorm = [2.6e5/5e8,2.34e8/5e8];
elseif strcmp(fName,'ferx')
    fcNorm = [4.8e4/5e8,2.58e8/5e8];
else 
    fcNorm = p.Results.fcNorm;
end
nthOrder = p.Results.nthOrder;
fType = p.Results.fType;
verbose = p.Results.verbose;

% Obtain filte's coeff.
fName = lower(fName);

switch fName
    case 'butter'
        % Butterworth
        [num, den] = butter(nthOrder,fcNorm,fType);
    case 'fetx'
        % Tx's front-end
        [num, den] = getFrontEndTx();
     case 'ferx'
        % Rx's front-end
        [num, den] = getFrontEndRx();
    case 'fetxtamas'
        % Tx's front-end
        [num, den] = getFrontEndTxTamas();
     case 'ferxtamas'
        % Rx's front-end
        [num, den] = getFrontEndRxTamas();
    otherwise
        error('%s is not implemented yet!',fName)

end

% Generate output
filterStruct.name = fName;
filterStruct.order = nthOrder;
filterStruct.num = num;
filterStruct.den = den;
filterStruct.fcNorm = fcNorm;
filterStruct.grpDelay = grpdelay(num,den,1);

% Freq. response as a function of f, assuming
% f has been normalized by the corresponding freq. 
% sampling, i.e., f/(fs/2)
filterStruct.H = @(f) freqz(...
    num,den,2*pi*f);%.*phasez(num,den,2*pi*f);


% CIR
maxMemoryLength = 2^10; % maximum memory length
threshold = 1-1e-6;      % energy threshold above which impulse respone is neglected

if den == 1 % FIR
    filterStruct.h = num;
else
    x = zeros(1, maxMemoryLength+1);
    x(1) = 1;
    y = filter(filterStruct.num, filterStruct.den, x);
    E = cumsum(abs(y).^2)/sum(abs(y).^2);
    y(E > threshold) = [];
    % y = y/abs(sum(y)); % normalize to have unit gain at DC
    filterStruct.h = y;
end

% Draw
if verbose

    if strcmp(filterStruct.name,'fetx') | ...
        strcmp(filterStruct.name,'ferx') 

        f_bw = 5e8; % Reference bandwidth (Hz)

        f = linspace(0,0.5,2^14); % in order to match the plot in 11-18-1574-04-00bb-lc-frontend-models
        freqResp = filterStruct.H(f);

        figure
        subplot(2,1,1); 
        semilogx(f*f_bw*2, pow2db(abs(freqResp).^2),'DisplayName', filterStruct.name)
        hold on;
        aa = axis;
        for idx = 1:length(filterStruct.fcNorm)
            h = semilogx(f_bw*[filterStruct.fcNorm(idx) filterStruct.fcNorm(idx)],[aa(3) aa(4)]);
            hasbehavior(h, 'legend', false);   % line will not be in legend
        end
        xlabel('Freq. [Hz]')
        ylabel('Magnitude (dB)')
        legend('-DynamicLegend')
        axis([5e4 5e8 -60 10])

        subplot(2,1,2);
        semilogx(f*f_bw*2, (unwrap(angle(freqResp))), 'DisplayName', filterStruct.name)
        xlabel('Freq. [Hz]')
        ylabel('Phase (rad)')
        legend('-DynamicLegend')
        % axis([1e5 1e9 -360 360])
        xlim([5e4 5e8])

    elseif strcmp(filterStruct.name,'fetxtamas') | ...
        strcmp(filterStruct.name,'ferxtamas')

        f_bw = 20e6*50; % Reference bandwidth (Hz)

        f = linspace(0,0.5,2^14); % in order to match the plot in 11-18-1574-04-00bb-lc-frontend-models
        freqResp = filterStruct.H(f);

        figure
        subplot(2,1,1); 
        semilogx(f*f_bw*2, pow2db(abs(freqResp).^2),'DisplayName', filterStruct.name)
        hold on;
        aa = axis;
        for idx = 1:length(filterStruct.fcNorm)
            h = semilogx(f_bw*[filterStruct.fcNorm(idx) filterStruct.fcNorm(idx)],[aa(3) aa(4)]);
            hasbehavior(h, 'legend', false);   % line will not be in legend
        end
        xlabel('Freq. [Hz]')
        ylabel('Magnitude (dB)')
        legend('-DynamicLegend')
        axis([5e4 5e8 -60 10])

        subplot(2,1,2);
        semilogx(f*f_bw*2, (unwrap(angle(freqResp))), 'DisplayName', filterStruct.name)
        xlabel('Freq. [Hz]')
        ylabel('Phase (rad)')
        legend('-DynamicLegend')
        % axis([1e5 1e9 -360 360])
        xlim([5e4 5e8])

    else

        % f = linspace(-0.5,0.5); % 0.5 is due to w = 2*pi*f and I want to plot it from -pi to pi
        f = -0.5:0.01:0.5; % 0.5 is due to w = 2*pi*f and I want to plot it from -pi to pi
        freqResp = filterStruct.H(f);
        figure
        subplot(2,1,1), hold on, box on;
        plot(2*f, pow2db(abs(freqResp).^2),'DisplayName', filterStruct.name)
        aa = axis;
        for idx = 1:length(filterStruct.fcNorm)
            h = plot([filterStruct.fcNorm(idx) filterStruct.fcNorm(idx)],[aa(3) aa(4)],'k--');
            hasbehavior(h, 'legend', false);   % line will not be in legend
            h = plot([-filterStruct.fcNorm(idx) -filterStruct.fcNorm(idx)],[aa(3) aa(4)],'k--');
            hasbehavior(h, 'legend', false);   % line will not be in legend
        end
        xlabel('Normalized frequency f/(fs/2) [x \pi rad/sample]')
        ylabel('Magnitude (dB)')
        legend('-DynamicLegend')
        axis([-1 1 -60 10])
            
        subplot(2,1,2), hold on, box on
        % plot(2*f, rad2deg(unwrap(angle(freqResp))), 'DisplayName', filterStruct.name)
        plot(2*f, (unwrap(angle(freqResp))), 'DisplayName', filterStruct.name)
        xlabel('Normalized frequency f/(fs/2) [x \pi rad/sample]')
        ylabel('Phase (rad)')
        legend('-DynamicLegend')
        % axis([-1 1 -360 360])

    end

end

end

function [num,den] = getFrontEndTx() 
    % see: https://mentor.ieee.org/802.11/dcn/18/11-18-1574-05-00bb-lc-frontend-models.pptx

    f_bw = 5e8;     % Reference bandwidth [Hz]
 
    %% Highpass filter
    n_hi               = 2;         % Filter order
    f_c_hi             = 2.6e5 ;    % cut-off frequency [Hz]
    [z_hi, p_hi, k_hi] = butter(n_hi, f_c_hi/f_bw, 'high');
    [sos_hi, g_hi]     = zp2sos(z_hi, p_hi, k_hi); 
     
    %% Lowpass filter
    n_lo               = 8;             % Filter order
    f_c_lo             = 2.34e8 ;   % Cut-off frequency [Hz]
    [z_lo, p_lo, k_lo] = butter(n_lo, f_c_lo/f_bw);
    [sos_lo, g_lo]     = zp2sos(z_lo, p_lo, k_lo);
     
    %% Combined bandpass filter
    passband_gain      = -23.17;    % Passb. gain [dB]
    sos                = [sos_hi; sos_lo];
    g                  = g_hi*g_lo*10^(passband_gain/20);
    H                  = dfilt.df2sos(sos, g);

    [num,den] = sos2tf(sos,g);

end

function [num,den] = getFrontEndRx() 
    % see: https://mentor.ieee.org/802.11/dcn/18/11-18-1574-05-00bb-lc-frontend-models.pptx

    f_bw = 5e8 ;    % Reference bandwidth (Hz)
 
    %% Highpass filter
    n_hi                = 4;        % Filter order
    f_c_hi          = 4.8e4;    % Highpass cut-off frequency (Hz)
    [z_hi, p_hi, k_hi]  = butter(n_hi, f_c_hi/f_bw, 'high');
    [sos_hi, g_hi]      = zp2sos(z_hi, p_hi, k_hi); 
     
    %% Lowpass filter
    n_lo                = 4;        % Filter order
    f_c_lo          = 2.58e8;   % Lowpass cut-off frequency (Hz)    
    [z_lo, p_lo, k_lo]  = butter(n_lo, f_c_lo/f_bw);
    [sos_lo, g_lo]      = zp2sos(z_lo, p_lo, k_lo);
     
    %% Combined bandpass filter
    passband_gain       = 4.6;  % Passb. gain (dB)
    sos                 = [sos_hi; sos_lo];
    g               = g_hi*g_lo*10^(passband_gain/20);
    H               = dfilt.df2sos(sos, g);

    [num,den] = sos2tf(sos,g);
    
end

function [num,den] = getFrontEndTxTamas() 


    % From Tamas
    f_bw = 20e6*50/2; % Reference bandwidth (Hz)
    %% Highpass filter
    n_hi               = 2;         % Filter order
    f_c_hi             = 2.6e5 ;    % cut-off frequency [Hz]
    [z_hi, p_hi, k_hi] = butter(n_hi, f_c_hi/f_bw, 'high');
    [sos_hi, g_hi]     = zp2sos(z_hi, p_hi, k_hi); 
     
    %% Lowpass filter
    n_lo               = 8;             % Filter order
    f_c_lo             = 20e6 ;     % Cut-off frequency [Hz]
    [z_lo, p_lo, k_lo] = butter(n_lo, f_c_lo/f_bw);
    [sos_lo, g_lo]     = zp2sos(z_lo, p_lo, k_lo);
     
    %% Combined bandpass filter
    passband_gain      = -23.17;    % Passb. gain [dB]
    sos                = [sos_hi; sos_lo];
    g                  = g_hi*g_lo*10^(passband_gain/20);
    H                  = dfilt.df2sos(sos, g);

    [num,den] = sos2tf(sos,g);

end

function [num,den] = getFrontEndRxTamas() 

    % From Tamas
    f_bw = 20e6*50/2; % Reference bandwidth (Hz)
    %% Highpass filter
    n_hi                = 4;        % Filter order
    f_c_hi          = 30e4;     % Highpass cut-off frequency (Hz)
    [z_hi, p_hi, k_hi]  = butter(n_hi, f_c_hi/f_bw, 'high');
    [sos_hi, g_hi]      = zp2sos(z_hi, p_hi, k_hi); 

    %% Lowpass filter
    n_lo                = 4;        % Filter order
    f_c_lo          = 2e8;  % Lowpass cut-off frequency (Hz)    
    [z_lo, p_lo, k_lo]  = butter(n_lo, f_c_lo/f_bw);
    [sos_lo, g_lo]      = zp2sos(z_lo, p_lo, k_lo);

    %% Combined bandpass filter
    passband_gain       = 3;    % Passb. gain (dB)
    sos                 = [sos_hi; sos_lo];
    g               = g_hi*g_lo*10^(passband_gain/20);
    H               = dfilt.df2sos(sos, g);
    [b,a] = sos2tf(sos,g);

    [num,den] = sos2tf(sos,g);
    
end