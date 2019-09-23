function cir = getCIRFromFreqRespTest(freqResp,varargin)
    %Obtain the time-domain CIR from a frequency response
    % This function is useful especially if the given frequency response is 
    % generated with sampling frequency that is not the same as the simulated one.

    p = inputParser;
    addOptional(p,'isReal',true,@(x) validateattributes(x,{'logical'},{'scalar'}))
    parse(p,varargin{:});
    isReal = p.Results.isReal;

    if isReal
        cir = filteringReal(ones(size(freqResp)),freqResp);
    else
        cir = filtering(ones(size(freqResp)),freqResp);
    end

    % Choose the significant ones
    threshold = 1-1e-6;      % energy threshold above which impulse respone is neglected
    E = cumsum(abs(cir).^2)/sum(abs(cir).^2);
    cir = cir(E < threshold);

end

function y = filteringReal(x,freqResp)
    %The output is real. Hence, the symmetric argument.
    y = ifft( fft(x).* ifftshift( freqResp ), 'symmetric' );
end