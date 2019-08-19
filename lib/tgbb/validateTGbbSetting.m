function isCorrect = validateTGbbSetting(cfgFormat,tgbbChannel)
%validateTGbbSetting checks whether the simulation setup is correct or not
%   ISCORRECT = validateTGbbSetting(CFGFORMAT,TGBBCHANNEL) checks whether
%   the simulation setup is correct or not. 
%
%

%   Logs:
%   16 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Created

assert(tgbbChannel.isIncludeDelay,'Delay should be considered!')
assert(tgbbChannel.isApplyNoise,'Noise must be injected before the Rx FE!')


end
