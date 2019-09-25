function listMatDir = getListMatDir(keywords)
% getListMatDir returns a list of directory of mat files
%
%   LISTMATDIR = getListMatDir(KEYWORDS) returns a list of 
%   directories of mat files according to KEYWORDS. KEYWORDS
%   must be a cell of strings. LISTMATDIR is a cell of strings.
%   The order in KEYWORDS matters.
%
%   KEYWORDS must be one of:
%       'all', 'optical', 'empty', 'enterprise-conference', 
%       'enterprise-office', 'hospital', 'industrial', 
%       'residential', 'individual', 'overall'

%   Logs:
%   21 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Created
%   30 Aug 19, ardimas (ardimasandipurwita@outlook.com):
%       Adding the 'empty' and 'enterprise' keywords 

allKeywords = {...
    'all',
    'optical',
    'empty',
    'enterprise',
    'enterprise-conference',
    'enterprise-office',
    'hospital',
    'industrial',
    'residential',
    'individual',
    'overall'
};

if ~all(ismember(keywords,allKeywords))
    disp('keywords must contain:')
    disp(allKeywords)
    error('keywords are not correct!')
end

% Get all dirs of all mat files
listMatFiles 

% Filters out listMatDir according to keywords
for idx = 1:length(keywords)
    idxFind = find(contains(listMatDir,keywords{idx},'IgnoreCase',true));
    listMatDir = listMatDir(idxFind);
end


end