function oriIndex = C6togene(varargin)
% e.g. C6togene('./data/cancerGeneList.mat', 'MYC_UP.V1_DN', 'MYC_UP.V1_UP')
% (1st argument is patients' gene list; 2nd:end is target oncogenic gene sets)
% Find all the relevant genes from the target oncogenic gene sets (2:end argument),
% and return the corresponding index of patients data (1st argument). 

load('../matdata/Gene2Oncogenic.mat'); % Gene2ONCO stringArray
load(varargin{1}); %'cancerGeneList.mat'; %tumorGA

% check unexisted oncogenic genes in tumorGA
unexist = Gene2ONCO(find(ismember(Gene2ONCO(:,1), tumorGA) ==0), 1);

%disp("Number of input arguments: " + nargin)
%celldisp(varargin)

lastUnionIndex = [];
% find the union gene of all the given C6
for i = 2:nargin
     tmpIndex = find(sum(strcmp(varargin{i}, Gene2ONCO(:,3:end)),2));
     lastUnionIndex = union(tmpIndex, lastUnionIndex);
end

oriIndex = zeros(size(lastUnionIndex));
[isMem, oriIndex] = ismember(Gene2ONCO(lastUnionIndex,1), tumorGA);
% remove unmatched/unexisted oncogenic genes in tumorGA
oriIndex = nonzeros(oriIndex);

end