function oriIndex = hallmark2gene(varargin)
% e.g. hallmark2gene('./data/cancerGeneList.mat', 'ESTROGEN_RESPONSE_LATE', 'ESTROGEN_RESPONSE_EARLY')
% (1st argument is patients' gene list; 2nd:end is target hallmark gene sets)
% Find all the relevant genes from the target hallmark gene sets (2:end argument),
% and return the corresponding index of patients data (tumorGA, 1st argument). 

load('../matdata/Gene2CancerHallmarks.mat'); % Gene2CH stringArray
load(varargin{1}); %'cancerGeneList.mat'; % tumorGA

% check unexisted hallmark genes in tumorGA
unexist = Gene2CH(find(ismember(Gene2CH(:,1), tumorGA) ==0), 1);

%disp("Number of input arguments: " + nargin)
%celldisp(varargin)

lastUnionIndex = [];
% find the union gene of all the given hallmark
for i = 2:nargin
     tmpIndex = find(sum(strcmp(varargin{i}, Gene2CH(:,3:end)),2));
     lastUnionIndex = union(tmpIndex, lastUnionIndex);
end

oriIndex = zeros(size(lastUnionIndex));
[isMem, oriIndex] = ismember(Gene2CH(lastUnionIndex,1), tumorGA);
% remove unmatched/unexisted hallmark genes in tumorGA
oriIndex = nonzeros(oriIndex);

end