function [ MutantsMat1,MutantsMat2 ] = MutantsCount3FastBPTI( WT,AAseqMatrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
MutantsMat1{1}='WT';
MutantsMat1{1,2}=0;
MutantsMat2=cell(1);
pos1=1;
pos2=0;
Flageq=0;
flagWT=0;
%size(diffMatrix,1)
for i=1:size(AAseqMatrix,1)
    if i>1&&i<size(AAseqMatrix,1)
        if isempty(find(~(AAseqMatrix{i}==AAseqMatrix{i-1}),1))
            Flageq=1;
        else
            Flageq=0;
        end
    end
    if ~Flageq
        indices=find(~(WT==AAseqMatrix{i}));
        if isempty(indices)
            flagWT=1;
        else
            flagWT=0;
        end
        if ~flagWT
            for j=1:length(indices)
                W=WT(indices(j));
                M=AAseqMatrix{i}(indices(j));
                if j>1
                    name{j}=sprintf('_%c%d%c',W,indices(j),M);
                else
                    name{j}=sprintf('%c%d%c',W,indices(j),M);
                end
            end
        end
    end
    if flagWT
        MutantsMat1{1,2}=MutantsMat1{1,2}+1;
    elseif length(indices)==1
        if ~Flageq
            pos1=pos1+1;
            FullName{1}=strcat(name{1:length(indices)});
            MutantsMat1{pos1,1}=FullName{1};
            MutantsMat1{pos1,2}=1;
        else
            MutantsMat1{pos1,2}=MutantsMat1{pos1,2}+1;
        end
        
    elseif length(indices)==2
        if ~Flageq
            pos2=pos2+1;
            FullName{1}=strcat(name{1:length(indices)});
            MutantsMat2{pos2,1}=FullName{1};
            MutantsMat2{pos2,2}=1;
        else
            MutantsMat2{pos2,2}=MutantsMat2{pos2,2}+1;
        end
        
    end
end
end

