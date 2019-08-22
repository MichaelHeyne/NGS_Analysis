function [ bin ] = Alignmentfast( WT,ThreshOld,AAseq )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if length(AAseq)~=length(WT)
    bin=0;
else
    diff=WT==AAseq;
    k=sum(diff);
    k=length(WT)-k;
    if k>ThreshOld
        bin=0;
    else
        k= strfind(AAseq,'*');
        if ~isempty(k)
            bin=1;
        else
            bin=1;
        end
    end
end
end

