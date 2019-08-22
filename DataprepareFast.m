
%%
p=1;
b=1;
for i=1:length(Data)
    pos=StartPosition(Data(i).Sequence,primer);
    if ~size(pos)
        continue;
    end
    AAseq=SeqTranslationFast(Data(i).Sequence,pos(1),length(WT));
    [bin]=Alignmentfast(WT,ThreshOld,AAseq);
    if bin
        if ~isempty(strfind(AAseq,'*'))
            AAStopCodon{b,1}=AAseq;
            b=b+1;
        else
            AAseqMatrix{p,1}=AAseq;
            p=p+1;
        end
    end
end