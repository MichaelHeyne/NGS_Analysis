function [ AAseq ] = SeqTranslationFast( sequence,pos,WTsize)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
k=length(sequence)-WTsize*3+1;
if pos>k
    AAseq(1)='A';    
else
  code=sequence(pos:((pos+WTsize*3)-1));
  AAseq=nt2aa(code,'ACGTOnly', false);
end
    
    
end

