function [ Table1,Table2 ] = BuildingTableTHBPTI( Presortnum, Filenum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long
na{1}=sprintf('Michael%d.mat',Presortnum);
load(na{1});
Freqpresort1=cell(size(MutantsMat1,1),2);
avglib=length(AAseqMatrix);
min=0/avglib;
for inc=1:size(Freqpresort1,1)
    Freqpresort1{inc,1}=MutantsMat1{inc,1};
    Freqpresort1{inc,2}=double(cell2mat(MutantsMat1(inc,2)))/avglib;
    if Freqpresort1{inc,2}<min 
       Freqpresort1{inc,2}=min;
    end
end

Freqpresort2=cell(size(MutantsMat2,1),2);
for inc=1:size(Freqpresort2,1)
    Freqpresort2{inc,1}=MutantsMat2{inc,1};
    Freqpresort2{inc,2}=double(cell2mat(MutantsMat2(inc,2)))/avglib;
    if Freqpresort2{inc,2}<min
       Freqpresort2{inc,1}=min;
    end
end


%%frequencies sorted and enrich
na{1}=sprintf('Michael%d.mat',Filenum);
load(na{1});
avglib=length(AAseqMatrix);
Frequency1=zeros(size(MutantsMat1,1),1);
Frequency2=zeros(size(MutantsMat2,1),1);
Enrich1=zeros(size(MutantsMat1,1),1);
Enrich2=zeros(size(MutantsMat2,1),1);
Average1=zeros(size(MutantsMat1,1),1);
Average2=zeros(size(MutantsMat2,1),1);
Mutants1=cell(size(MutantsMat1,1),1);
Mutants2=cell(size(MutantsMat2,1),1);


pos1=1;
for inc=1:size(Frequency1,1)
    F1=double(cell2mat(MutantsMat1(inc,2)))/avglib;
    loc=find(ismember(Freqpresort1(:,1),MutantsMat1(inc,1)));
    if F1<(0/avglib)
        if isempty(loc) || double(cell2mat(Freqpresort1(loc,2)))==min
            continue;
        else
            Frequency1(pos1)=0/avglib;
            Average1(pos1)=double(cell2mat(MutantsMat1(inc,2)));
            Mutants1(pos1)=MutantsMat1(inc,1);
        end
    else
        Frequency1(pos1)=F1;
        Average1(pos1)=double(cell2mat(MutantsMat1(inc,2)));
        Mutants1(pos1)=MutantsMat1(inc,1);
    end
    if isempty(loc)
        Enrich1(pos1)=Frequency1(pos1)/min;
    else
        Enrich1(pos1)=Frequency1(pos1)/double(cell2mat(Freqpresort1(loc,2)));
    end
    pos1=pos1+1;
end

pos2=1;
for inc=1:size(Frequency2,1)
    F2=double(cell2mat(MutantsMat2(inc,2)))/avglib;
    loc=find(strcmp(Freqpresort2(:,1),MutantsMat2(inc,1)));
    if F2<(0/avglib)
        if isempty(loc) || double(cell2mat(Freqpresort2(loc,2)))==min
            continue;
        else
            Frequency2(pos2)=(0/avglib);
            Average2(pos2)=double(cell2mat(MutantsMat2(inc,2)));
            Mutants2(pos2)=MutantsMat2(inc,1);
        end
    else
        Frequency2(pos2)=F2;
        Average2(pos2)=double(cell2mat(MutantsMat2(inc,2)));
        Mutants2(pos2)=MutantsMat2(inc,1);
    end
    if isempty(loc)
        Enrich2(pos2)=Frequency2(pos2)/min;
    else
        Enrich2(pos2)=Frequency2(pos2)/double(cell2mat(Freqpresort2(loc,2)));
    end
    pos2=pos2+1;
end

%% Table
Table1=table(Mutants1,Average1,Frequency1,Enrich1);
Table2=table(Mutants2,Average2,Frequency2,Enrich2);
na{1}=sprintf('EnrichTable%d.mat',Filenum);
save(na{1});


end

