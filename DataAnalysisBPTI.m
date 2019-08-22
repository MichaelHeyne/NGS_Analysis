
for n = 1:1
    WT = 'RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRT';
    primer = 'TAGC';
    ThreshOld = 1; %Singles - 1, doubles - 2
    kk = 2;
    kkk = 1;
    filename = sprintf('2Mix_doubles_Si.xlsx',n);% output
    na{1} = sprintf('2Mix.combined.extendedFrags.fastq',n); % input (fastq file)
    Data = fastqread(na{1});
    DataprepareFast;
    AAseqMatrix = sortcell(AAseqMatrix,1);
    [MutantsMat1,MutantsMat2]=MutantsCount3FastBPTI(WT,AAseqMatrix);
    na{1} = sprintf('Michael%d.mat',kk);
    save(na{1});
    if ThreshOld == 1
        xlswrite(filename,MutantsMat1);
    elseif ThreshOld == 2
        xlswrite(filename,MutantsMat2);
    else
        disp('Threshold exceeds range');
    end
    fprintf('%s was created.\n',filename);
    clear all;
end 
     %building table
%     na{1}=sprintf('Michael%d.mat',kk);
%     load(na{1});
%     [Table1,Table2]=BuildingTableTHBPTI(2,2);
%     na{1}=sprintf('EnrichtableBPTI%d.mat',kk);
%     save(na{1});
%     clear Table1 Table2
