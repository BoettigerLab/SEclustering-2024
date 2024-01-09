
GenomeFolder = 'U:\';
NAS02_Vol4 = 'N:\';
saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\';
mm10coords = readtable([NAS02_Vol4,'Derek\','20210722_L10_SE_MultiplexedCells\Elist_coords_mm10.txt']);
c= 19;
damID_chr = ReadTableFile([GenomeFolder,'GenomeData\ByPublication\vanSteensel2022\GSE181693_damid_average_replicates.tsv']);
disp(damID_chr(1,:))

nE = height(mm10coords);
lad_score = nan(nE,1);
for e=1:nE
    isChr = strcmp(mm10coords.Var1{e},damID_chr.seqnames);
    damID_c = damID_chr(isChr,:);
    inRange = damID_c.start >= mm10coords.Var2(e)-10e3 & damID_c.xEnd <= mm10coords.Var3(e)+10e3;
    lad_score(e) = nanmean(cellfun(@str2double,damID_c.PT_0h(inRange)));
end

save([saveFolder,'LADscore.mat'],'lad_score');