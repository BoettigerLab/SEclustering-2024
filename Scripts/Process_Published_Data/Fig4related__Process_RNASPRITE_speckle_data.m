
chrTableHyb = readtable([saveFolder,'chrTableHyb_mm10.csv']) ; % mm10  

fname = 'U:\GenomeData\ByPublication\Guttman2021\GSE186264_DMSO.RDSPRITE.clusters.txt';
% mm10
fid = fopen(fname);
ftxt = textscan(fid,'%s','Delimiter','\n');  % '\r\n' % any combination of \r and \n carraige return and newline
fclose(fid);
ftxt

hasMalat = contains(ftxt{1},'Malat1'); sum(hasMalat)
hasU1 = contains(ftxt{1},'U1.snRNA'); sum(hasU1)
hasU2 = contains(ftxt{1},'U2.snRNA'); sum(hasU2)
hasU4 = contains(ftxt{1},'U4.snRNA'); sum(hasU4)
hasU6 = contains(ftxt{1},'U6.snRNA'); sum(hasU6)
hasSpeck = hasMalat | hasU1  | hasU2 | hasU4 | hasU6;

speckClusters = ftxt{1}(hasSpeck);

hasDNA = contains(speckClusters,'DPM');
speckleDNA = speckClusters(hasDNA);
nS = length(speckleDNA);
dnaMap = cell(nS,1);

res =0.5e6;  
spriteMat = zeros(24*0.25e9/res,1); % each chromosome gets a 250 Mb domain in a linear array.  
for s=1:nS
    disp(['processing cluster s=',num2str(s),' of ',num2str(nS)]);
    contacts = strsplit(speckleDNA{s},'\t');
    if length(contacts) < 100 % toss super large clusers
        isDNAcontacts = contains(contacts,'DPM');
        sum(isDNAcontacts)/length(isDNAcontacts);
        dnaLocs = cellfun(@(x) x(8:end),contacts(isDNAcontacts),'UniformOutput',false);
        [chrs,chrStart,chrStops] = ParseLocusName(dnaLocs);
        chrNums = str2double(regexprep(chrs,'chr',''));
        chrNums(isnan(chrNums)) = 23;
        % record the precise mapped coordinates
        dnaMap{s} = chrNums*1e9+chrStart;
        
        % bin the data now. 
        % this way each cluster can only contribute 1 read per bin, and not
        % create artificial peaks
        genomeIdx = floor( (chrNums*0.25e9 + chrStart)/res)+1; % genome ID, in Mb.  
        spriteMat(genomeIdx) = spriteMat(genomeIdx)+1;
       
    end
end

%%
allSpeckMap = cat(1,dnaMap{:});  % allMalatMap


%% chrC

c=11; 
xs =  floor(c*0.25e9/res)+1 : floor((c+1)*0.25e9/res)+1 ; % the separate floors are critical 
b = spriteMat(xs);
xe = find(b>0,1,'last')+1;
x = (0:xe-1)*res;
figure(2); clf;  bar(x,b(1:xe));
chrC_table = chrTableHyb(chrTableHyb.chrNum==c,:);
cs = chrC_table.start;
figure(2); hold on; plot(cs,zeros(length(cs),1),'r.')
contTable = chrTableHyb([91,257,291,296,333],:)
chrC_table = contTable(contTable.chrNum==c,:);
cs = chrC_table.start;
figure(2); hold on; plot(cs,zeros(length(cs),1),'r+')

%% chrC

stp = .25e6;
chrCMap = allSpeckMap(allSpeckMap>c*1e9 & allSpeckMap<(c+1)*1e9);
[b,x] = hist(chrCMap-c*1e9,0:stp:250e6);
figure(1); clf; bar(x,b,1);
chrC_table = chrTableHyb(chrTableHyb.chrNum==c,:);
cs = chrC_table.start;
figure(1); hold on; plot(cs,zeros(length(cs),1),'r.')
csr = stp*round(cs/stp);
figure(1); hold on; plot(csr,zeros(length(cs),1),'go');

%% all chr
figure(10); clf;
stp = 0.5e6; % bin size
speckSE = cell(23,1);
for c=1:23 % c=17
    try
    % map (bin) speckles to genome at desired resolution 
    chrCMap = allSpeckMap(allSpeckMap>c*1e9 & allSpeckMap<(c+1)*1e9);
    [b,x] = hist(chrCMap-c*1e9,0:stp:250e6);
    % map SEs to genome at same resolution 
    chrC_table = chrTableHyb(chrTableHyb.chrNum==c,:);
    cs = chrC_table.start;
    csr = stp*round(cs/stp);
    nSE = length(csr);
    speckSE{c} = nan(nSE,1);
    for s=1:nSE
        speckSE{c}(s) = b(x==csr(s));
    end
    k=c;
    xe  = x(find(b>0,1,'last'));
    if c==23; k=k-2; end
    figure(10); subplot(7,3,k); cla; bar(x,b,1);
    hold on; plot(cs,zeros(length(cs),1),'r.')
    xlim([1,xe]); title(c); ylim([0,150])
    catch
    end
end
speckleSE=  cat(1,speckSE{:});
chrTableHyb(speckleSE > 30,:)
%%
saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\'
% save([saveFolder,'speckleSE.mat'],'speckleSE');
% save([saveFolder,'chrTableHyb.mat'],'chrTableHyb');

%%