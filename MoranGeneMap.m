function [ output_args ] = MoranGeneMap(Data, sample_name)              % Data is the merged data imported after Total process, output is local indication of spatial heteorogeneity and mean values of IF and FISH scores
 mkdir('I:\experiment\IHC FISH\Joint project with Daniel\FINAL MERGED-AND-FILTERED TABLES TO ANALYZE',char(sample_name));
 folder_name=['I:\experiment\IHC FISH\Joint project with Daniel\FINAL MERGED-AND-FILTERED TABLES TO ANALYZE\', char(sample_name)];

rows=Data.Y<-1100; Data(rows,:)=[];
h=HetMapOfficial(Data);

LHcells=NaN;
HLcells=NaN;
HHcells=NaN;
LLcells=NaN;

Xmax=max(Data.X);
Ymax=min(Data.Y);

OEvar2='OEHER2';
Data.OEHER2=Data.HER2Expression;
radius=100;

vars={'X','Y',OEvar2}; 

Data.MoranOEHER2=LocalMoransMapV3(Data{:,vars},radius,1,1);

hOEH=figure('Name','OverexpressionHER2');
gscatter(Data.X,Data.Y,Data.MoranOEHER2,'kgcbr','.*+xo',12);
xlim([0,Xmax]);
ylim([Ymax,0]);
saveas(hOEH,[folder_name '\OE.tif']);


OEvar4='OEratio';
Data.OEratio=(Data.HER2Expression-Data.backgroundHER2)./Data.CytokeratinExpression;
vars={'X','Y',OEvar4}; 
Data.MoranOERatio=LocalMoransMapV3(Data{:,vars},radius,1,4);

hOER=figure('Name','OverexpressionRatio');
gscatter(Data.X,Data.Y,Data.MoranOERatio,'kgcbr','.*+xo',12);
xlim([0,Xmax]);
ylim([Ymax,0]);

saveas(hOER,[folder_name '\Overexpression.tif']);
LHcells=Data(and(and(Data.OEratio(:,1)==0,(Data.MoranOERatio(:,2)==1)), (Data.MoranOERatio(:,3)==1)),{'Tilename','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression','OEratio'}); 
LHcellNo=size(LHcells,1);
HLcells=Data(and(and(Data.MoranOERatio(:,1)==1,(Data.MoranOERatio(:,2)==0)), (Data.MoranOERatio(:,3)==1)),{'Tilename','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression','OEratio'});
HLcellNo=size(HLcells,1);
HHcells=Data(and(and(Data.MoranOERatio(:,1)==1,(Data.MoranOERatio(:,2)==1)), (Data.MoranOERatio(:,3)==1)),{'Tilename','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression','OEratio'});
HHcellNo=size(HHcells,1);
LLcells=Data(and(and(Data.MoranOERatio(:,1)==0,(Data.MoranOERatio(:,2)==0)), (Data.MoranOERatio(:,3)==1)),{'Tilename','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression','OEratio'});
LLcellNo=size(LLcells,1);
notDefCell=Data((Data.MoranOERatio(:,3)==0),{'Tilename','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression','OEratio'});
notDefCellNo=size(notDefCell,1);

writetable(HHcells,[folder_name '\OE_HHtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
writetable(LLcells,[folder_name '\OE_LLtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
writetable(HLcells,[folder_name '\OE_HLtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
writetable(LHcells,[folder_name '\OE_LHtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);

saveas(h,[folder_name '\HetMap.tif']);
OEvar1='AMratio';
Data.AMratio=Data.HER2Loci./Data.CEP17Loci;
vars={'X','Y',OEvar1}; 
Data.MoranAMratio=LocalMoransMapV3(Data{:,vars},radius,1,2);
hAMR=figure('Name','AMratio');
gscatter(Data.X,Data.Y,Data.MoranAMratio,'kgcbr','.*+xo',12);

xlim([0,Xmax]);
ylim([Ymax,0]);


saveas(hAMR,[folder_name '\AMR.tif']);
OEvar3='NFISH';     %Different than HER2 Loci
Data.Nnuclei=floor(Data.NucleusSize/mean(Data.NucleusSize));            %Number of nuclei per cluster of nuclei
Data.Nnuclei(Data.Nnuclei<1)=1;
Data.NFISH=Data.HER2Loci./Data.Nnuclei;
vars={'X','Y',OEvar3}; 
Data.MoranAMNFISH=LocalMoransMapV2(Data{:,vars},radius,1,3);

%h4 = figure;
hNFISH=figure('Name','NFISH');
gscatter(Data.X,Data.Y,Data.MoranAMNFISH,'kgcbr','.*+xo',12);
xlim([0,Xmax]);
ylim([Ymax,0]);
saveas(hNFISH,[folder_name '\Amplification.tif']);

Data.HNucleuspositionXFISH=Data.NucleuspositionX*2;
Data.NucleuspositionYFISH=Data.NucleuspositionY*2;
TotalCellNo=size(Data,1);

LHcells=Data(and(and(Data.MoranAMNFISH(:,1)==0,(Data.MoranAMNFISH(:,2)==1)), (Data.MoranAMNFISH(:,3)==1)),{'Tilename','HNucleuspositionXFISH','NucleuspositionYFISH','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression'}); 
LHcellNo=size(LHcells,1);
HLcells=Data(and(and(Data.MoranAMNFISH(:,1)==1,(Data.MoranAMNFISH(:,2)==0)), (Data.MoranAMNFISH(:,3)==1)),{'Tilename','HNucleuspositionXFISH','NucleuspositionYFISH','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression'});
HLcellNo=size(HLcells,1);
HHcells=Data(and(and(Data.MoranAMNFISH(:,1)==1,(Data.MoranAMNFISH(:,2)==1)), (Data.MoranAMNFISH(:,3)==1)),{'Tilename','HNucleuspositionXFISH','NucleuspositionYFISH','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression'});
HHcellNo=size(HHcells,1);
LLcells=Data(and(and(Data.MoranAMNFISH(:,1)==0,(Data.MoranAMNFISH(:,2)==0)), (Data.MoranAMNFISH(:,3)==1)),{'Tilename','HNucleuspositionXFISH','NucleuspositionYFISH','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression'});
LLcellNo=size(LLcells,1);
notDefCell=Data((Data.MoranAMNFISH(:,3)==0),{'Tilename','HNucleuspositionXFISH','NucleuspositionYFISH','HER2Loci','CEP17Loci','HER2Expression','CytokeratinExpression'});
notDefCellNo=size(notDefCell,1);

meanNFISH=mean(Data.NFISH);
meanCEP17=mean(Data.CEP17Loci);
HERCEPRatio=meanNFISH/meanCEP17;
stdNFISH=std(Data.NFISH);
Data.HCR1=Data.OEratio;
meanHCR1=mean(Data.HCR1);
STDHCR1=std(Data.HCR1);
HER2status=NaN;
Data.HERCEPRATIO=Data.NFISH./Data.CEP17Loci;
meanHERCEPRATIOCELL=mean(Data.HERCEPRATIO);
STDHERCEPRATIOCELL=std(Data.HERCEPRATIO);

STDNFISH=std(Data.NFISH);
if HERCEPRatio>2 HER2status=2;
end

if meanNFISH>6    HER2status=2;
end

if and(and(meanNFISH<6, meanNFISH>4), (HERCEPRatio<=2)) 
    HER2status=1;
end    

if and(meanNFISH<4, HERCEPRatio<=2) HER2status=0;
end 

if HER2status<2 infiltratingHETrate=HLcellNo/TotalCellNo;
    clusterHETrate=HHcellNo/TotalCellNo;
end

 if   (HER2status==2) infiltratingHETrate=LHcellNo/TotalCellNo;
clusterHETrate=LLcellNo/TotalCellNo;

 end
writetable(HHcells,[folder_name '\AM_HHtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
writetable(LLcells,[folder_name '\AM_LLtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
writetable(HLcells,[folder_name '\AM_HLtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
writetable(LHcells,[folder_name '\AM_LHtable'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);

writetable(Data,[folder_name '\Data'],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
output_args=[meanHCR1;meanNFISH; HERCEPRatio; HHcellNo; HLcellNo; LHcellNo;  LLcellNo; infiltratingHETrate; clusterHETrate;STDNFISH; meanCEP17; mean(Data.HER2Loci); meanHERCEPRATIOCELL; STDHERCEPRATIOCELL;STDHCR1;stdNFISH];

end


