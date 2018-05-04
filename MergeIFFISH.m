function [ output_args ] = MergeIFFISH( IF, FISH )                 %import IF and FISH data, then run total process, output is a merged table with all FISH an
ok=0;

Type={'Tissue biopsy','HER2+ cell line','HER2- cell line'};

sample_name=inputdlg('Enter SAMPLE name:','Input',1);
batch_name=inputdlg('Enter BATCH name:','Input',1);
choice = menu('Choose the type of sample',Type);
if choice==1
    IHC_score=inputdlg('Enter IHC score:','Input',1);
    FISH_score=inputdlg('Enter FISH score:','Input',1);
    CEP17_score=inputdlg('Enter CEP17 score:','Input',1);
 end
folder_name=uigetdir('I:\experiment\IHC FISH\Joint project with Daniel','Choose where to save the Merge results');
ntilesX=char(inputdlg('Enter tile number in X:','Input',1));

IF.Properties.RowNames=IF.Cellcode;
IF.Cellcode=[];
IF.Code=[];
IF.signalHER2=[];
IF.SignalCYTO=[];
IF=IF(FISH.Cellcode,:);
FISH.Imagename=[];

Merged=[IF FISH];
Merged.Tilename=Merged.Code;
ntilesX=floor(str2double(ntilesX));
Merged.AbsolutePosX=mod(Merged.Code,(ntilesX))*150*0.645+Merged.NucleuspositionX*0.645;
Merged.AbsolutePosY=-(floor(Merged.Code/(ntilesX))*150*0.645+Merged.NucleuspositionY*0.645);

writetable(Merged,fullfile(folder_name,'Merged'),'WriteRowNames',false,'WriteVariableNames',true,'FileType','text'); 

Data=Merged;


%Filtering
IHC_score=NaN;
FISH_score=NaN;
CEP17_score=NaN;
n_low=0.5; 
n_high=10; 

n1=size(Data,1);

moran_distance=100;
moran_alpha=0.05;
moran_iterations=1;


Data.Properties.VariableNames(4)=cellstr('HER2Expression');
Data.Properties.VariableNames(6)=cellstr('CytokeratinExpression');
Data.Properties.VariableNames(15)=cellstr('HER2Loci');
Data.Properties.VariableNames(14)=cellstr('CEP17Loci');
Data.Properties.VariableNames(16)=cellstr('NucleusSize');
Data.Properties.VariableNames(25)=cellstr('X');
Data.Properties.VariableNames(26)=cellstr('Y');

if choice==1
    
    CK_median=median(Data.CytokeratinExpression);
    CK_spread=1.4826*mad(Data.CytokeratinExpression,1);
    min=CK_median-n_low*CK_spread;
    max=CK_median+n_high*CK_spread;
    rows=Data.CytokeratinExpression<min; Data(rows,:)=[]; rows=Data.CytokeratinExpression>max; Data(rows,:)=[]; % Histogram filterging
end
rows=isnan(Data.HER2Loci); Data(rows,:)=[];
rows=Data.HER2Loci<1; Data(rows,:)=[];
rows=Data.HER2Loci>25; Data(rows,:)=[];
rows=Data.CEP17Loci<1; Data(rows,:)=[];
rows=Data.NucleusSize<460; Data(rows,:)=[];
rows=Data.greenContrast<0.75; Data(rows,:)=[];
%% Adding sample name
for i=1:size(Data,1)
    Data.Sample(i)=sample_name;
    Data.Type(i)=Type(choice);
    Data.Batch(i)=batch_name;
    Data.IHCscore(i)=IHC_score;
    Data.FISHscore(i)=FISH_score;
    Data.CEP17score(i)=CEP17_score;
end

%% Saving

filteredDATA=FilterOverlap(Data, ntilesX);
folder_name='I:\experiment\IHC FISH\Joint project with Daniel\Result from batch analysis';
writetable(filteredDATA,[folder_name '\' char(sample_name)],'FileType','text','WriteVariableNames',true,'WriteRowNames',false);
output_args=MoranGeneMap(filteredDATA, sample_name);
save([folder_name '\analysisResult.mat']);
msgbox('Operation Completed');


end

