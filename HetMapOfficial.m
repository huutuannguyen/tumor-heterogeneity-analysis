function [output_args] = HetMapOfficial(Data) %2D plan of IF and FISH in the same image, convert to micrometer

Data.HER2Loci(Data.HER2Loci>30)=30;


Data.HER2CYTORatio=(Data.HER2Expression-Data.backgroundHER2)./Data.CytokeratinExpression;
Data.AMratio=Data.HER2Loci./Data.CEP17Loci;
FISH_min=min(Data.AMratio);
FISH_min=0;
FISH_max=max(Data.AMratio);
FISH_max=8;

FISH_range=FISH_max-FISH_min;

IF_min=0;
IF_max=max(Data.HER2CYTORatio);
IF_max=1;
IF_range=IF_max-IF_min;

xmax=max(Data.X);
ymax=min(Data.Y);

hAMR=figure('Name','coIFFISH');
hold on
for i=1:size(Data,1)
    if (Data.HER2Loci>0)
        if (Data.AMratio(i)>FISH_max) Data.AMratio(i)=FISH_max;
        end   
        if (Data.HER2CYTORatio(i)>IF_max) Data.HER2CYTORatio(i)=IF_max;
        end
         if (Data.HER2CYTORatio(i)<IF_min) Data.HER2CYTORatio(i)=IF_min;
        end
        Overexpression=(Data.HER2CYTORatio(i)-IF_min)/IF_range;
        Amplification=(Data.AMratio(i)-FISH_min)/FISH_range;
       plot(Data.X(i),Data.Y(i),'o','MarkerFaceColor',[Amplification 0 0],'MarkerEdgeColor',[Overexpression Overexpression 0],'LineWidth',1,'MarkerSize',10); 

    end
end
axis equal
xlim([0,xmax]);
ylim([ymax,0]);
output_args=hAMR;
end