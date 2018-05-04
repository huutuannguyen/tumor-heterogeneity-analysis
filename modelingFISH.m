close all;

totalNuc=20000;
Nmax=50;
m=0:totalNuc;
R=7.01;
to=4;
AMcellrateArray=[];
NFISHArray=[];
NFISHMatrix=[];
AMcellrateMatrix=[];
figure
  hold on
for noSimulation=0:10
for N=1:Nmax
    
Xarray=[];
parray=[];
    for i=0:totalNuc
       x= random('unif',-R,R);
       if (x>(R-to/2))
           Vt=2*pi*R^3/3-pi*R^2*(x-to/2)+pi*(x-to/2)^3/3;
       end;
       if (x<(-R+to/2))
           Vt=2*pi*R^3/3-pi*R^2*(abs(x)-to/2)+pi*(abs(x)-to/2)^3/3;
       end;
       if and((x>=(-R+to/2)), (x<=(R-to/2)))
       Vt=pi*to*(3*R*R-to*to/4-3*x*x)/3;
       end;
       p=Vt/(pi*4*R^3/3);
    cdf = random('unif',0,1)                %cdf cumulative densitty function
    xinv=binoinv(cdf,N,p);
    Xarray=[Xarray;xinv];
    parray=[parray;p];
    end;
%     figure
%    hist(Xarray);
    NFISH=mean(Xarray);
NFISHArray=[NFISHArray;NFISH];
if NFISH<=6 HetCell=Xarray(Xarray>6);
end;
if NFISH>6 HetCell=Xarray(Xarray<4);   
end;
   if isempty(HetCell) Hetcellrate=0;
   else     Hetcellrate=numel(HetCell);
   end
    AMcellrateArray=[AMcellrateArray;Hetcellrate];
    end;
    AMcellrateMatrix=horzcat(AMcellrateMatrix,AMcellrateArray);
    NFISHMatrix=horzcat(NFISHMatrix,NFISHArray);
    plot (NFISHArray,AMcellrateArray/totalNuc);
    AMcellrateArray=[];
    NFISHArray=[];
end;