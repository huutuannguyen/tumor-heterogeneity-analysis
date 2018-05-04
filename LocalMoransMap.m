function [Group]=LocalMoransMapV3(Data,max_distance,iterations, mode) %mode= 3 for HER2/cell 2 for HER2/CEP17 and 1 for IF 
global_std=std(Data(:,3));
if mode==2 
    Low=1.8;
    High=2.2;
   global_threshold_center=2;
     k=0;    
end
if mode==3
    Low=4;
    High=6;
    global_threshold_center=5;
    k=0; 
end 

if mode==1
       global_threshold_center=mean(Data(:,3));
    k=1;
    distanceToCenter=k*global_std;

    Low=global_threshold_center- distanceToCenter;
    High=global_threshold_center+ distanceToCenter;
end 
if mode==4
   global_threshold_center=0.2512;

    k=1;
    distanceToCenter=0;
    Low=global_threshold_center- distanceToCenter;
    High=global_threshold_center+ distanceToCenter;
    
end 
Group=zeros(size(Data,1),3);

num_obs=size(Data,1);
Distances=pdist(Data(:,1:2),'euclidean');
W=zeros(size(Distances));
for i=1:numel(W)
    if Distances(i)<=max_distance
        W(i)=1;
    end
end
W=squareform(W);

z=NaN(num_obs,1);
Wz=NaN(num_obs,1);

LowThreshold=(Low-global_threshold_center)/global_std;
HighThreshold=(High-global_threshold_center)/global_std;

    for i=1:size(Data,1)
    values=W(:,i).*Data(:,3);
    values(values==0)=[];
    standardized_values=(values-global_threshold_center)/global_std;
    z(i)=(Data(i,3)-global_threshold_center)/global_std;
    Wz(i)=sum(standardized_values)/numel(values);
    if z(i)<LowThreshold && Wz(i)<0,Group(i,1)=0;Group(i,2)=0; Group(i,3)=1;end  
    if z(i)<LowThreshold && Wz(i)>=0,
        Group(i,1)=0;Group(i,2)=1; Group(i,3)=1;end 
    if z(i)>=HighThreshold && Wz(i)<0, Group(i,1)=1;Group(i,2)=0; Group(i,3)=1;end 
    if z(i)>=HighThreshold && Wz(i)>=0, Group(i,1)=1;Group(i,2)=1; Group(i,3)=1;end 

    end

end


