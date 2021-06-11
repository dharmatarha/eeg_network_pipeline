function []=get_net_param(dirname)

cd(dirname)
a=dir;
a=regexpi({a.name},'^[atlag].*[xls]', 'match');
subjlist=[a{:}];

out_data=zeros(length(subjlist),36);

for x =1:length(subjlist)
    
    temp_data=xlsread(subjlist{x});
    
    subjname=regexpi(subjlist{x},'\d{2,}', 'match');
    
    MST=get_MST(temp_data); %PLI ertekek MST-be rendezve, szimmetrikus matrixban
    degree=degrees_und(MST)/14;
    leaf=sum(degrees_und(full(MST))==1)/14;
    %[lambda,efficiency,ecc,radius,diameter] = charpath(distance_bin(MST));
    [~,~,~,~,diameter] = charpath(distance_bin(MST));
    diameter=diameter/14;
    MST_bin=double(full(MST)~=0);
    BC=(betweenness_bin(MST_bin)/14)/14;
    Th=(leaf*14)/(2*14*max(BC));
    
    out_data(x,:)=[str2num(subjname{:}), leaf, diameter, Th, max(degree), max(BC), degree, BC];
    
    
    
end

cd('..')

deg_label=[];
BC_label=[];

for k=1:15
    deg_label=[deg_label {['degree Ch. ' num2str(k)]}];
    BC_label=[BC_label {['BC Ch. ' num2str(k)]}];
end

header = [{'subject', 'leaf', 'diameter', 'Th', 'max degree', 'max BC'}, deg_label, BC_label];

xlswrite([dirname '.xls'],header);
xlswrite([dirname '.xls'],out_data,['A2:AJ' num2str(length(subjlist)+1)]);



