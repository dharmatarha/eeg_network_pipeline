load test
MST=get_MST(test); %PLI ertekek MST-be rendezve, szimmetrikus matrixban
degree=degrees_und(MST)/14;
leaf=sum(degrees_und(full(MST))==1)/14;
[lambda,efficiency,ecc,radius,diameter] = charpath(distance_bin(MST));
diameter=diameter/14;
MST_bin=double(full(MST)~=0);
BC=(betweenness_bin(MST_bin)/14)/14;
Th=(leaf*14)/(2*14*max(BC));