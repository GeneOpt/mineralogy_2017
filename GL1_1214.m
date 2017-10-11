
%%% Glacier 1 12_14



minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi']);
allMins = 1:length(minerals);

a = [0,12,20,25,45,54,63,78,91,100];
b = [0,8,15,32,40,45,50,52,65,80,100];
c = [0,10,15,25,30,35,45,55,68,85,100];
d = [0,35,50,60,68,90,92,100];
e = [0,3,12,18,25,40,50,60,74,78,90,100];

aa  = [1,1,1,1,6,1,2,2,2];
bb  = [1,2,2,7,1,1,4,2,2,2];
cc  = [2,7,5,2,2,1,2,2,2,2];
dd  = [2,2,7,2,1,1,1];
ee  = [2,5,1,1,2,6,7,2,4,2,2];


lines = 5;

aveLengthMtx = zeros(lines,length(allMins));
numLengthMtx = zeros(lines,length(allMins));


for k = 1:lines
    

    if k == 1
        ind = a;
        indMin = aa;
    end
    if k == 2
        ind = b;
        indMin = bb;
    end
    if k == 3
        ind = c;
        indMin = cc;
    end
    if k == 4
        ind = d;
        indMin = dd;
    end
    if k == 5
        ind = e;
        indMin = ee;
    end

    totLength = zeros(1,length(allMins));
    numLength = zeros(1,length(allMins));

    for i = 1:(length(ind)-1)
        lgth = ind(i+1)-ind(i);
        j = indMin(i);
        totLength(j) = totLength(j)+lgth;
        numLength(j) = numLength(j)+1;
    end

    aveLength = (totLength./numLength)*0.05;
    in = isnan(aveLength);
    for p = 1:length(allMins)
        if in(p) == 1
            aveLength(p) = 0;
        end
    end
    
   
    aveLengthMtx(k,:) = aveLength;
    numLengthMtx(k,:) = numLength;
    
end

aveLengthMtx
totLength = sum(aveLengthMtx)
allLength = sum(totLength);
modalComp14 = sum(aveLengthMtx)/allLength*100;

minCount = sum(1 - (aveLengthMtx == 0));
aveLengthTot = (totLength./minCount)
alt14 = aveLengthTot;

    in = isnan(aveLengthTot);
    for p = 1:length(aveLengthTot)
        if in(p) == 1;
            alt14(p) = 0;
        end
    end

 


