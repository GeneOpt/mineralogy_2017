%%% Glacier 1 12_17



minerals = cellstr(['Qtz';'Plg';'Mcr';'Sph';'Chl';'Act';'Bti';'Hbl';'Epi']);
allMins = 1:length(minerals);

a = [0,25,35,45,60,70,72,85,100];
b = [0,20,30,4055,65,68,80,90,95,100];
c = [0,5,15,30,40,48,55,70,75,95,100];
d = [0,30,40,50,60,70,85,95,100];
e = [0,8,20,30,40,50,56,100];
f = [0,6,23,39,60,82,90,100];

aa  = [1,2,1,2,1,1,1,2];
bb  = [2,2,7,7,1,1,2,3,7,2];
cc  = [3,2,1,3,2,7,2,1,3,2];
dd  = [3,2,3,2,3,3,1,7];
ee  = [2,1,2,2,1,3,2];
ff = [7,1,1,2,1,2,2];


lines = 6;

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
    
    if k == 6
        ind = f;
        indMin = ff;
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

totLength = sum(aveLengthMtx)
allLength = sum(totLength)
modalComp17 = sum(aveLengthMtx)/allLength*100

minCount = sum(1 - (aveLengthMtx == 0))
aveLengthTot = (totLength./minCount)
alt17 = aveLengthTot

    in = isnan(aveLengthTot);
    for p = 1:length(aveLengthTot)
        if in(p) == 1;
            alt17(p) = 0;
        end
    end





