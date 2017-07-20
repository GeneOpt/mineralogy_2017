function [newI] = keepOnly(Ir,Ig,Ib,mineral,black,repCol)

elr = (abs(Ir - mineral(1))<1e-2);
elg = (abs(Ig - mineral(2))<1e-2);
elb = (abs(Ib - mineral(3))<1e-2);

notM = (logical(1-(elr.*elg.*elb)));
Ir(notM)=1;
Ig(notM)=1;
Ib(notM)=1;

if black == 1
    M = (logical((elr.*elg.*elb)));
    Ir(M)=repCol(1);
    Ig(M)=repCol(2);
    Ib(M)=repCol(3);
end

newI(:,:,1) = Ir;
newI(:,:,2) = Ig;
newI(:,:,3) = Ib;
