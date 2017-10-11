function [rc] = ind2rc(nRow,ind)
rt = ind/nRow;
cl = ceil(rt);
rw = (rt-(cl-1))*nRow;
rc = [rw,cl];