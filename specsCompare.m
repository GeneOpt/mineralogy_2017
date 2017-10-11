function [mlc, mfc, sz, tp] = specsCompare(data)

sz = 10;

%MSNS
if data == 1
    mlc = [0 0 0];
    mfc = [0.6 0.6 0.6];
    tp = '^';
end

%MSS
if data == 2
    mlc = [0 0 0];
    mfc = [0.6 0.6 0.6];
    tp = 'o';
end

%MXS
if data == 3
    mlc = [0 0 0];
    mfc = [1 0 0];
    tp = 'o';
end

%MXNS
if data == 4
    mlc = [0 0 0];
    mfc = [1 0 0];
    tp = '^';
end

%PG
if data == 5
    mlc = [0 0 0];
    mfc = [1 0.6 0.6];
    tp = '^';
end

%SG
if data == 6
    mlc = [0 0 0];
    mfc = [0.3 0.3 1];
    tp = 's';
end

%BH
if data == 7
    mlc = [0 0 0];
    mfc = [0.3 1 0.3];
    tp = 's';
end