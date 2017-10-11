filename = ['D:\Field_data\2013\Summer\Geochemistry\2013_clayTypes.xlsx']

ExcelData = cellstr(['B';'D';'F';'H';'J';'L';'N';'B';'D'])
wpArr = zeros(1,11)'

for h = 1:length(ExcelData)
v1 = cellstr('5')
v2 = cellstr('15')
FC = ExcelData(h)
FC1 = strcat(FC{1},v1{1})
FC2 = strcat(FC{1},v2{1})
FC12 = strcat(FC1,:,FC2)
wp = xlsread(filename,1,FC12)
wpArr = [wpArr,wp]
end
wpArr = wpArr(:,2:end)