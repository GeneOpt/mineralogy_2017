%% Qemscan for sediment part A

%% bulk mineralogy - mass percent 
[dat_BMMP_S,var_BMMP_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',2,'A1:AR21');
dat_BMMP_S =dat_BMMP_S(:,3:end)
%% bulk mineralogy - area percent 
[dat_BMAP_S,var_BMAP_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',3,'A1:AR21');

%% Detailed PSD (Area %)
[dat_PSD_S,var_PSD_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',4,'A1:H141');

%% Min by Size - Mass % in Sample
[dat_SMPS_S,var_SMPS_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',5,'A1:AT81');
%% Min by Size - Area % in Sample
[dat_SAPS_S,var_SAPS_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',6,'A1:AT81');
%% Min by Size - Mass % in Group
[dat_SMPG_S,var_SMPG_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',7,'A1:AT81');
%% Min by Size - Area % in Group
[dat_SAPG_S,var_SAPG_S] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_suspended_seds.xlsx',8,'A1:AT81');


%% dont forget about glacier 1 powder

% 
% %% Here we have the Fragmented rock
% %% bulk mineralogy - mass percent...OKAY, I HAVE TAKEN A BREAK HERE, BECAUSE I THINK THAT THE RESULTS NEED TO BE REASSESED. SKIPPING THE WHOLE ROCK
% [dat1,var1] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\Fragmented Rocks Merged.xlsx',2,'A1:AS33')
% %% bulk mineralogy - area percent
% [dat1,var1] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8152-AUG16 ARQS Datasheet SFU Crompton - Fragmented Rock.xlsx',3,'A1:AR33')
% %% Detailed PSD (Area %)
% [dat1,var1] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan\MI8152-AUG16 ARQS Datasheet SFU Crompton - Fragmented Rock.xlsx',4,'A1:I113')


%% Qemscan for fragmented rock 
%% bulk mineralogy - mass percent
[dat_BMMP_F,var_BMMP_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',2,'A1:AT27');
dat_BMMP_F = dat_BMMP_F(:,3:end)
%% bulk mineralogy - area percent
[dat_BMAP_F,var_BMAP_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',3,'A1:AR27');
%% Detailed PSD (Area %)
[dat_PSD_F,var_PSD_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',4,'A1:I183');
%% Min by Size - Mass % in Sample
[dat_SMPS_F,var_SMPS_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',5,'A1:AU105');
%% Min by Size - Area % in Sample
[dat_SAPS_F,var_SAPS_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',6,'A1:AU105');
%% Min by Size - Mass % in Group
[dat_SMPG_F,var_SMPG_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',7,'A1:AU105');
%% Min by Size - Area % in Group
[dat_SAPG_F,var_SAPG_F] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',8,'A1:AU105');


%% sample properties
[dat_ss,var_ss] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',1,'A1:G27');
[dat_sed,var_sed] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\sample_specs.xlsx',2,'A1:F21');


%% Mineral digest for sediment 

[datMD,varMD] = xlsread('D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\QS_fragRock.xlsx',9,'A1:BI31');

