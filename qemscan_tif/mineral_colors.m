% this is the old list
% min_col = ...  
%     [[255 255 255];
%      [234 225 125];
%      [255 200 0];
% 	 [255 192 128];
%      [255 128 0];
% 	 [128 64 0];
% 	 [205 133 63];
% 	 [139 0 0];
% 	 [192 64 0];
% 	 [174 0 0];
% 	 [64 0 0];
% 	 [51 153 102];
%      [154 205 50];
% 	 [128 128 0];
% 	 [99 150 78];
% 	 [128 128 128];
% 	 [192 255 192];
% 	 [143 188 139];
% 	 [192 255 255];
% 	 [0 255 255];
% 	 [0 0 255];
% 	 [255 255 0];
% 	 [255 128 255];
% 	 [255 0 0];
%      [192 0 0];
%      [128 128 255];
%      [138 43 226];
% 	 [128 0 128];
%      [0 255 0];
%      [255 128 128];
% 	 [255 0 153]]./255;
% 
% 
% mins = {'blank','Qtz','Ksp','Ab','An25','An75','Mscv','Bti_low','Bti_Mg','Bti_int',...
%     'Bti_Fe','Kao','Chl_Fe','Chl_Mg','Mg_Clay','Mg_Si','Ill_Smec','Ill_Smec_Fe','Cal',...
%     'Dol','Dol_Fe','Prt','Gyp_Anh','Rtl_Ilm','Ttn','Lmt','Cpx',...
%     'Fe_Amph','Epd_Zo','Apt','Zrc'};
% 
% minsN = {'blank','Qtz','Ksp','Ab','An25','An75','Mscv','Bti low','Bti Mg','Bti int',...
%     'Bti Fe','Kao','Chl Fe','Chl Mg','Mg Clay','Mg Si','Ill Smec','Ill Smec Fe','Cal',...
%     'Dol','Dol Fe','Prt','Gyp Anh','Rtl Ilm','Ttn','Lmt','Cpx',...
%     'Fe Amph','Epd Zo','Apt','Zrc'};


mins = {'blank','Qtz','Ksp','Ab','An25','An50','An75','Mscv','Bti_low','Bti_Mg','Bti_int',...
    'Bti_Fe','Kao','Chl_Fe','Chl_Mg','Mg_Clay','Mg_Si','Ill_Smec','Ill_Smec_Fe','Cal',...
    'Dol','Dol_Fe','Fe_ox_si','Prt','Gyp_Anh','Hlt','Rtl_Ilm','Ilm','Ttn','Lmt','Cpx',...
    'Fe_Amph','Epd_Zo','Apt','Trml','Zrc','AlOx'};

minsN = {'blank','Qtz','Ksp','Ab','An25','AN50','An75','Mscv','Bti low','Bti Mg','Bti int',...
    'Bti Fe','Kao','Chl Fe','Chl Mg','Mg Clay','Mg Si','Ill Smec','Ill Smec Fe','Cal',...
    'Dol','Dol Fe','Fe ox si','Prt','Gyp Anh','Hlt','Rtl Ilm','Ilm','Ttn','Lmt','Cpx',...
    'Fe Amph','Epd Zo','Apt','Trml','Zrc','AlOx'};

minNFull = ...
{'Background', 'Quartz', 'K Feldspar', 'Plagioclase Ab', 'Plagioclase An25',...
'Plagioclase An50', 'Plagioclase An75','Muscovite', 'Biotite (Low Fe & Mg)',...
'Biotite (Mg-rich)', 'Biotite (Intermediate)', 'Biotite (Fe-rich)', 'Kaolinite',...
'Fe Chlorite', 'Mg Chlorite', 'Mg Clays', 'Mg Silicate', 'Illite & illite-smectite',...
'Fe-Illite & illite-smectite','Calcite', 'Dolomite', 'Ferroan Dolomite',...
'Fe Oxide & siderite', 'Pyrite', 'Gypsum / Anhydrite', 'Halite', 'Rutile & Ilmenite',...
'Ilmenite', 'Titanite', 'Laumontite', 'Clinopyroxene', 'Fe Amphibole', 'Epidote / Zoisite',...
'Apatite', 'Tourmaline', 'Zircon', 'Aluminium Oxide'}

min_col = ...
[[255,255,255];
 [234,225,125];
 [255, 200, 0];
 [255,192,128];
 [255,128,0];
 [192,64,0];
 [128,64,0];
 [205,133,63];
 [139,0,0];
 [220,64,0];
 [174, 0,0];
 [125,0,0];
 [51,153,102];
 [154,205,50];
 [128,128,0];
 [99,150,78];
 [0, 100, 0];
 [192,255,192];
 [143,188,139];
 [192,255,255];
 [0,255,255];
 [0,0,255];
 [0,0,0];
 [255,255,0];
 [255,128,255];
 [255,192,255];
 [255,0,0];
 [192,0,0];
 [220,0,0];
 [128,128,255];
 [138,43,226];
 [128,0,128];
 [0,255,0];
 [255,128,128];
 [189,183,107];
 [255,0,153];
 [128,128,128]];
]