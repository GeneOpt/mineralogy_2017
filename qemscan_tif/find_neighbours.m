clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\revisedColors\rock_frags\MI8152-AUG16\1 - 25 um\';
[nms] = dir([folder '\*.TIF'])
matNm = {nms.name}

% for H = 1:length(matNm)
for H = 1
    fname = [folder matNm{H}];
    
    nm = matNm{H};
    for h = 1:length(nm)
    
        if strcmp(nm(h),'.')
            img = nm(1:h-1);
        end
        
    end
    
    [I,cmap] = imread(fname);

    Irgb = ind2rgb(I,cmap);
    figure
    imshow(Irgb)
    return
    %% grid search
    %%% here you cruise through each pixel of interest, and see what minerlas
    %%% are touching the other minerals. 

    lenMin = length(mins);
    numI = zeros(lenMin,lenMin);
    fullI = zeros(lenMin,lenMin);
    num_sg = zeros(1,lenMin-1);
    num_MinAll = zeros(1,lenMin-1);

    for l = 2:lenMin

        Irgb_loop = Irgb;
        % mineral = 'Ill_Smec'
        % minRGB = min_col(find(strcmp(mins,mineral)),:)
        minRGB = min_col(l,:); %grab the mineral that you want
        elr = (abs(Ir - minRGB(1))<1e-3); %find all the pixels that have the same r, g, and b as your mineral
        elg = (abs(Ig - minRGB(2))<1e-3);
        elb = (abs(Ib - minRGB(3))<1e-3);
        M = find(logical((elr.*elg.*elb)));  %now actually find those pixels
        num_Min(l-1) = length(M);

        ar = [-1 0 1];
        elM = [];

        for k = 1:length(M)   % for each pixel loop around the pixel to figure out what the mineral is, and stor it in an array

        c = 1;
        Mk = M(k);
        rc = ind2rc(size(Irgb,1),Mk);  % here is a function you created to turn the index into row column

        rw = int32(rc(1));
        cl = int32(rc(2));

            for j = 1:3  % here you run a cricle around your pixel of choice

                cl_j = (cl+ar(j));

                for i = 1:3

                    rw_i = (rw+ar(i));
                    rgbAr = [Irgb_loop(rw_i,cl_j,1),Irgb_loop(rw_i,cl_j,2),Irgb_loop(rw_i,cl_j,3)]; % get the rgb of the ith and jth pixel
                    diffMtx = min_col-rgbAr;
                    el_k = find(sum(abs(diffMtx),2)==0); % see which one it matches with on the list

                    if isempty(el_k)==0
                        elM(k,c) = find(sum(abs(diffMtx),2)==0); % if there actually is a mineral that matches, grab it
                        c = c+1;
                    end

                end

            end

    %     Irgb_loop(rw,cl,1) = -1; %% turn your pixel black so that it doesn't get counted twice when another pixel of the same composition is neighbouring
    %     Irgb_loop(rw,cl,2) = -1; 
    %     Irgb_loop(rw,cl,3) = -1;

        end

        %% analyze the neighbour matrix

        if isempty(elM)==0
            nM = [elM(:,1:4),elM(:,6:9)];
            num_sg(l-1) = sum(sum(nM==1,2)==8);


            for p = l:length(mins)

                numI(l,p) = sum(sum(nM==p));
        %         keyboard

            end
        end
    end

    numI(1,2:end)=num_sg;
    save(['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\' img '\' f_img],'numI','num_Min')

end
return

%% load and add numI
% in this section you can load the numI from the various images, and get a
% total count. 

fileDir = ['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\' img]
[nms] = dir([fileDir '\*.mat'])
matNm = {nms.name}

for i = 1:length(matNm)
    load([fileDir '\' matNm{i}]);
    fullI = numI+fullI;
    num_MinAll = num_MinAll+num_Min;
end
M = fullI(2:end,2:end);
All_M = M+M';

nma = repmat(num_MinAll,[lenMin-1,1]);
fullI_R = (All_M./nma)*100;
n = lenMin-1;
full2 = All_M;
full2(1:n+1:n*n)=0;
fullI_N = (full2./repmat(sum(full2),[n,1]))*100;


%%% write it out to an excel file
xlswrite(['D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\neighbourPixels\' img '.xlsx'],All_M,1,'B2:AE31')
xlswrite(['D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\neighbourPixels\' img '.xlsx'],fullI_R,2,'B2:AE31')
xlswrite(['D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\neighbourPixels\' img '.xlsx'],fullI_N,3,'B2:AE31')