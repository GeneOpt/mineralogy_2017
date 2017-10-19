clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\revisedColors\sediment\MI8118-APR16\1 - 25 um\';
[nms] = dir([folder '\*.TIF'])
matNm = {nms.name}

for H = 1:length(matNm)

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
    

    %% Here you repeat the same process but with just the mineral that you care about
    for m = 2:length(mins)
        M_col = min_col(m,:)
        % turn all non M minerals to white
        I_Mr = (Irgb(:,:,1)==M_col(1));
        I_Mg = (Irgb(:,:,2)==M_col(2));
        I_Mb = (Irgb(:,:,3)==M_col(3));
        I_M =  1-(I_Mr.*I_Mg.*I_Mb);

        %% grid search for all grains of mineral M
        %%% here you cruise through each pixel and find the isolated clusters
        %%% of mineral M

        [num_r,num_c] = size(I);
        ar_i = [-1 0 1 1 1 0 -1 -1];
        ar_j = [-1 -1 -1 0 1 1 1 0];

        I_mtxM = zeros(num_r,num_c);
        I_ackM = zeros(num_r,num_c);
        count = 1;

        % in this first loop you go through and row by column and assign colored
        % pixels a tag. First, you find the colored pixel, then cirlce to see if
        % there is a neighbouring colored pixel with a previous tag, and if not, give
        % a new one
        for j = 2:num_c-1

            for i = 2:num_r-1

                if I_M(i,j) == 0

                    for k = 1:8

                        i_in = i+ar_i(k);
                        j_in = j+ar_j(k);
                        pix(k,:) = [i_in,j_in];
                        I_mtxM_k(k) = I_mtxM(i_in,j_in);

                    end

                    elG = find(I_mtxM_k>0);

                    if elG
        %                 I_mtx_k(elG(1))
                        I_mtxM(i,j) = I_mtxM_k(elG(1));
                    elseif isempty(elG)
        %                 keyboard
                        I_mtxM(i,j) = count;
                    end

                    count = count+1;

                end
            end
        end

        num_u = []
        num_u(1) = length(unique(I_mtxM));

        % because the first loop can't get it right, you go through each tagged
        % pixel to see if a neighbouring pixel has a lower value, and grab it. Do
        % this enough times such that there is no longer a change in the number of
        % tags
        d_numu(1) = 10
        u = 1
        while d_numu > 1

            for j = 2:num_c-1

                for i = 2:num_r-1

                    I_ind = I_mtxM(i,j);

                    if I_ind~=0

                        for k = 1:8

                            i_in = i+ar_i(k);
                            j_in = j+ar_j(k);
                            pix(k,:) = [i_in,j_in];
                            a(k) = I_mtxM(i_in,j_in);

                        end

                        new_ind = min(a(find(a>0)));

                        if new_ind
                            I_mtxM(i,j) = new_ind;
                        end

                    end
                end
            end

        num_u(u+1) = length(unique(I_mtxM));
        d_numu = abs(num_u(u+1)-num_u(u));
        u = u+1;
        end


        % here you can give each tagged pixel a sequentially lower tag with no
        % missing values
        xUM = unique(I_mtxM);
        xArr = 0:length(xUM)-1;
        for i = 1:length(xUM)
            I_mtxM(I_mtxM==xUM(i))=xArr(i);
        end
        I_mtxM_all.(mins{m}) = I_mtxM;

    end
    
    save(['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\sed\' img(1:5) '\' img '.mat'],'I_mtxM_all','-append')

end












