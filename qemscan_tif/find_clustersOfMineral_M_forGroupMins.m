%% concatenate mnrlMtx
% in this section the mnrlMtx is operated on, and all images of the same
% sample are concatenated and placed in the 'rock' or 'seds' folder
clear all
run mineral_colors.m;
scf = 25/38;
    
MAbb_i = [4 6 8 11 13 18]
combMin = {'Plag_A','Biot_A','Chlor_A','Ill_smec_A','Dol_A','Ilm_A'}
combMinP = {'Plag A','Biot A','Chlor A','Ill smec A','Dol A','Ilm A'}
% MAbb.(combMin{1}) = [4 5 6 7]
% MAbb.(combMin{2}) = [9 10 11 12]
% MAbb.(combMin{3}) = [14 15]
% MAbb.(combMin{4}) = [18 19]
% MAbb.(combMin{5}) = [21 22]
% MAbb.(combMin{6}) = [27 28]

% get all the folder names (each folder is a rock sample)
folder = 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_imDat_revisedCol\grain_basics\rock\';
[nms] = dir([folder]);
matNm = {nms.name}

for F = 3:28 %rock
% for F = 3:23 %sed

    fName = dir([folder matNm{F} '\*.mat'])
    fN = {fName.name}
    
    for N = 1:length(fN)
        
        fullName = [folder matNm{F} '\' fN{N}]
        load(fullName);
    
return
    %% Here you repeat the same process but with just the mineral that you care about
        for m = 1:length(MAbb_i)
            M_col = m+39
            I_m = I_mtx_C_abbv.(combMin{m});
            % turn all non M minerals to white

            %% grid search for all grains of mineral M
            %%% here you cruise through each pixel and find the isolated clusters
            %%% of mineral M

            [num_r,num_c] = size(I_m);
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

                    if I_m(i,j) == M_col

                        for k = 1:8

                            i_in = i+ar_i(k);
                            j_in = j+ar_j(k);
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
            I_mtxM_all.(combMin{m}) = I_mtxM;
        end
        save(fullName, 'I_mtxM_all','-append')
    end
    
end












