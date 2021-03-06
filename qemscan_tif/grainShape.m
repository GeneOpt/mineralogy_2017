clc
close all
clear variables

run mineral_colors

folder = 'D:\Field_data\2013\Summer\Geochemistry\qemscan_edited\images\revisedColors\rock_frags\MI8152-AUG16_needs to ne redone\1 - 25 um\';
[nms] = dir([folder '\*.TIF'])
matNm = {nms.name}

% for H = 1:length(matNm)-1
for H = 2:3
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


    
    %% grid search for all grains
    %%% here you cruise through each pixel of interest, and see what minerlas
    %%% are touching the other minerals. Here you don't really care what
    %%% the minerals are, you're just finding isolated grains

    lenMin = length(mins);
    [num_r,num_c] = size(I);
    ar_i = [-1 0 1 1 1 0 -1 -1];
    ar_j = [-1 -1 -1 0 1 1 1 0];

    I_mtx = zeros(num_r,num_c);
    I_ack = zeros(num_r,num_c);
    count = 1;

    % in this first loop you go through and row by column and assign colored
    % pixels a tag. First, you find the colored pixel, then cirlce to see if
    % there is a neighbouring colored pixel with a previous tag, and if not, give
    % a new one
    for j = 2:num_c-1

        for i = 2:num_r-1

            rgbA = [Irgb(i,j,1),Irgb(i,j,2),Irgb(i,j,3)]; % get the rgb of the ith and jth pixel

            if sum(rgbA)<3

                for k = 1:8

                    i_in = i+ar_i(k);
                    j_in = j+ar_j(k);
                    pix(k,:) = [i_in,j_in];
                    I_mtx_k(k) = I_mtx(i_in,j_in);

                end

                elG = find(I_mtx_k>0);

                if elG
    %                 I_mtx_k(elG(1))
                    I_mtx(i,j) = I_mtx_k(elG(1));
    %                 keyboard
                elseif isempty(elG)
                    I_mtx(i,j) = count;
                end

                count = count+1;

            end
        end
    end

    num_u = []
    num_u(1) = length(unique(I_mtx));

    % because the first loop can't get it right, you go through each tagged
    % pixel to see if a neighbouring pixel has a lower value, and grab it. Do
    % this enough times such that there is no longer a change in the number of
    % tags
    d_numu(1) = 10
    u = 1
    while d_numu > 1

        for j = 2:num_c-1

            for i = 2:num_r-1

                I_ind = I_mtx(i,j);

                if I_ind~=0

                    for k = 1:8

                        i_in = i+ar_i(k);
                        j_in = j+ar_j(k);
                        pix(k,:) = [i_in,j_in];
                        a(k) = I_mtx(i_in,j_in);

                    end

                    new_ind = min(a(find(a>0)));

                    if new_ind
                        I_mtx(i,j) = new_ind;
                    end

                end
            end
        end

    num_u(u+1) = length(unique(I_mtx));
    d_numu = abs(num_u(u+1)-num_u(u));
    u = u+1;
    end

    figure
    plot(1:u,num_u,'-o')


    % here you can give each tagged pixel a sequentially lower tag with no
    % missing values
    xU = unique(I_mtx);
    xArr = 0:length(xU)-1;
    for i = 1:length(xU)
        I_mtx(I_mtx==xU(i))=xArr(i);
    end
    
    %% size and mineralogy
    clc
    mnrlMtx = zeros(length(xU),length(mins));

    for i = 1:length(xU)
    %     figure
    %     testRGB(I_mtx,i)
        sizeMtx(i) = sum(sum(I_mtx==i));

        M = find(I_mtx==i);

        M_Arr = zeros(1,length(mins));
        for k = 1:length(M)

            Mk = M(k);
            rc = ind2rc(size(Irgb,1),Mk);
            rw = int32(rc(1));
            cl = int32(rc(2));   
            rgbAr = [Irgb(rw,cl,1),Irgb(rw,cl,2),Irgb(rw,cl,3)];
            diffMtx = min_col-rgbAr;
            el_k = find(sum(abs(diffMtx),2)==0);
            M_Arr(el_k) = M_Arr(el_k)+1;

        end
        mnrlMtx(i,:)=M_Arr;

    end



    %% here you start to compute grain shape
    fs4 = 12
    close all
    nRow = size(I_mtx,1);
    minSizeLim = 25
    fracMLim = 0.8    
    countP = 1;
    TOS = 0;

    for i = 1:length(xU)-1

        G_ind = find(I_mtx==i); % isolate the grain that you are looking for

        if length(G_ind)>minSizeLim

        %         for CM = 1:2

                [rc] = ind2rc(nRow,G_ind); % find the row column from the index
                rc = int32(rc);
                n_p = length(rc);
        %             if CM == 1
                x_bar = sum(rc(:,1))/n_p; % find the center of mass of x
                y_bar = sum(rc(:,2))/n_p; % find the center of mass of y
                mP = [x_bar,y_bar];
        %             elseif CM == 2
        %                 x_bar = x_bar2
        %                 y_bar = y_bar2
        %                 mP = [x_bar2,y_bar2];
        %             end
                B = (I_mtx==i);   % turn the particle that you want black
                BP = bwperim(B);    % get matlba to help your laziness outline the ptcl perimeter
                bpInd = find(BP);  % find the perimeter elements
                [rcP] = ind2rc(nRow,bpInd);  % turn perimeter elemets into row column


        %             close all
        %             figure(1)
        %             subplot(2,2,1)
        %     %         plot(rc(:,1),rc(:,2),'k.')
        %             hold on 
        %             p1 = plot(x_bar,y_bar,'ks','markerfacecolor','k')
        %             hold on
        %             grid on
        %         plot(rcP(:,1),rcP(:,2),'r+')


            %     hold on
            %     plot(rcP(ordArr(1:28),1),rcP(ordArr(1:28),2),'b')

                px = rcP(:,1); % get the x coords of perimeters
                py = rcP(:,2); % get the y coords of perimeters
                pX = repmat(px,[1,length(px)]);  % make a matrix of all x values
                pY = repmat(py,[1,length(py)]);  % make a matrix of all y values
                ordArr = 1;
                dMtx = sqrt((pX-pX').^2+(pY-pY').^2);  % compute distance between all perimeter points
                for j = 1:length(px)-1  % make a loop that goes through all points and finds the closest neighboring perimeter point that is not already in the new set

                    nexPix = ordArr(j);
                    dj_col = dMtx(:,nexPix);
                    clst = (min(dj_col(dj_col>0)));
                    elMin = find(dj_col==clst);

                    if length(elMin>1)

            %            elMY =  min(py(elMin));
                       elMin = max(elMin);

                    end

                    ordArr(j+1) = elMin;
                    dMtx(ordArr,ones(1,length(j))*elMin) = 0;
                end

        %             hold on
                rcP_c = rcP(ordArr,:);
                rcP_c(end+1,:) = rcP_c(1,:);
        %             p2 = plot(rcP_c(:,1),rcP_c(:,2),'b+-');
        %             xlabel('x position','fontsize',fs4)
        %             ylabel('y position','fontsize',fs4)
        %             axis square
        %             xl1 = get(gca,'xlim');
        %             yl1 = get(gca,'ylim');
        %             lmax = max([(xl1(2)-xl1(1)), (yl1(2)-yl1(1))]);
        %             xlim([xl1(2)-lmax xl1(2)]);
        %             ylim([yl1(2)-lmax yl1(2)]);

                d_mP = mP-rcP_c;  % find the distance from the center of mass to each perimeter point
                R_j = sqrt(sum(d_mP.^2,2));
                theta_j = (atan((rcP_c(:,2)-y_bar)./(rcP_c(:,1)-x_bar)));  % find the angle from the center of mass to each perimeter point
                elz_n = find((diff(theta_j))<-2);  % find the differences in angles that jump from theta to pi-theta as x<0
                elz_p = find((diff(theta_j))>2);   % and repeat
        %         elz_s = [elz_n;elz_p]
        %         elz = sort(elz_s)

        %             figure(2)
        %             subplot(2,2,1)
        %             plot(theta_j)
        %             subplot(2,2,3)
        %             plot(diff(theta_j))

                for dt = 1:length(elz_n)  % correct for the difference in pi
                    theta_j([elz_n(dt)+1:end]) = theta_j([elz_n(dt)+1:end])+pi;
                end

                for dt = 1:length(elz_p)
                    theta_j([elz_p(dt)+1:end]) = theta_j([elz_p(dt)+1:end])-pi;
                end

        %         theta_j = abs(theta_j-theta_j(1))    % reset all values so that the first point is theta=0
                if sum(theta_j>0)==0
                    theta_j = abs(theta_j);
                    g = 1;
                else
                    g = 0;
                end


                dt1 = theta_j([2:end])-theta_j(1);
                if sum(dt1>0)==0
                    theta_j = theta_j(end:-1:1);
                    R_j = R_j(end:-1:1);
                    g = 2;
                end

                elD = 1;
                el2 = 1;
                count = 1;
                while el2 < length(theta_j)-1  % remove double valued points from data set.
                    dt = theta_j([el2+1:end])-theta_j(el2);
                    el1 = find(dt>0);

                    if isempty(el1) ~= 1            
                        el2 = el1(1)+el2;
                        elD(count+1) = el2;
                        count = count+1;
                    else
                        el2 = length(theta_j)+1;
                    end

                end

        %             subplot(2,2,2)
        %             plot(theta_j)
        %             subplot(2,2,4)
        %             plot(diff(theta_j))

        %             figure(1)
        %             subplot(2,2,3)
        %             p3 = plot(theta_j,R_j,'b+');
        %             hold on
                tjD = theta_j(elD);
                if length(tjD)>2
                RjD = R_j(elD);
        %         plot(tjD,RjD,'mo')
        %             xlabel('Angle (radian)','fontsize',fs4)
        %             ylabel('Radius (pixels)')
        %             grid on

        %         subplot(2,2,3)
        %         hold on
        %         plot(rcP_c(elD,1),rcP_c(elD,2),'m','linewidth',4)


                yy = [];
                xx = linspace(tjD(1)+pi/1000,tjD(end),100); % here you are doing a linear interpolation of the points to get evenly spaced theta valus
                for s = 1:length(xx)
                    x_i = xx(s);
                    firstEl = find((x_i-tjD)>0);
                    el1 = firstEl(end);
                    dx = (tjD(el1+1)-tjD(el1));
                    dy = RjD(el1+1)-RjD(el1);
                    m = dy/dx;
                    b = RjD(el1) - m*tjD(el1);
                    yy(s) = m*(x_i)+b;
                end




        %             figure(1)
        %             subplot(2,2,3)
        %             p4 = plot(xx,yy,'r.');

                dxe = (xx(end)-xx(1));
                dxeLim = (2*pi)-(pi/4);

                if dxe < dxeLim
                     TOS = TOS+1
        %                 keyboard
                else

        %                 xIntp = yy.*cos(xx)+x_bar  %% then you put the radial coords back back into cartesian coords
        %                 yIntp = yy.*sin(xx)+y_bar
            %         
        %                 x_bar2 = sum(xIntp)/length(xIntp); % find the center of mass of x
        %                 y_bar2 = sum(yIntp)/length(yIntp); % find the center of mass of y
        %                 mP2 = [x_bar2,y_bar2];

        %                 figure(1)
        %                 subplot(2,2,1)
        %                 hold on
        %                 plot(x_bar2,y_bar2,'rs')

        %                 if isempty(find(tjD<0)) && g == 0  %% and do some weird stuff to correct the points based on the angles. trial and error
        % 
        %                     subplot(2,2,1)
        %                     p5 =  plot(2*x_bar-xIntp,2*y_bar-yIntp,'r.')
        % 
        %                 elseif isempty(find(tjD<0)) && g == 1
        % 
        %                     subplot(2,2,1)
        %                     p5 =  plot(xIntp,2*y_bar-yIntp,'r.')
        % 
        % 
        %                 elseif isempty(find(tjD<0)) && g == 2
        % 
        %                     subplot(2,2,1)
        %                     p5 =  plot(2*x_bar-xIntp,yIntp,'r.')
        % 
        %                 else
        % 
        %                     subplot(2,2,1)
        %                     p5 =  plot(xIntp,yIntp,'r.')
        % 
        %                 end

        %                 legend([p1 p2 p5], {'center of mass','grain outline','interpolated outline'},'location','southeast')

        %                 figure(3)
        %                 th = linspace(0,2*pi,50)
        %                 T = th(2)-th(1)
        %                 Fs = 1/T
                        lx = length(xx);
        %                 jf = 0:lx-1; 
        %                 kArr = [1:20]

        %                 YF = fft(yy)
        %                 P2 = abs(YF/lx)
        %                 P1 = P2(1:lx/2+1)
        %                 P1(2:end-1) = 2*P1(2:end-1)
        %                 fF = Fs*(0:(lx/2))/lx;
            %         subplot(2,2,3)
        %                 hold on
        %                 plot(fF*2*pi,P1)
        % 
        %                 figure(1)
        %                 subplot(2,2,4)
        %                 plot(1:length(yy), yy, 'r.','markersize',10)
        % 
        %                 NLArr = [1,5,10,50]
        %                 colr = {'r','m','b','k'}

        %                 for NL = 1:length(NLArr)
        %                     NLi = NLArr(NL)
            %             
            %             Y = zeros(1,length(th))
            %             for N = 1:length(kArr)
            %                 Y = (A_k(N)*cos(N*th) + B_k(N)*sin(N*th))+Y
            %             end

                        NFFT = lx;
                        g = fft(yy,NFFT);
                        g = g(1:length(g)/2);
        %                     realg = abs(g);
        %                     angleg = unwrap(angle(g));
                        modes = [1:length(g)];
                        NLi = 50;
                        g( modes > NLi ) = 0;
                        filg = ifft(g,NFFT,'symmetric');

        %                     figure(1)
        %                     subplot(3,1,1)
        %                     plot(1:length(yy),filg,'k-')
        %                     xlabel('j')
        %                     grid on
        % 
        %                     xIntpF = filg.*cos(xx)+x_bar  %% then you put the radial coords back back into cartesian coords
        %                     yIntpF = filg.*sin(xx)+y_bar
        % % 
        %                     subplot(3,1,2)
        %                     plot(xIntpF,yIntpF,'k-')
        %                     grid on
        %                     
        %                     subplot(3,1,3)
        %                     plot(modes,g,'k-')
        %                     keyboard

                        %                 end

        %                 f1 = figure(1)
        %                 subplot(2,2,1)
        %                 title('Raw & intp. data')
        %                 axis square
        %                 xlim([xl1(2)-lmax xl1(2)])
        %                 ylim([yl1(2)-lmax yl1(2)])
        % 
        %                 subplot(2,2,2)
        %                 axis square
        %                 xlim([xl1(2)-lmax xl1(2)])
        %                 ylim([yl1(2)-lmax yl1(2)])
        %                 title('Spectral decompostion')
        %                 xlabel('x position','fontsize',fs4)
        %                 ylabel('y position','fontsize',fs4)
        %                 legend(ph,{'K = 1','K = 5','K = 10','K = 50',},'location','southeast')
        % 
        %                 subplot(2,2,4)
        %                 ylabel('Radius (pixels)')
        %                 f50 = figure(50)
        %                 hold on
        %                 plot(modes, log(abs(filterg)./max(abs(filterg))),'linewidth',1,'color','k')
        %                 xlim([0 50])
        %                 grid on
        %                 ylabel('Normalized power','fontsize',16)
        %                 xlabel('Harmonic (N)','fontsize',16)

                    fSpec(countP,:)= g;
                    mnrlMtxP(countP,:) = mnrlMtx(i,:);
                    countP = countP+1;
                end
             end

        end

    end
    save(['D:\Code\Summer_2013_data\mineral_data\qemscan_tif\sample_image_data\shape\rock\' img '.mat'],'fSpec','mnrlMtxP','I_mtx')
end
