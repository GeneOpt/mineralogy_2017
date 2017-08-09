%% here you should just make some triangles
% when you're setting up these test trials be sure that none of the
% elements are repeated in your shape
close all
clear variables
f2 = figure(2)

hArr = [0.1 0.5 1 2 10]
for hA = 1:length(hArr)
    w = 2;
    h = hArr(hA);
    np = 100
    x1 = linspace(0,w,np);
    y1 = zeros(1,np)+1;
    x2 = linspace(w,w/2,np/2);
    y2 = (2/3)*((-h.*x2/w)+h)+1;
    x3 = linspace(w/2,0,np/2);
    y3 = y2(end:-1:1);
    x = [x1,x2(2:end),x3(2:end-1)];
    y = [y1,y2(2:end),y3(2:end-1)];

    figure(2)
%     subplot(2,1,1)
    mX = mean(x)
    mY = mean(y)
%     hold on
%     
%     subplot(2,1,1)
%     plot(mX,mY,'ro')
    
    dx = x-mX;
    dy = y-mY;
    
    thetaD = atan(dy./dx);
    thetaD(dx<0) = thetaD(dx<0)+pi;
    thetaD(dy<0 & dx>0) = thetaD(dy<0 & dx>0)+2*pi;
    R = sqrt(dx.^2+dy.^2);
    subplot(2,1,1)
%     plot(thetaD,R)

    [theta,sT] = sort(thetaD);
    R = R(sT);
    xs = x(sT);
    ys = y(sT);
    
%     subplot(2,1,1)
%     hold on
%     figure
%     plot(xs,ys,'o')
%     hold on
%     plot(xs(15:16),ys(15:16),'+')

    dthetaO2 = (theta(2)-theta(1))/2
    thetaN = linspace(0+dthetaO2,2*pi+dthetaO2,np);
    
    for i = 1:length(thetaN)
    
        tN = thetaN(i);
        elB = find((theta-tN)>0);
        
        if isempty(elB)
            elL = length(x);
            elB = 1;
            msh = (ys(elB)-ys(elL))/(xs(elB)-xs(elL));
            bsh = (ys(elL)-msh*(xs(elL)));
            ml = tan(tN);
            bl = mY-ml*mX;
            xI(i) = ((bl-bsh)/(msh-ml));
            yI(i) = msh*xI(i)+bsh;
            elB = NaN;
        
        elseif elB(1) == 1
            
            elB = elB(1)
            elL = length(x);
            msh = (ys(elB)-ys(elL))/(xs(elB)-(xs(elL)));
            bsh = ys(elB)-msh*xs(elB);
            ml = tan(tN);
            bl = mY-ml*mX;
            xI(i) = (bl-bsh)/(msh-ml);
            yI(i) = msh*xI(i)+bsh;

        else 
            
            elB = elB(1);
            elL = elB-1;
            msh = (ys(elB)-ys(elL))/(xs(elB)-xs(elL));
            bsh = ys(elB)-msh*xs(elB);
            ml = tan(tN);
            bl = mY-ml*mX;
            xI(i) = (bl-bsh)/(msh-ml);
            yI(i) = msh*xI(i)+bsh;
   
        end


    end
    
%     subplot(2,1,1)
%     hold on
%     plot(xI,yI,'k+')
    
    dx = xI-mX;
    dy = yI-mY;
    
    thetaD = atan(dy./dx);
    thetaD(dx<0) = thetaD(dx<0)+pi;
    thetaD(dy<0 & dx>0) = thetaD(dy<0 & dx>0)+2*pi;
    R = sqrt(dx.^2+dy.^2);


    [theta,sT] = sort(thetaD);
    R = R(sT);
    xIs = xI(sT);
    yIs = yI(sT);
    
    subplot(2,1,1)
    hold on
    plot(theta,R);
    hold on
    plot(theta(1),R(1),'ro')
    hold on
    plot(theta(end),R(end),'ro')
    

    NFFT = length(theta);
    g = fft(R,NFFT);
    g = g(1:length(g)/2);
    modes = [1:length(g)];

    % NLi = NFFT;
    % g( modes > NLi ) = 0;
    filg = ifft(g);
    % return
    % 
    % figure(2)
    % hold on
    % plot(theta,filg,'ko-')

    D = abs(g)
    phi = unwrap(angle(g))

    y = log(D./max(D))
    subplot(2,1,2)
    hold on
    p(hA) = plot(modes,D/max(D),'-','linewidth',1)
 
end

ylabel('Normalized power')
xlabel('Harmonic (n)')
% legend([p(1) p(2) p(3) p(4)],{'AR = 0.05', 'AR = 0.25','AR = 0.5','AR = 1'})
grid on
% savePDFfunction(f2, 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\grain shape\rectangles')
