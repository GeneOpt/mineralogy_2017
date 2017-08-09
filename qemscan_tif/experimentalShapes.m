clear all
close all

A_k = [5 3 1 0.1]
B_k = [5 3 1 0.1]
f_k = [0 1 2 3]
f_k = f_k/(2*pi)

theta = linspace(0,2*pi,400);
yT = 0

for i = 1:length(A_k)
    figure(1)
    subplot(3,2,1)
    yA = A_k(i)*cos(2*pi*theta*f_k(i))
    hold on
    plot(theta,yA,'-')
    
    subplot(3,2,2)
    hold on
    xxA = yA.*cos(theta)
    yyA = yA.*sin(theta)
    plot(xxA,yyA)
    axis equal
    
    subplot(3,2,3)
    yB = B_k(i)*sin(2*pi*theta*f_k(i))
    hold on
    plot(theta,yB,'-')
    
    subplot(3,2,4)
    hold on
    xxB = yB.*cos(theta)
    yyB = yB.*sin(theta)
    plot(xxB,yyB)
    axis equal
    
    subplot(3,2,5)
    hold on
    yT = yA+yB+yT
    hold on
    plot(theta,yT,'-')
    
    subplot(3,2,6)
    hold on
    xxT = yT.*cos(theta)
    yyT = yT.*sin(theta)
    plot(xxT,yyT)
    axis equal
    
end


%% here you should just make a square and some recatngles
close all
clear variables
f2 = figure(2)

hArr = [0.1 0.5 1 2]
hArr = 1
for hA = 1:length(hArr)
    w = 2;
    h = hArr(hA);
    np = 50
    x1 = zeros(1,np);
    y1 = linspace(h,0,np);
    x2 = linspace(0,w,np);
    y2 = zeros(1,np);
    x3 = ones(1,np)*w;
    y3 = linspace(0,h,np);
    x4 = linspace(w,0,np);
    y4 = ones(1,np)*h;
    x = [x1,x2,x3,x4];
    y = [y1,y2,y3,y4];

    figure(10)
    plot(x,y,'-')

    mX = mean(x)
    mY = mean(y)
    hold on

%     plot(mX,mY,'ro')
    R = sqrt((x-mX).^2+(y-mY).^2)
    theta = atan((y-mY)./(x-mX))
    dtheta = theta(2:end)-theta(1:end-1)
    elD = find(dtheta<0)
    thetaO = theta
    for dth = 1:length(elD)
        theta([elD(dth)+1:end]) = theta([elD(dth)+1:end])+pi;
    end
%            
%     plot(theta,R,'o')
    yy = [];
    xx = linspace(theta(1)+pi/1000,theta(end),50); % here you are doing a linear interpolation of the points to get evenly spaced theta valus
    for s = 1:length(xx)
        x_i = xx(s);
        firstEl = find((x_i-theta)>0);
        el1 = firstEl(end);
        dx = (theta(el1+1)-theta(el1));
        dy = R(el1+1)-R(el1);
        m = dy/dx;
        b = R(el1) - m*theta(el1);
        yy(s) = m*(x_i)+b;
    end
    
    figure(2)
    subplot(2,1,1)
    hold on
    plot(theta,R,'-','linewidth',1)
    hold on
    plot(xx,yy,'m+')
    
    xInt = yy.*cos(xx)+mX
    yInt = yy.*sin(xx)+mY
    figure(10)
    hold on
    plot(xInt,yInt,'k+')
    
    
    return
    ylabel('Radius')
    grid on
    xlabel('\theta')
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
    p(hA) = plot(modes(1:10),D(1:10)/max(D),'-','linewidth',1)
 
end

ylabel('Normalized power')
xlabel('Harmonic (n)')
legend([p(1) p(2) p(3) p(4)],{'AR = 0.05', 'AR = 0.25','AR = 0.5','AR = 1'})
grid on
% savePDFfunction(f2, 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\grain shape\rectangles')

%% here you should just make some ellipses
% close all
% clear variables
% f2 = figure(2)

hArr = [0.1 0.25 0.75 1]
for hA = 1:length(hArr)
    a2 = 1;
    b2 = hArr(hA);
    np = 60
    x1 = linspace(0,1,np)
    x2 = -x1
    y1 = sqrt(b2*(1-((x1.^2)/(a2))))
    y2 = -y1
    x = [x2,x1,x1,x2]
    y = [y1,y1,y2,y2]

%     plot(x,y,'o-')

    mX = mean(x)
    mY = mean(y)
    hold on

%     plot(mX,mY,'ro')
    R = sqrt((x-mX).^2+(y-mY).^2)
    theta = atan((y-mY)./(x-mX))
    [theta,iT] = sort(theta)
    R = R(iT)

    theta = [theta,theta+pi]
    R = [R,R]
    
    figure(2)
    subplot(2,1,1)
    hold on
    plot(theta,R,'-','linewidth',2)
    ylabel('Radius')
    grid on
    xlabel('\theta')
% return
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
    p2 = plot(modes(1:10),D(1:10)/max(D),'-','linewidth',2)
 
end

subplot(2,1,1)
set(gca,'fontsize',16)
subplot(2,1,2)
ylabel('Normalized power')
xlabel('Harmonic (n)')
legend([p1 p2], {'rectangles','ovals'})
% legend([p(1) p(2) p(3) p(4)],{['b^2 = ' num2str(hArr(1))],['b^2 = ' num2str(hArr(2))],['b^2 = ' num2str(hArr(3))],['b^2 = ' num2str(hArr(4))]})
grid on
set(gca,'fontsize',16)
% savePDFfunction(f2, 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\grain shape\oANDr')

%% here you should just make some noise over the ellipse
close all
clear variables
f2 = figure(2)

hArr = [1]

for hA = 1:length(hArr)
    a2 = 1;
    b2 = 0.1;
    np = 60
    x1 = linspace(0,1,np)
    x2 = -x1
    y1 = sqrt(b2*(1-((x1.^2)/(a2))))
    y2 = -y1
    x = [x2,x1,x1,x2]
    y = [y1,y1,y2,y2]

    figure(1)
    plot(x,y,'.')

    mX = mean(x)
    mY = mean(y)
    hold on

    plot(mX,mY,'ro')
    R = sqrt((x-mX).^2+(y-mY).^2)
    theta = atan((y-mY)./(x-mX))
    [theta,iT] = sort(theta)
    R = R(iT)

    theta = [theta,theta+pi]
    R = [R,R]
    
    figure(2)
    subplot(2,1,1)
    hold on
    plot(theta,R,'k-','linewidth',2)
    ylabel('Radius')
    grid on
    xlabel('\theta')
    
    nn = 0.05
    rR = (rand(1,length(R))*nn)-nn/2. + R
    hold on
    plot(theta,rR,'b-','linewidth',1)
    
    xx = rR.*cos(theta)+mX;
    yy = rR.*sin(theta)+mY;
    
    figure(1)
    hold on
    plot(xx,yy,'b-')

    NFFT = length(theta);
    g = fft(R,NFFT);
    g = g(1:length(g)/2);
    modes = [1:length(g)];
    
    gR = fft(rR,NFFT);
    gR = gR(1:length(gR)/2);
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
    Dr = abs(gR)
    phi = unwrap(angle(g))

    y = log(D./max(D))
    figure(2)
    subplot(2,1,2)
    hold on
    p2 = plot(modes(1:20),log(D(1:20)/max(D)),'k-','linewidth',2)
    hold on
    p3 = plot(modes(1:20),log(Dr(1:20)/max(Dr)),'b-','linewidth',1)
 
end

ylabel('log(Normalized power)')
xlabel('Harmonic (n)')
% legend([p1 p2], {'rectangles','ovals'})
% legend([p(1) p(2) p(3) p(4)],{['b^2 = ' num2str(hArr(1))],['b^2 = ' num2str(hArr(2))],['b^2 = ' num2str(hArr(3))],['b^2 = ' num2str(hArr(4))]})
grid on
savePDFfunction(f2, 'D:\Code\Summer_2013_data\mineral_data\qemscan_tif\figures\grain shape\noise')

%% here you should just make some triangles
close all
clear variables
f2 = figure(2)

hArr = [0.1 0.5 1 4 10]
hArr = 2
for hA = 1:length(hArr)
    w = 2;
    h = hArr(hA);
    np = 40
    x1 = linspace(0,w,np);
    y1 = zeros(1,np)+1;
    x2 = linspace(w,w/2,np/2);
    y2 = (2/3)*((-h.*x2/w)+h)+1;
    x3 = linspace(w/2,0,np/2);
    y3 = y2(end:-1:1);
    x = [x1,x2,x3];
    y = [y1,y2,y3];

    figure(2)
    subplot(2,1,2)
    mX = mean(x)
    mY = mean(y)
    hold on

    plot(mX,mY,'ro')

    for i = 1:length(x)

        xi = x(i) 
        yi = y(i)
        dx = xi-mX
        dy = yi-mY


        if dy>0 & dx> 0 
            thetaD(i) = atan(dy/dx)
            mc = 'k'
            'FIRST'
        elseif dy>0 & dx<0
            thetaD(i) = pi +  atan(dy/dx)  
            mc = 'r'
            'SECOND'
        elseif dy<0 & dx<0
            thetaD(i) = atan(dy/dx) + pi
            mc = 'b'
            'THIRD'
        elseif dy<0 & dx>0
            thetaD(i) = 2*pi + atan(dy/dx) 
            mc = 'g'
            'FOURTH'
        end
        
        subplot(2,1,2)
        hold on
        plot(xi,yi,'o-','color',mc)
        axis equal
        
        R(i) = sqrt(dx.^2+dy.^2);
        subplot(2,1,1)
        hold on
        plot(thetaD(i),R(i),'o','color',mc)

    end

    [theta,sT] = sort(thetaD)
    R = R(sT)
    xs = x(sT)
    ys = y(sT)
    
    thetaN = linspace(0+pi/200,2*pi+pi/200,100)
    
    for i = 1:length(thetaN)
    
        tN = thetaN(i)
        elB = find((theta-tN)>0)
        
        if isempty(elB)
            elL = length(x)
            elB = 1
            msh = (ys(elB)-ys(elL))/(xs(elB)-xs(elL))
            bsh = (ys(elL)-msh*(xs(elL)))

            ml = tan(tN)
            bl = mY-ml*mX
            xI(i) = ((bl-bsh)/(msh-ml))
            yI(i) = msh*xI(i)+bsh
            'end'

            elB = NaN
        
        elseif elB(1) == 1

            elL = length(x)
            msh = (ys(elB)-ys(elL))/(xs(elB)-(xs(elL)-2*pi))
            bsh = ys(elB)-msh*xs(elB) 

            ml = tan(tN)
            bl = mY-ml*mX
            xI(i) = (bl-bsh)/(msh-ml)
            yI(i) = msh*xI(i)+bsh
            'first'


        else 
            
            elB = elB(1)
            elL = elB-1
            msh = (ys(elB)-ys(elL))/(xs(elB)-xs(elL))
            bsh = ys(elB)-msh*xs(elB) 
            
            ml = tan(tN)
            bl = mY-ml*mX
            xI(i) = (bl-bsh)/(msh-ml)
            yI(i) = msh*xI(i)+bsh
   
        end


    end
    
    subplot(2,1,2)
    plot(xI,yI,'k+')
    
    for i = 1:length(xI)

        xi = xI(i) 
        yi = yI(i)
        dx = xi-mX
        dy = yi-mY


        if dy>0 & dx> 0 
            thetaD(i) = atan(dy/dx)
            mc = 'k'
            'FIRST'
        elseif dy>0 & dx<0
            thetaD(i) = pi +  atan(dy/dx)  
            mc = 'r'
            'SECOND'
        elseif dy<0 & dx<0
            thetaD(i) = atan(dy/dx) + pi
            mc = 'b'
            'THIRD'
        elseif dy<0 & dx>0
            thetaD(i) = 2*pi + atan(dy/dx) 
            mc = 'g'
            'FOURTH'
        end
        
        subplot(2,1,2)
        hold on
        plot(xi,yi,'^-','color',mc)
        axis equal
        
        R(i) = sqrt(dx.^2+dy.^2);
        subplot(2,1,1)
        hold on
        plot(thetaD(i),R(i),'^','color',mc)

    end
    
    
    figure(20)
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

%% here's another partially failed attempt at fourier.

A_k = D.*cos(phi);
B_k = D.*sin(phi);
figure(4)
subplot(2,1,1)
stem(A_k)
ylabel('A_k')
subplot(2,1,2)
stem(B_k)
ylabel('B_k')

figure(5)
subplot(2,1,1)
stem(real(g))
ylabel('real(g)')
subplot(2,1,2)
stem(imag(g))
ylabel('imag(g)')

lx = length(x)
kArr = 1:(lx/2+1)
for kI = 1:length(kArr)
    A(kI) = (2/lx)*sum(R.*cos(2*pi*theta*kArr(kI)/lx));
    B(kI) = (2/lx)*sum(R.*sin(2*pi*theta*kArr(kI)/lx));
end

figure(6)
subplot(2,1,1)
stem(A)
ylabel('A')
subplot(2,1,2)
stem(B)
ylabel('B')

for t = 1:length(theta)
    tt = theta(t)
    y(t) = sum(A.*cos(2*pi*tt*kArr/lx)+B.*sin(2*pi*tt*kArr/lx))
end

figure(2)
hold on
plot(theta,y,'ko-')

return

for kI = 1:length(kArr)
    k = kArr(kI)
    B_k(kI) = (2/k)*sum(yy.*sin(2*pi*jf*k/lx))
    A_k(kI) = (2/k)*sum(yy.*cos(2*pi*jf*k/lx))
end

return
for tt = 1:length(theta)
    ttA = theta(tt)
    yt(tt) = sum((A_k.*cos(modes*ttA*(2*pi)))+(B_k.*sin(modes*ttA*(2*pi))))
end

figure(2)
hold on
plot(theta,yt,'r+-')

return
%%
clear variables
close all
W = 8
H = 3
t1 = linspace(-pi/4,pi/4,50);
del_t = t1(2)-t1(1)
R1 = W./(2*cos(t1))
plot(t1,R1)
t2 = linspace(pi/4,3*pi/4,50) 
R2 = H./(2*sin(t2))
hold on
plot(t2,R2)



%% see if you can find angles
x1 = [-1,0]
y1 = [1,2]
m1 = (y1(2)-y1(1))/(x1(2)-x1(1))

theta = -pi/6
m2 = tan(theta)
x2 = [0,2]
y2 = [0,m2*x2(2)]

close all
plot(x1',y1')
hold on
plot(x2',y2')
axis equal

intA = atan(m1)-theta
intA_d = abs(intA*180/pi)
text(0.8,0.9,num2str(intA_d),'units','normalized')











