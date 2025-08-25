%% ROUGHNESS AND NOISE REMOVAL
clc
clear all
close all
load('input_temp.mat');
Round_spheri_rough= [0 0 0 0 0 0 0];
save('Round_spheri_rough', 'Round_spheri_rough')
h= '.' ;
ho= '.';
res=glob(h,filename);
nit=size(res,1);

for P= 1:nit
    h= '.';
ho= '.';
    save('h','h')
    save('ho','ho')
    save('P','P')
    result=glob(h,filename);
    clr = lines(20);
    he = imread(result{P});% to read the current image
    level = graythresh(he);
    bw = im2bw(he,level);
    bw = bwareaopen(bw, 50);
    BW_filled = imfill(bw,'holes');
    im_final = bw;
    graindata = regionprops(im_final,'Area','Perimeter','Centroid','MajorAxisLength','MinorAxisLength','ConvexHull','Orientation','Eccentricity'); % using regionprops all the information is extracted from the image and stored in graindata matrix.
    [B,L] = bwboundaries(im_final,'noholes');
    figure, imshow(bw)
    imwrite(bw, 'output1.jpg');

    %  hold on
    %  plot(graindata.ConvexHull(:,1),graindata.ConvexHull(:,2), 'blue', 'LineWidth', 2)
    % impoint(gca, graindata.Centroid(1,1),graindata.Centroid(1,2))
    boundary =  B{1};
    %figure, plot(boundary(:,2), boundary(:,1), 'blue', 'LineWidth', 2)
    
    
   %% minimum circumscribing circle
    % finding the maximum length of chord within the grain to find the maximum distance of boundary point from centroid
    tic
    x = B{1}(:,2); % x component of coordinate
    y = B{1}(:,1); % y component of coordinate
    d=[];
    mo=[];
    dist_p=[];
    dist_pl=[];
    for i = 1:length(x) - 1
        d = [d; x(i) y(i) sqrt((graindata.Centroid(1,2) - x(i))^2 + (graindata.Centroid(1,2) - y(i))^2) i];% maximum distance from the centroied
    end
    [b,ind] = max(d(:,3));                   %Index of the max d
    x1= d(ind,1);                            % x component corresponds to max d
    y1=d(ind,2);                             % y component corresponds to max d
    index= d(ind,4);                         %Index of the max d
    
    % to find the maximum distance of a boundary point from the obtained boundary point
    
    for j = 1:length(x) - 1
        distance = sqrt((x1 - x(j))^2 + (y1 - y(j))^2);% distance from the max d to all points in boundary
        mo = [mo; x(j) y(j) distance j];                 % all points and corresponding distances
    end
    [c,ind1] = max(mo(:,3));                            %Index of the max m
    x2= mo(ind1,1);                                     %x component corresponds to max m
    y2=mo(ind1,2);                                      %y component corresponds to max m
    index1= mo(ind1,4);                                 %Index of the max m
    % figure, imshow(bw)
    % hold on
    % plot([x1, x2],[y1, y2],  'blue', 'LineWidth', 2)
    
    xc= (x1+x2)/2;
    yc= (y1+y2)/2;
    rc= sqrt((x1 - xc)^2 + (y1 - yc)^2);
    xnc= [x1;x2];
    ync=[y1;y2];
    xt=[];
    yt=[];
    dnc=[];
    % loop to check if boundary points are outside the circle and including
    % the farthest point in forming a circle until all boundary points are within the circle
    
    for i=1: size(graindata.ConvexHull,1)
        for j=1: size(graindata.ConvexHull,1)
            dfc= sqrt((xc - graindata.ConvexHull(j,1))^2 + (yc - graindata.ConvexHull(j,2))^2);
            if dfc > rc
                xt= [xt; graindata.ConvexHull(j,1)];
                yt= [yt; graindata.ConvexHull(j,2)];
                dnc=[dnc; dfc];
            end
        end
        if size(xt,1) == 0
            break
        end
        [c2,ind2] = max(dnc(:,1));
        xnc= [xnc; xt(ind2,1)];
        ync= [ync; yt(ind2,1)];
        [xc1,yc1,Re,a] = circfit(xnc,ync); % function to fit a circle to a given set of points
        xc=xc1;
        yc=yc1;
        rc=Re;
        xt=[];
        yt=[];
        dnc=[];
    end
    
    % syms xi yi
    % f(xi,yi)=(xi-xc)^2 + (yi-yc)^2 - rc^2;
    % h=ezplot(f,[-1000,4000,0,4000]);
    % coor_cir = get(h,'contourMatrix');
    % set(h, 'Color', 'r');
    radius_cir= rc;
    dc= rc*2;
    save('dc', 'dc')
    
    %% to extract boundary points at constant degree intervals
    
    points=[];
    lineLength = graindata.MajorAxisLength(1,1)*0.8;
    bp=[];
    thint=0.0017; % 0.1 degree upto 4th decimal place
    for angle = 0:thint:2*pi
        x1 = graindata.Centroid(1,1);
        y1 = graindata.Centroid(1,2);
        x2 = x1 + lineLength * cos(angle);
        y2 = y1 + lineLength * sin(angle);
        xdp=[x1 x2];
        ydp=[y1 y2];
        if (abs(x1- x2)) > (abs(y1- y2))
            if x1 < x2
                Xdp= x1: x2;
            else
                Xdp= x1: -1: x2;
            end
            
            Xdp= Xdp.';
            Xdp= round(Xdp);
            Ydp=interp1(xdp,ydp,Xdp);
            Ydp= round(Ydp);
            
        else
            if y1 < y2
                Ydp= y1: y2;
            else
                Ydp= y1: -1: y2;
            end
            Ydp= Ydp.';
            Ydp= round(Ydp);
            Xdp=interp1(ydp,xdp,Ydp);
            Xdp= round(Xdp);
        end
        if isnan(Xdp(end,1))== 1
            Xdp(end,1)= Xdp(end-1,1);
        end
        if isnan(Ydp(end,1))== 1
            Ydp(end,1)= Ydp(end-1,1);
        end
        if isnan(Ydp(1,1))== 1
            Ydp(1,1)= Ydp(2,1);
        end
        if isnan(Xdp(1,1))== 1
            Xdp(1,1)= Xdp(2,1);
        end
        points=[Ydp Xdp];
        for i= 1: (size(points,1))
            p= points(i,2);
            q= points(i,1);
            if q> size(BW_filled,1) || p> size(BW_filled,2) || q<1 || p<1
                points(i,3) = 0;
            else
                ind = BW_filled(q,p);
                points(i,3) = ind;
            end
        end
        for i= 1: (size(points,1)-1)
            if points(i,3) ~= points(i+1,3)
                break
            end
        end
        for k=1: size(B{1},1)
            if points(i,2) == B{1}(k,2) && points(i,1) == B{1}(k,1)
                bp=[ bp; points(i,2) points(i,1)];
                break
            else
                if abs(points(i,2)-B{1}(k,2))<= 1 && abs(points(i,1)-B{1}(k,1))<= 1
                    bp=[ bp; B{1}(k,2) B{1}(k,1)];
                    break
                end
            end
        end
        points=[];
    end
    bp(end,:)= bp(1,:);
    % figure,plot(bp(:,1), bp(:,2), 'r', 'LineWidth', 2)
    
    %% to find polar coordinates
    radius=[];
    angle=[];
    for i=1:size(bp,1);
        x=graindata.Centroid(1,1)-bp(i,1);
        y=graindata.Centroid(1,2)-bp(i,2);
        r=sqrt(x^2+y^2);
        radius=[radius; r bp(i,2) bp(i,1)];
        the=atan(y/x);
        if (the <0 && x<0)
            the= 2*pi+the;
        end
        if the <0 && y<0 || (the>0 && x>0 && y>0)
            the= pi+the;
        end
        angle =[angle; the];
    end
    radius(end,1)= radius(1,1);
    ar=[angle radius];
    ar= sortrows(ar);
    ar(end,2)= ar(1,2);
    p1p=ar(1,2);
    ar(:,2)=ar(:,2)-p1p;
    radius= ar(:,2);
%     radius= [radius; radius(1:20,1)];
    save('radius', 'radius')
    angle= ar(:,1);
    N= size(radius,1);
    T= thint*N;
    angle1 = linspace(angle(1,1),T,N);
    
    %% to understand transmission characterics and filter specific wavelength
    
    % to plot the data in time domain
    
    n=N;% number of profile points
    dx = thint; % spacing in mm
    x = (0:1:n-1).'*dx; % generate x axis data from 0 to 8 mm
    % figure, plot(x,radius,'k');
    % xlabel('Distance (mm)');
    % ylabel('Amplitude in time domain (\mum)');
    
    % to plot the data in frequency domain
    
    X_mags=fft(radius); % perform FFT
    j=(2:1:floor(n/2)+1); % generate the wavelength array
    lambda=n*dx./(j-1)'; % generate the wavelength array
    amp=abs(X_mags(2:floor(n/2)+1,1));
    har = (1:1:n/2).';
    har= (1./har);
    % figure,
    % %plot(lambda,amp,'Color',clr(P,:),'LineWidth', 2); % plot half of the % FFT array
    % plot(har,amp,'Color','red','LineWidth', 2); % plot half of the % FFT array
    % % % hold on
    % % % legendInfo{P} = ['image' num2str(P)];
    % % % end
    % xlabel('(1/ Descriptor no.)');
    % ylabel('Scaled amplitude (pixels)');
    % title('Frequency spectrum');
    % legend(legendInfo)
    
    
    
    
    
    %% cutoff
    amp= amp*2/n;
    noise_cut = 0.1* max(abs(amp))/100;
    rough_cut= 1.2* max(abs(amp))/100;
    amp=amp*n/2;
    noise_cut_amp= noise_cut*n/2;
    rough_cut_amp= rough_cut *n/2;
    for i= size(har,1):-1:1
        if amp(i,1)> noise_cut_amp
            harc_n= har(i,1); % in pixels
            lambdac_n= lambda(i,1);
            break
        else if i ==1
                IDX = knnsearch(har,noise_cut_amp) ;
                harc_n= har(IDX,1);
                lambdac_n = lambda(IDX,1);
            end
        end
    end
    for i= size(har,1):-1:1
        if amp(i,1)> rough_cut_amp
            harc_r= har(i,1); % in pixels
            lambdac_r= lambda(i,1);
            break
        else if i == 1
                harc_r=harc_n;
                lambdac_r= lambdac_n;
            end
        end
    end
    harmonic_cutoff_noise = 1/harc_n;
    harmonic_cutoff_roughness = 1/harc_r;
    fid = fopen(strcat('results_', filename, '.txt'),'a');
    if fid<0
        error('Couldn''t open results.txt for writing!');
    end
    
    fprintf(fid, '%d\n', harmonic_cutoff_roughness);
    disp(harmonic_cutoff_roughness);
    
    
    
     %% FOR REMOVING NOISE
    % to plot weighing function and transmission characteriss of gaussian regression filter
    
    % weighing function of gaussian regression filter
    
    load('radius')
    lambdac= lambdac_n;
    dx=thint; % in mm
    x1=(-lambdac:dx:lambdac)';
    alpha=0.4697;
    S=(1/(alpha*lambdac)).*exp(-pi*(x1/(alpha*lambdac)).^2);
    % generate the Gaussian filter
    S=S/sum(S); % normalize to unit sum
    % plot(x1,S);
    % xlabel('Wavelength component of outline ');
    % ylabel('weighting function');
    
    % transmission characteriss
    
    m=size(S,1); % length of Gaussian filter
    l=n*thint; % length of a profile is assumed to be 8 mm
    n=l/dx; % number of profile points
    S=[zeros(n/2-floor(m/2),1); S; zeros(n/2-floor(m/2)-1,1)];
    % center the filter and zero-pad to 8000 long
    Sf=fft(S); % DFT of S
    j=(2:1:floor(n/2)+1)';% generate wavelength array for X axis of
    % transmission plot
    wave=n*dx./(j-1);
    % figure,
    % semilogx(wave,100*abs(Sf(2:floor(n/2)+1,1)));
    % xlabel('Wavelength component of outline ');
    % ylabel('Amplitude (pixels');
    
    
    % low pass filter specified wavelength using gaussian regression filter of second order for 50% amplitude attenuation at cutoff
    const = sqrt(log(2)/2/pi/pi);
    for k = 1:n
        p = (1:1:n)'; % For each position k (center of the filter), generate the filter over the entire profile length
        S = (1/sqrt(2*pi*pi)/const/lambdac).*exp(-0.5*((k-p)*dx/const/lambdac).^2);% generate weighting function
        S = S/sum(S); % normalize to unit sum  and three variables given in Eq. 9.9. Ignore the dx term in Eq. 9.9 because it is a constant.
        x= (k-p)*dx;
        M=[];
        Q=[];
        M(1,1) = sum(S.*x.^4);
        M(1,2) = sum(S.*x.^3);
        M(1,3) = sum(S.*x.^2);
        M(2,1) = M(1,2);
        M(2,2) = M(1,3);
        M(2,3) = sum(S.*x);
        M(3,1) = M(1,3);
        M(3,2) = M(2,3);
        M(3,3) = sum(S);
        Q(1,1) = sum(radius.*S.*x.^2);
        Q(2,1) = sum(radius.*S.*x);
        Q(3,1) = sum(radius.*S);
        P = inv(M) * Q; % determine array P containing the values of % three unknown
        A = P(1); B = P(2); C = P(3);
        rapprox(k,1) = C; % determine waviness profile value at that
        % location
    end
%        figure, plot(x,radius,'r','LineWidth', 1.5);
%         hold on
%         plot(x,rapprox,'k','LineWidth', 1)
%         xlabel('Angle (radians)');
%         ylabel('Radial distance (pixels)');
%         title('Original Profile and Noise removed Profile');
    xxx=x;
    
    % reconstruction
    radius= radius + p1p;
    rapprox= rapprox + p1p;
    angle1=angle1.';
    % radius= [radius; radius(1,1)];
    % rapprox= [rapprox; rapprox(1,1)];
    x_new=[];
    y_new=[];
    for i=1:size(angle1,1)
        x_new=[x_new; rapprox(i,1)*cos(angle1(i,1))];
        y_new=[y_new; rapprox(i,1)*sin(angle1(i,1))];
    end
    
    x_new= x_new+ graindata.Centroid(1,1);
    y_new= y_new+ graindata.Centroid(1,2);
    x_new(end,1)=x_new(1,1);
    y_new(end,1)=y_new(1,1);
    x_old=[];
    y_old=[];
    for i=1:size(angle1,1)
        x_old=[x_old; radius(i,1)*cos(angle1(i,1))];
        y_old=[y_old; radius(i,1)*sin(angle1(i,1))];
    end
    x_old= x_old+ graindata.Centroid(1,1);
    y_old= y_old+ graindata.Centroid(1,2);
    y_new(end,1)= (y_new(end-1,1)+y_new(end-2,1)+y_new(1,1)+y_new(2,1)+y_new(3,1))/5;
    x_new(end,1)= (x_new(end-1,1)+x_new(end-2,1)+x_new(1,1)+x_new(2,1)+x_new(3,1))/5;
    y_new(1,1)= (y_new(end-1,1)+y_new(end-2,1)+y_new(1,1)+y_new(2,1)+y_new(3,1))/5;
    x_new(1,1)= (x_new(end-1,1)+x_new(end-2,1)+x_new(1,1)+x_new(2,1)+x_new(3,1))/5;
    % figure,
    % plot(x_new(:,1), y_new(:,1), 'green', 'LineWidth', 2)
    % figure,imshow(bw)
    % hold on
    % plot(x_new(:,1), y_new(:,1), 'red', 'LineWidth', 2)
    % hold on
    % plot(x_old(:,1), y_old(:,1), 'green', 'LineWidth', 2)
    % % plot(x_v(:,1),y_v(:,1), 'k*');
    % % plot(x_p(:,1), y_p(:,1), 'r*');
    % title('original figure and noise removed figure');
    xy_n=[y_new x_new];
    
    %% FOR REMOVING ROUGHNESS
    % to plot weighing function and transmission characteriss of gaussian filter
    
    % weighing function of gaussian filter
    
    % providing lambda cutoff
    load('radius')
    rapprox=rapprox-p1p;
    lambdac= lambdac_r;
    dx=thint; % in mm
    x1=(-lambdac:dx:lambdac)';
    alpha=0.4697;
    S=(1/(alpha*lambdac)).*exp(-pi*(x1/(alpha*lambdac)).^2);
    % generate the Gaussian filter
    S=S/sum(S); % normalize to unit sum
    % plot(x1,S);
    % xlabel('Distance (mm)');
    % ylabel('weighting function');
    
    % transmission characteriss
    
    m=size(S,1); % length of Gaussian filter
    l=n*thint; % length of a profile
    n=l/dx; % number of profile points
    S=[zeros(n/2-floor(m/2),1); S; zeros(n/2-floor(m/2)-1,1)];
    % center the filter and zero-pad to 8000 long
    Sf=fft(S); % DFT of S
    j=(2:1:floor(n/2)+1)';% generate wavelength array for X axis of
    % transmission plot
    wave=n*dx./(j-1);
    % figure,
    % semilogx(wave,100*abs(Sf(2:floor(n/2)+1,1)));
    % xlabel('Wavelength (mm)');
    % ylabel('Amplitude (\mum)');
    
    % low pass filter specified wavelength using gaussian regression filter of second order for amplitude attenuation at cutoff)
    
    const = sqrt(log(2)/2/pi/pi);
    for k = 1:n
        p = (1:1:n)'; % For each position k (center of the filter), generate the filter over the entire profile length
        S = (1/sqrt(2*pi*pi)/const/lambdac).*exp(-0.5*((k-p)*dx/const/lambdac).^2);% generate weighting function
        S = S/sum(S); % normalize to unit sum  and three variables given in Eq. 9.9. Ignore the dx term in Eq. 9.9 because it is a constant.
        x= (k-p)*dx;
        M=[];
        Q=[];
        M(1,1) = sum(S.*x.^4);
        M(1,2) = sum(S.*x.^3);
        M(1,3) = sum(S.*x.^2);
        M(2,1) = M(1,2);
        M(2,2) = M(1,3);
        M(2,3) = sum(S.*x);
        M(3,1) = M(1,3);
        M(3,2) = M(2,3);
        M(3,3) = sum(S);
        Q(1,1) = sum(radius.*S.*x.^2);
        Q(2,1) = sum(radius.*S.*x);
        Q(3,1) = sum(radius.*S);
        P = inv(M) * Q; % determine array P containing the values of % three unknown
        A = P(1); B = P(2); C = P(3);
        rapprox_r(k,1) = C; % determine waviness profile value at that location
    end
  %% to plot the weighing function of gaussian regression filter
% x = (0:1:n-1)'*dx;
% figure,
% for k=1:48:462 % for one eighth of the profile length
% p=(1:1:n)'; % For each position k (center of the filter),
% S = (1/sqrt(2*pi*pi)/const/lambdac).*exp(-0.5*((k-p)*dx/const/lambdac).^2);% generate weighting function
% S = S/sum(S);
% plot(l(1:462),S(1:462));
% hold on;
% end
    %%
    
%         figure, plot(xxx,radius,'r','LineWidth', 1.5);
%         hold on
%         plot(xxx,rapprox_r,'g','LineWidth', 1)
%         xlabel('Angle (radians)');
%         ylabel('Radial distance (pixels)');
%         title('Original profile and roughness removed profile');
%         figure, plot(xxx,rapprox,'k','LineWidth', 1.5);
%         hold on
%         plot(xxx,rapprox_r,'g','LineWidth', 1)
%         xlabel('Angle (radians)');
%         ylabel('Radial distance (pixels)');
%         title('Noise removed profile and Roughness removed profile');
    % rough_mean= linspace(mean(abs(rapprox-rapprox_r)),mean(abs(rapprox-rapprox_r)),size(rapprox,1));
    % %  figure, plot(x, abs(rapprox-rapprox_r),'k',x, rough_mean,'r')
    % %  xlabel('Distance (mm)');
    % % ylabel('Amplitude of roughness (\mum)');
    % % title('roughness removed profile and its mean');
    
    % reconstruction
    radius= radius + p1p;
    rapprox_r= rapprox_r + p1p;
    
    % end point correction
    
%         y= [rapprox_r(end-5,1),rapprox_r(end-4,1),rapprox_r(end-3,1),rapprox_r(end-2,1),rapprox_r(end-1,1),rapprox_r(end,1),rapprox_r(1,1),rapprox_r(2,1),rapprox_r(3,1),rapprox_r(4,1),rapprox_r(5,1),rapprox_r(6,1)];
%         yy = smooth(y,'sgolay') ;
%         figure, plot(y)
%         hold on
%         plot(yy,'r')
%         rapprox_r(end-3: end,1)= yy(1:4,1);
%         rapprox_r(1: 3,1)= yy(5:7,1);
x = 1:16;
yy= [rapprox_r(end-8,1),rapprox_r(end-7,1),rapprox_r(end-6,1),rapprox_r(end-5,1),rapprox_r(end-4,1),rapprox_r(end-3,1),rapprox_r(end-2,1),rapprox_r(end-1,1),rapprox_r(end,1),rapprox_r(1,1),rapprox_r(2,1),rapprox_r(3,1),rapprox_r(4,1),rapprox_r(5,1),rapprox_r(6,1),rapprox_r(7,1)];
% figure,
% plot(x, yy, 'b*-', 'LineWidth', 2, 'MarkerSize', 15);
coeffs = polyfit(x, yy, 3);
% Get fitted values
fittedX = linspace(min(x), max(x), 16);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
% hold on;
% plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
%         rapprox_r(end-3: end,1)= fittedY(1,1:4);
%         rapprox_r(1: 3,1)= fittedY(1,5:7);

x = 1:16;
% y= [rapprox_r(end-5,1),rapprox_r(end-4,1),rapprox_r(end-3,1),rapprox_r(end-2,1),rapprox_r(end-1,1),rapprox_r(end,1),rapprox_r(1,1),rapprox_r(2,1),rapprox_r(3,1),rapprox_r(4,1),rapprox_r(5,1),rapprox_r(6,1)];
%     figure,
%     plot(x, fittedY, 'b*-', 'LineWidth', 2, 'MarkerSize', 15);
coeffs = polyfit(x, fittedY, 2);
% Get fitted values
fittedX = linspace(min(x), max(x), 16);
fittedY1 = polyval(coeffs, fittedX);
% Plot the fitted line
%     hold on;
%     plot(fittedX, fittedY1, 'r-', 'LineWidth', 3);
rapprox_r(end-8: end,1)= fittedY1(1,1:9);
rapprox_r(1:7,1)= fittedY1(1,10:16);
fittedY1=fittedY1.';
x= 1:9;
yy= [fittedY1(5,1),fittedY1(6,1),fittedY1(7,1),fittedY1(8,1),fittedY1(9,1),fittedY1(10,1),fittedY1(11,1),fittedY1(12,1),fittedY1(13,1)];
% figure,
% plot(x, yy, 'b*-', 'LineWidth', 2, 'MarkerSize', 15);
coeffs = polyfit(x, yy, 3);
% Get fitted values
fittedX = linspace(min(x), max(x), 9);
fittedY2 = polyval(coeffs, fittedX);
% Plot the fitted line
% hold on;
% plot(fittedX, fittedY2, 'r-', 'LineWidth', 3);
rapprox_r(end-4: end,1)= fittedY2(1,1:5);
rapprox_r(1:4,1)= fittedY2(1,6:9);

%     radius(end,1)= radius(1,1);
%     radius=[radius; radius(2:10,1)];
%      rapprox_r(end,1)= rapprox_r(1,1);
%     rapprox_r=[rapprox_r; rapprox_r(2:10,1)];
x_new=[];
y_new=[];
for i=1:size(angle1,1)
    x_new=[x_new; rapprox_r(i,1)*cos(angle1(i,1))];
    y_new=[y_new; rapprox_r(i,1)*sin(angle1(i,1))];
end
x_new= x_new+ graindata.Centroid(1,1);
y_new= y_new+ graindata.Centroid(1,2);
x_new(end,1)=x_new(1,1);
y_new(end,1)=y_new(1,1);
x_old=[];
y_old=[];
for i=1:size(angle1,1)
    x_old=[x_old; radius(i,1)*cos(angle1(i,1))];
    y_old=[y_old; radius(i,1)*sin(angle1(i,1))];
end
x_old= x_old+ graindata.Centroid(1,1);
y_old= y_old+ graindata.Centroid(1,2);
%     y_new(end,1)= (y_new(end-1,1)+y_new(end-2,1)+y_new(1,1)+y_new(2,1)+y_new(3,1))/5;
%     x_new(end,1)= (x_new(end-1,1)+x_new(end-2,1)+x_new(1,1)+x_new(2,1)+x_new(3,1))/5;
%     y_new(1,1)= (y_new(end-1,1)+y_new(end-2,1)+y_new(1,1)+y_new(2,1)+y_new(3,1))/5;
%     x_new(1,1)= (x_new(end-1,1)+x_new(end-2,1)+x_new(1,1)+x_new(2,1)+x_new(3,1))/5;
%     figure,imshow(bw)
%     hold on
%     plot(x_new(:,1), y_new(:,1), 'red', 'LineWidth', 2)
%     hold on
%     plot(x_old(:,1), y_old(:,1), 'green', 'LineWidth', 2)
% title('original figure and roughness removed figure');
xy_r=[y_new x_new];
save('newb','xy_r')
save('open_fig','rapprox_r')
save('angle1','angle1')
mean_dif= mean(abs(rapprox_r - radius));
Rq = sqrt(mean((abs(rapprox_r - radius)- mean_dif).^2));
save('Rq', 'Rq');
toc
    
    %% CORNER IDENTIFICATION
    clear all
    load('input_temp.mat');
    load('h')
    load('P')
    result=glob(h,filename);
    he = imread(result{P});
    level = graythresh(he);
    bw = im2bw(he,level);
    bw = bwareaopen(bw, 50);
    BW_filled = imfill(bw,'holes');
    im_final = bw;
    graindata = regionprops(im_final, 'Area','Centroid','MajorAxisLength','MinorAxisLength','ConvexHull','Orientation','Eccentricity'); % using regionprops all the information is extracted from the image and stored in graindata matrix.
    [B,L] = bwboundaries(im_final,'noholes');
    load('newb')
    boundary =  B{1};
    % figure, imshow(bw)
    % hold on
    % plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 1.5)
    % title('Smoothened Boundary')
    
    %% polygon approximation
    tic
    load('dc')
    n= 0.09*dc/100; % maximum cutoff
    x = xy_r(:,2); % x component of coordinate
    y = xy_r(:,1); % y component of coordinate
    g= pdist2(xy_r,xy_r,'euclidean');
    nm = max(g(:));
    save('nm', 'nm')
    [ro co] = find(g== nm);
    
    % d=[];
    % m=[];
    % dist_p=[];
    % dist_pl=[];
    % figure, imshow (bw)
    % hold on
    % impoint(gca, graindata.Centroid(1,1),graindata.Centroid(1,2))
    % % to find the maximum distance of boundary point from centroid
    % for i = 1:length(x) - 1
    %     d = [d; x(i) y(i) sqrt((graindata.Centroid(1,2) - x(i))^2 + (graindata.Centroid(1,2) - y(i))^2) i];% maximum distance from the centroied
    % end
    % [b,ind] = max(d(:,3));%Index of the max d
    % x1= d(ind,1);% x component corresponds to max d
    % y1=d(ind,2);% y component corresponds to max d
    % index= d(ind,4);%Index of the max d
    % impoint(gca, x1,y1)
    % figure, imshow(bw)
    % hold on
    % impoint(gca, x1,y1)
    % % to find the maximum distance of a boundary point from the obtained boundary point
    % for j = 1:length(x) - 1
    %     distance = sqrt((x1 - x(j))^2 + (y1 - y(j))^2);% distance from the max d to all points in boundary
    %     m = [m; x(j) y(j) distance j];% all points and corresponding distances
    % end
    % [c,ind1] = max(m(:,3));%Index of the max m
    % x2= m(ind1,1);%x component corresponds to max m
    % y2=m(ind1,2);%y component corresponds to max m
    % index1= m(ind1,4);%Index of the max m
    % figure, plot(xy_r(:,2),xy_r(:,1), 'm', 'LineWidth', 2)
    % plot([x1, x2],[y1, y2],  'blue', 'LineWidth', 2)
    % hold on
    % impoint(gca, x2,y2)
    % line([x1 x2],[y1 y2], 'Color','g','LineWidth',3)
    
    index = ro(1,1);
    index1 = ro(2,1);
    %     figure, imshow(bw)
    %     hold on
    %     plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 1.5)
    % h=impoint(gca, xy_r(ro(1,1),2),xy_r(ro(1,1),1));
    % setColor(h,'g');
    % h1=impoint(gca, xy_r(ro(2,1),2),xy_r(ro(2,1),1));
    % setColor(h1,'g');
    % line([xy_r(ro(1,1),2) xy_r(ro(2,1),2)],[xy_r(ro(1,1),1) xy_r(ro(2,1),1)], 'Color','b','LineWidth',2)
    c= sqrt((xy_r(ro(1,1),2)-xy_r(ro(2,1),2))^2 + (xy_r(ro(1,1),1)-xy_r(ro(2,1),1))^2);
    % the boundary is divided and operated on in 2 parts
    % for first half
    
    if max(index,index1) <= (length(x))/2 % index is location of the two extreme points in xy_r matrix
        v=[ (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0; (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0 ]; % v stores coordinates of the two extreme points
        im=[ min(index,index1); max(index,index1)]; % index of extreme points are stored in im
        dsp=c; % c is Maximum distance
        dsp1=c;
        sd=c;
        p=1;
        p1=1;
        im_new=[];
        v_new=[];
        xh=[];
        yh=[];
        indexh=[];
        dist_p=[];
        while dsp>n
            for i= 1: (size(im,1)-1)
                for j= im(i) : im(i+1) % to run for each segment
                    pt = [x(j),y(j),0];
                    a = v(i,:)-v(i+1,:);%Vector
                    b = pt-v(i+1,:);%Vector
                    dist_p = [dist_p; x(j) y(j) (norm(cross(a,b)) / norm(a)) j];
                end
                [dsp,indh] = max(dist_p(:,3));
                xh= [xh; dist_p(indh,1)];
                yh= [yh; dist_p(indh,2)];
                indexh= [indexh; dist_p(indh,4)];
                %             impoint(gca, dist_p(indh,1),dist_p(indh,2))
                dist_p=[];
            end
            
            for l= 1: ((size(im,1))+(size(im,1)-1))
                if mod(l,2) == 1
                    im_new(l,:)= im(p,:);
                    v_new(l,:) = v(p,:);
                    p=p+1;% to check how many times the loop is executed
                else
                    im_new(l,:)= indexh(p1,:);
                    v_new(l,1) = xh(p1,1);
                    v_new(l,2) = yh(p1,1);
                    v_new(l,3) = 0;
                    p1=p1+1;% to check how many times the loop is executed
                end
            end
            im= im_new;
            v= v_new;
            xh=[];
            yh=[];
            indexh=[];
            p=1;
            p1=1;
            if dsp>n
                im_new=[];
                v_new=[];
            else
                for i= 1: (size(im,1)-1)
                    for j= im(i) : im(i+1)
                        pt = [x(j),y(j),0];
                        a = v(i,:)-v(i+1,:);%Vector
                        b = pt-v(i+1,:);%Vector
                         dist_p = [dist_p; x(j) y(j) (norm(cross(a,b)) / norm(a)) j];
                    end
                    [dsp1,indh] = max(dist_p(:,3));
                    if dsp1 > n
                        xh= [xh; dist_p(indh,1)];
                        yh= [yh; dist_p(indh,2)];
                        indexh= [indexh; dist_p(indh,4)];
                    else
                        xh= [xh; 0];
                        yh= [yh;0];
                        indexh= [indexh; 0];
                    end
                    dist_p=[];
                end
              for l= 1: ((size(im,1))+(size(im,1)-1))
                        
                 if mod(l,2) == 1
                    im_new(l,:)= im(p,:);
                    v_new(l,:) = v(p,:);
                    p=p+1;% to check how many times the loop is executed
                else
                    im_new(l,:)= indexh(p1,:);
                    v_new(l,1) = xh(p1,1);
                    v_new(l,2) = yh(p1,1);
                    v_new(l,3) = 0;
                    p1=p1+1;% to check how many times the loop is executed
                end
              end
            end
        end
        %         hold on
        %         for t= 1: (size(v_new,1)-1)
        %             plot([v_new(t,1),v_new(t+1,1)],[v_new(t,2),v_new(t+1,2)],'--gs','LineWidth',2,...
        %     'MarkerEdgeColor','b',...
        %     'MarkerFaceColor','y',...
        %     'MarkerSize',3)
        %             hold on
        %         end
        v_new1=[];
        im_new1=[];
        for i= 1: size(v_new,1)
            if im_new(i,1) ~=0
                v_new1=[v_new1; v_new(i,:)];
                im_new1=[im_new1; im_new(i,:)];
            end
        end
        v_new= v_new1;
        im_new= im_new1;
        
    else
        pp =[];
        ss= max(index,index1)- min(index,index1) +1;
        ssd= length(x)- max(index,index1);
        pp(1:ss,2)= x( min(index,index1): max(index,index1));
        pp(1:ss,1)= y( min(index,index1): max(index,index1));
        pp(ss+1:ss+ssd,2)= x( max(index,index1)+1: length(x));
        pp(ss+1:ss+ssd,1)= y( max(index,index1)+1: length(x));
        pp(ss+ssd+1:length(x),2)= x( 1: min(index,index1)-1);
        pp(ss+ssd+1:length(x),1)= y( 1: min(index,index1)-1);
        v=[ (pp(1,2)),(pp(1,1)),0; (pp(ss,2)),(pp(ss,1)),0 ];
        im=[ 1; ss];
        dsp=c;
        dsp1=c;
        sd=c;
        p=1;
        p1=1;
        im_new=[];
        v_new=[];
        xh=[];
        yh=[];
        indexh=[];
        dist_p=[];
        while dsp>n
            for i= 1: (size(im,1)-1)
                for j= im(i) : im(i+1)
                    pt = [pp(j,2),pp(j,1),0];
                    a = v(i,:)-v(i+1,:);%Vector
                    b = pt-v(i+1,:);%Vector
                    dist_p = [dist_p; pp(j,2) pp(j,1) (norm(cross(a,b)) / norm(a)) j];
                end
                [dsp,indh] = max(dist_p(:,3));
                xh= [xh; dist_p(indh,1)];
                yh= [yh; dist_p(indh,2)];
                indexh= [indexh; dist_p(indh,4)];
                %             h= impoint(gca, dist_p(indh,1),dist_p(indh,2))
                %             setColor(h,'m');
                dist_p=[];
            end
            
            for l= 1: ((size(im,1))+(size(im,1)-1))
                if mod(l,2) == 1
                    im_new(l,:)= im(p,:);
                    v_new(l,:) = v(p,:);
                    p=p+1;% to check how many times the loop is executed
                else
                    im_new(l,:)= indexh(p1,:);
                    v_new(l,1) = xh(p1,1);
                    v_new(l,2) = yh(p1,1);
                    v_new(l,3) = 0;
                    p1=p1+1;% to check how many times the loop is executed
                end
            end
            im= im_new;
            v= v_new;
            xh=[];
            yh=[];
            indexh=[];
            p=1;
            p1=1;
            if dsp> n
                im_new=[];
                v_new=[];
            else
                for i= 1: (size(im,1)-1)
                    for j= im(i) : im(i+1)
                        pt = [pp(j,2),pp(j,1),0];
                        a = v(i,:)-v(i+1,:);%Vector
                        b = pt-v(i+1,:);%Vector
                        dist_p = [dist_p; pp(j,2) pp(j,1) (norm(cross(a,b)) / norm(a)) j];
                    end
                    [dsp1,indh] = max(dist_p(:,3));
                    if dsp1 > n
                        xh= [xh; dist_p(indh,1)];
                        yh= [yh; dist_p(indh,2)];
                        indexh= [indexh; dist_p(indh,4)];
                    else
                        xh= [xh; 0];
                        yh= [yh;0];
                        indexh= [indexh; 0];
                    end
                    dist_p=[];
                end
              for l= 1: ((size(im,1))+(size(im,1)-1))
                        
                 if mod(l,2) == 1
                    im_new(l,:)= im(p,:);
                    v_new(l,:) = v(p,:);
                    p=p+1;% to check how many times the loop is executed
                else
                    im_new(l,:)= indexh(p1,:);
                    v_new(l,1) = xh(p1,1);
                    v_new(l,2) = yh(p1,1);
                    v_new(l,3) = 0;
                    p1=p1+1;% to check how many times the loop is executed
                end
              end
            end
        end
        
        %         hold on
        %         for t= 1: (size(v_new,1)-1)
        %             plot([v_new(t,1),v_new(t+1,1)],[v_new(t,2),v_new(t+1,2)],'--gs','LineWidth',2,...
        %     'MarkerEdgeColor','b',...
        %     'MarkerFaceColor','y',...
        %     'MarkerSize',3)
        %             hold on
        %         end
        v_new1=[];
        im_new1=[];
        for i= 1: size(v_new,1)
            if im_new(i,1) ~=0
                v_new1=[v_new1; v_new(i,:)];
                im_new1=[im_new1; im_new(i,:)];
            end
        end
        v_new= v_new1;
        im_new= im_new1;
    end
    
    % for 2nd half
    
    if max(index,index1) <= (length(x))/2
        vn=[ (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0; (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0 ];
        imn=[ max(index,index1); length(x)];
        dspn=c;
        dspn1=c;
        sd=c;
        pn=1;
        p1n=1;
        im_newn=[];
        v_newn=[];
        xhn=[];
        yhn=[];
        indexhn=[];
        dist_pn=[];
        while dspn >n
            for in= 1: (size(imn,1)-1)
                if imn(in) > 0
                    for jn= imn(in) : imn(in+1)
                        if imn(in+1) > 0
                            ptn = [x(jn),y(jn),0];
                            an = vn(in,:) - vn(in+1,:);
                            bn = ptn - vn(in+1,:);
                            dist_pn = [dist_pn; x(jn) y(jn) (norm(cross(an,bn)) / norm(an)) jn];
                        end
                    end
                    [dspn,indh] = max(dist_pn(:,3));
                    xhn= [xhn; dist_pn(indh,1)];
                    yhn= [yhn; dist_pn(indh,2)];
                    indexhn= [indexhn; dist_pn(indh,4)];
                    %                  impoint(gca, dist_pn(indh,1),dist_pn(indh,2))
                    dist_pn=[];
                end
            end
            
            for ln= 1: ((size(imn,1))+(size(imn,1)-1))
                if mod(ln,2) == 1
                    im_newn(ln,:) = imn(pn,:);
                    v_newn(ln,:) = vn(pn,:);
                    pn=pn+1;
                else
                    im_newn(ln,:) = indexhn( p1n,:);
                    v_newn(ln,1) = xhn( p1n,1);
                    v_newn(ln,2) = yhn( p1n,1);
                    v_newn(ln,3) = 0;
                    p1n=p1n+1;
                end
            end
            imn= im_newn;
            vn= v_newn;
            xhn=[];
            yhn=[];
            indexhn=[];
            pn=1;
            p1n=1;
            if dspn> n
                im_newn=[];
                v_newn=[];
            else
                for in= 1: (size(imn,1)-1)
                    if imn(in) > 0
                        for jn= imn(in) : imn(in+1)
                            if imn(in+1) > 0
                                ptn = [x(jn),y(jn),0];
                                an = vn(in,:) - vn(in+1,:);
                                bn = ptn - vn(in+1,:);
                                dist_pn = [dist_pn; x(jn) y(jn) (norm(cross(an,bn)) / norm(an)) jn];
                            end
                        end
                        [dspn1,indh] = max(dist_pn(:,3));
                        if dspn1 > n
                            xhn= [xhn; dist_pn(indh,1)];
                            yhn= [yhn; dist_pn(indh,2)];
                            indexhn= [indexhn; dist_pn(indh,4)];
                        else
                            xhn= [xhn; 0];
                            yhn= [yhn; 0];
                            indexhn= [indexhn; 0];
                        end
                        dist_pn=[];
                    end
                end
                     for ln= 1: ((size(imn,1))+(size(imn,1)-1))
                if mod(ln,2) == 1
                    im_newn(ln,:) = imn(pn,:);
                    v_newn(ln,:) = vn(pn,:);
                    pn=pn+1;
                else
                    im_newn(ln,:) = indexhn( p1n,:);
                    v_newn(ln,1) = xhn( p1n,1);
                    v_newn(ln,2) = yhn( p1n,1);
                    v_newn(ln,3) = 0;
                    p1n=p1n+1;
                end
            end
                
            end
            
        end
        %         for tn= 1: (size(v_newn,1)-1)
        %             plot([v_newn(tn,1),v_newn(tn+1,1)],[v_newn(tn,2),v_newn(tn+1,2)],'--gs','LineWidth',2,...
        %     'MarkerEdgeColor','b',...
        %     'MarkerFaceColor','y',...
        %     'MarkerSize',3)
        %             hold on
        %         end
        %         hold on
         v_newn1=[];
        im_newn1=[];
        for i= 1: size(v_newn,1)
            if im_newn(i,1) ~=0
                v_newn1=[v_newn1; v_newn(i,:)];
                im_newn1=[im_newn1; im_newn(i,:)];
            end
        end
        v_newn= v_newn1;
        im_newn= im_newn1;
       
    else
        pp =[];
        ss= max(index,index1)- min(index,index1) +1;
        ssd= length(x)- max(index,index1);
        pp(1:ss,2)= x( min(index,index1): max(index,index1));
        pp(1:ss,1)= y( min(index,index1): max(index,index1));
        pp(ss+1:ss+ssd,2)= x( max(index,index1)+1: length(x));
        pp(ss+1:ss+ssd,1)= y( max(index,index1)+1: length(x));
        pp(ss+ssd+1:length(x),2)= x( 1: min(index,index1)-1);
        pp(ss+ssd+1:length(x),1)= y( 1: min(index,index1)-1);
        vn=[ (pp(ss,2)),(pp(ss,1)),0; (pp(1,2)),(pp(1,1)),0 ];
        imn=[ ss+1; length(pp)];
        dspn=c;
        dspn1=c;
        pn=1;
        p1n=1;
        im_newn=[];
        v_newn=[];
        xhn=[];
        yhn=[];
        indexhn=[];
        dist_pn=[];
        while dspn>n
            for in=1: (size(imn,1)-1)
                if imn(in) > 0
                    for jn= imn(in) : imn(in+1)
                        if imn(in+1) > 0
                            ptn = [pp(jn,2),pp(jn,1),0];
                            an = vn(in,:) - vn(in+1,:);
                            bn = ptn - vn(in+1,:);
                            dist_pn = [dist_pn; pp(jn,2) pp(jn,1) (norm(cross(an,bn)) / norm(an)) jn];
                        end
                    end
                    [dspn,indh] = max(dist_pn(:,3));
                    xhn= [xhn; dist_pn(indh,1)];
                    yhn= [yhn; dist_pn(indh,2)];
                    indexhn= [indexhn; dist_pn(indh,4)];
                    %                 h= impoint(gca, dist_pn(indh,1),dist_pn(indh,2))
                    %                 setColor(h, 'm')
                    dist_pn=[];
                end
            end
            
            for ln= 1: ((size(imn,1))+(size(imn,1)-1))
                if mod(ln,2) == 1
                    im_newn(ln,:) = imn(pn,:);
                    v_newn(ln,:) = vn(pn,:);
                    pn=pn+1;
                else
                    im_newn(ln,:) = indexhn( p1n,:);
                    v_newn(ln,1) = xhn( p1n,1);
                    v_newn(ln,2) = yhn( p1n,1);
                    v_newn(ln,3) = 0;
                    p1n=p1n+1;
                end
            end
            imn= im_newn;
            vn= v_newn;
            xhn=[];
            yhn=[];
            indexhn=[];
            pn=1;
            p1n=1;
            if dspn> n
                im_newn=[];
                v_newn=[];
               else
                for in= 1: (size(imn,1)-1)
                    if imn(in) > 0
                        for jn= imn(in) : imn(in+1)
                            if imn(in+1) > 0
                                ptn =  [pp(jn,2),pp(jn,1),0];
                                an = vn(in,:) - vn(in+1,:);
                                bn = ptn - vn(in+1,:);
                                dist_pn = [dist_pn; pp(jn,2) pp(jn,1) (norm(cross(an,bn)) / norm(an)) jn];
                            end
                        end
                        [dspn1,indh] = max(dist_pn(:,3));
                        if dspn1 > n
                            xhn= [xhn; dist_pn(indh,1)];
                            yhn= [yhn; dist_pn(indh,2)];
                            indexhn= [indexhn; dist_pn(indh,4)];
                        else
                            xhn= [xhn; 0];
                            yhn= [yhn; 0];
                            indexhn= [indexhn; 0];
                        end
                        dist_pn=[];
                    end
                end
             for ln= 1: ((size(imn,1))+(size(imn,1)-1))
                if mod(ln,2) == 1
                    im_newn(ln,:) = imn(pn,:);
                    v_newn(ln,:) = vn(pn,:);
                    pn=pn+1;
                else
                    im_newn(ln,:) = indexhn( p1n,:);
                    v_newn(ln,1) = xhn( p1n,1);
                    v_newn(ln,2) = yhn( p1n,1);
                    v_newn(ln,3) = 0;
                    p1n=p1n+1;
                end
            end
                
            end
            
        end
        %         for tn= 1: (size(v_newn,1)-1)
        %             plot([v_newn(tn,1),v_newn(tn+1,1)],[v_newn(tn,2),v_newn(tn+1,2)], '--gs','LineWidth',2,...
        %     'MarkerEdgeColor','b',...
        %     'MarkerFaceColor','y',...
        %     'MarkerSize',3)
        %             hold on
        %         end
        %         hold on
            v_newn1=[];
        im_newn1=[];
        for i= 1: size(v_newn,1)
            if im_newn(i,1) ~=0
                v_newn1=[v_newn1; v_newn(i,:)];
                im_newn1=[im_newn1; im_newn(i,:)];
            end
        end
        v_newn= v_newn1;
        im_newn= im_newn1;
    end
    v_new(:,3) = 0;
     v_newn(:,3) = 0;
    INDEX1= [];
    INDEX1= [im_new; im_newn];
    COOR_DP1 = [];
    COOR_DP1=[v_new ; v_newn];
    COOR_DP=[];
    INDEX= [];
    po=1;
    for i= 1:size(COOR_DP1,1)
        if i== size(COOR_DP1,1)
            COOR_DP(po,1:2)= COOR_DP1(i,1:2);
            INDEX(po,:)= INDEX1(i,:);
            break
        end
        if COOR_DP1(i,1) ~= COOR_DP1(i+1,1) || COOR_DP1(i,2) ~= COOR_DP1(i+1,2)
            COOR_DP(po,1:2)= COOR_DP1(i,1:2);
            INDEX(po,:)= INDEX1(i,:);
            po= po+1;
        end
    end
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 2)
%     hold on
%     plot(COOR_DP(:,1),COOR_DP(:,2),'--gs','LineWidth',2,...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor','y',...
%         'MarkerSize',3)
    
    %% to separate concave inwards and concave outwards portion
    
    % to find points on the joining line segments
    Xdp1=[];% dp stands for dominant points
    Xdp2=[];
    Ydp1=[];
    Ydp2=[];
    COOR_LP=[];
    gh = 1;
    
    for i= 1: (size(INDEX,1)-1)
        
        a1=COOR_DP(i,:); % coordinates of dominant points
        b1=COOR_DP(i+1,:);
        % Xdp & Ydp stores points of adjoining line segments joining 2 dominant points
        xdp=[a1(2) b1(2)];
        ydp=[a1(1) b1(1)];
        if (abs(a1(2)- b1(2))) > (abs(a1(1)- b1(1)))
            if a1(2) < b1(2)
                Xdp= a1(2): gh: b1(2);
            else
                Xdp= a1(2): -gh: b1(2); % Xdp stores poi
            end
            
            Xdp= Xdp.'; % transposed
            
            Xdp1=[Xdp1; Xdp]; % Xdp1 stores coordinates of all line segments
            Xdp2=[Xdp2; 0; Xdp];% Xdp2 stores coordinates of all line segments separated by zero
            Ydp=interp1(xdp,ydp,Xdp);% Ydp interpolates y coordinates of line segments
            %         Ydp= round(Ydp);
            Ydp1=[ Ydp1; Ydp];
            
            Ydp2=[ Ydp2; 0; Ydp];
        else
            if a1(1) < b1(1)
                Ydp= a1(1): gh:  b1(1);
            else
                Ydp= a1(1): -gh: b1(1);
            end
            %Ydp=min(a1(1),b1(1)):max(a1(1),b1(1));
            Ydp= Ydp.';
            Ydp1=[Ydp1; Ydp];
            Ydp2=[ Ydp2; 0; Ydp];
            Xdp=interp1(ydp,xdp,Ydp); % interpolation
            %         Xdp= round(Xdp);
            Xdp1=[ Xdp1; Xdp];
            
            Xdp2=[ Xdp2; 0; Xdp];
        end
    end
    COOR_LP=[Ydp2 Xdp2]; % all coordinates of line segments between dominant points is stored
    % index of dominant point coordinates is stored
    for i= 1:(size(COOR_DP,1))
        
        COOR_DP(i,3)= i;
    end
%     figure, imshow(bw)
%     hold on
%     plot(COOR_DP(:,1),COOR_DP(:,2),'--rs','LineWidth',2,...
%         'MarkerEdgeColor','y',...
%         'MarkerFaceColor','g',...
%         'MarkerSize',10)
%     
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'm', 'LineWidth', 2)
    % hold on
    
    % to find intersecting points of lines from white index to black index incase the line segment joining 2 dominant points crosses the boundary
    % to segregate a convex plus concave part into convex and concave part % separately for calculation
    % to find the corresponding indicies from BW_filled_n matrix of all points in the line segment joining dominant points
    
    % stored in 3rd column of COOR_LP
    ind=0;
    for i= 1: (size(COOR_LP,1))
        p= COOR_LP(i,2);
        q= COOR_LP(i,1);
        if p ==0 && q == 0 % if coordinate is 0, then index stored as 2
            ind =2;
            COOR_LP(i,3) = ind;
        else
            ind = inpolygon(COOR_LP(i,1),COOR_LP(i,2),xy_r(:,2),xy_r(:,1));
            COOR_LP(i,3) = ind;
        end
    end
    cn=0;
    % COOR_DP has intersecting points stored and indexed
    for i= 1: (size(COOR_LP,1)-1)
        if COOR_LP(i,3)==2
            cn= cn+1;
        else
            if COOR_LP(i,3) ~= COOR_LP(i+1,3) && COOR_LP(i+1,3)~=2 % if adjacent coordinate index is not same then store the last white index
                if COOR_LP(i,3) == 1;
                    
                    COOR_DP = [COOR_DP ; COOR_LP(i,1) COOR_LP(i,2) (COOR_DP(cn,3)+COOR_DP(cn+1,3))/2];% the third column is added to sort the rows in COOR_DP
                else
                    COOR_DP = [COOR_DP ; COOR_LP(i+1,1) COOR_LP(i+1,2) (COOR_DP(cn,3)+COOR_DP(cn+1,3))/2];
                end
            end
        end
    end
    c= COOR_DP(130:end,:);
    COOR_DP= sortrows(COOR_DP, 3); % to arrange the new set of dominant points including intersecting points
    COOR_DPF=[]; % removes repeated points in COOR_DP
    pf=1;
    for i= 1:size(COOR_DP,1)
        if i== size(COOR_DP,1)
            COOR_DPF(pf,1:2)= COOR_DP(i,1:2);
            break
        end
        if COOR_DP(i,1) ~= COOR_DP(i+1,1) || COOR_DP(i,2) ~= COOR_DP(i+1,2)
            COOR_DPF(pf,1:2)= COOR_DP(i,1:2);
            pf= pf+1;
        end
        
    end
    COOR_DPF = COOR_DPF(setdiff(1:size(COOR_DPF,1),[size(COOR_DPF,1)]),:); % to remove the last coordinate matching with the first one
    
    % to store position of coordinates in the matrix
    for i= 1:(size(COOR_DPF,1))
        
        COOR_DPF(i,3)= i;
    end
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 2)
%     hold on
%     plot(COOR_DPF(:,1),COOR_DPF(:,2),'--gs','LineWidth',2,...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor','y',...
%         'MarkerSize',3)
    % to get rid of closing error ignore te dominant points in that area
    COOR_DPPF=[];
    dis=[];
    for i= 1: size(COOR_DPF,1)
        dis=[dis; sqrt((COOR_DPF(i,1)- xy_r(1,2))^2 +(COOR_DPF(i,2)- xy_r(1,1))^2)];
    end
    [k,idx]= min(dis);
    COOR_DPPF= COOR_DPF(1: idx-2,:);
    COOR_DPPF=[COOR_DPPF; COOR_DPF(idx+3:end,:)];
    for i= 1:(size(COOR_DPPF,1))
        COOR_DPPF(i,3)= i;
    end
    COOR_DPFK=COOR_DPPF;
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 2)
%     hold on
%     plot(COOR_DPFK(:,1),COOR_DPFK(:,2),'--gs','LineWidth',2,...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor','y',...
%         'MarkerSize',3)
%     copy=[];
%     for i=1: size(COOR_DPF,1)
%         for j=1: size(COOR_DP1,1)
%             if j == size(COOR_DP1,1)
%                 copy=[copy; COOR_DPF(i,:)];
%             end
%             if COOR_DPF(i,1:2) ~= COOR_DP1(j,1:2)
%                 continue
%             end
%             if COOR_DPF(i,1:2) == COOR_DP1(j,1:2)
%                 break
%             end
%             
%         end
%     end
    % hold on
    % plot(copy(:,1),copy(:,2),'m*','linestyle','none',...
    %     'MarkerEdgeColor','m',...
    %     'MarkerFaceColor','y',...
    %     'MarkerSize',5)
    % title('Final Set of Dominants including intersecting points')
    
    %% to store potential corner regions
    
    kin=[];
    cp=0;
    hh=[];
    hh3=[];
    f=1;
    sf1=0;
    hh1=[];
    hh2=[];
    sf=0;
    aw=0;
    b1=1;
    cc=[];
    q1=[];
    ab=0;
    COOR_DPF1=[COOR_DPFK ; COOR_DPFK];
    for i= 1:(size(COOR_DPF1,1))
        
        COOR_DPF1(i,3)= i;
    end
    lim=[0 0 0]; % to store coordinates of corners
    CORNER=[];
    % to find smallest corner comprising 3 points
    for i= 1:(size(COOR_DPF1,1))-1
        
        for j= b1+1:(size(COOR_DPF1,1))
            
            if (COOR_DPF1(b1,1)== COOR_DPF1(j,1)) && (COOR_DPF1(b1,2)== COOR_DPF1(j,2))
                break
            end
            % to find adjacent points are inside the figure or not
            ap=COOR_DPF1(b1,1:2);
            bp=COOR_DPF1(j,1:2);
            
            xlp=[ap(2) bp(2)];
            ylp=[ap(1) bp(1)];
            if (abs(ap(2)- bp(2)))> (abs(ap(1)- bp(1)))
                if ap(2) < bp(2)
                    Xlp= ap(2): bp(2);
                else
                    Xlp= ap(2): -1: bp(2);
                end
                Xlp= Xlp.';
                Ylp=interp1(xlp,ylp,Xlp);
                %             Ylp= round(Ylp);
            else
                if ap(1) < bp(1)
                    Ylp= ap(1): bp(1);
                else
                    Ylp= ap(1): -1: bp(1);
                end
                
                Ylp= Ylp.';
                Xlp=interp1(ylp,xlp,Ylp);
                %             Xlp= round(Xlp);
            end
            kin=[]; % to store indicies of line segment joining adjacent dominant points
            for k= 1: size(Xlp,1)
                ind1 = inpolygon(Ylp(k,1), Xlp(k,1),xy_r(:,2),xy_r(:,1));
                kin= [kin ; ind1];% to store the indicies of line segment joining two dominant points
                
            end
            cp=0;
            for h= 1: size(kin,1)
                if kin(h,1)==1
                    cp=cp+1;
                end
            end
            % all white
            if ( (cp == (size(kin,1)))) || ((cp <= (size(kin,1)))&& (cp >= 3))
                aw=aw +1;
                hh1=[hh1; b1 j];
                hh2=[hh2; b1 j];
                
                lim=[lim; COOR_DPF1(b1,:);COOR_DPF1(j,:) ];
                break
            end
            % all black
            if cp == 1 || cp ==2
                ab=ab+1;
                break
            end
            
        end
        
        b1=j;
        if aw==1 && ab~=1
            
            continue
        end
        if (ab>0 && aw==0)
            lim=[lim; 0 0 0];
            ab=0;
            continue
        end
        if aw==1 && ab==1
            [r,c,v]= find(not(lim),1,'last');
            if r~=1
                [r1,c1,v1]= find(not(lim(1:r-1,1)),1,'last');
                if (r-r1)>= 4
                    ap=lim(r-4,1:2);
                    bp=lim(r+2,1:2);
                    xlp=[ap(2) bp(2)];
                    ylp=[ap(1) bp(1)];
                    if (abs(ap(2)- bp(2)))> (abs(ap(1)- bp(1)))
                        if ap(2) < bp(2)
                            Xlp= ap(2): bp(2);
                        else
                            Xlp= ap(2): -1: bp(2);
                        end
                        Xlp= Xlp.';
                        Ylp=interp1(xlp,ylp,Xlp);
                    else
                        if ap(1) < bp(1)
                            Ylp= ap(1): bp(1);
                        else
                            Ylp= ap(1): -1: bp(1);
                        end
                        
                        Ylp= Ylp.';
                        Xlp=interp1(ylp,xlp,Ylp);
                    end
                    kin=[];
                    for k= 1: size(Xlp,1)
                        ind1 = inpolygon(Ylp(k,1), Xlp(k,1),xy_r(:,2),xy_r(:,1));
                        kin= [kin ; ind1];
                        
                    end
                    cp=0;
                    for h= 1: size(kin,1)
                        if kin(h,1)==1
                            cp=cp+1;
                        end
                    end
                    
                    if ( (cp == (size(kin,1)))) && size(CORNER,1) > 0 
                        CORNER(end,:)= lim(r+2,:);
                    end
                    lim=[lim; 0 0 0];
                    aw=0;
                    ab=0;
                end
            else
                lim=[lim; 0 0 0 ];
            end
        end
        
        if aw==2
            [r,c,v]= find(not(lim),1,'last');
            lim=[lim; 0 0 0];
            ap=lim(r+1,1:2);
            bp=lim(end-1,1:2);
            xlp=[ap(2) bp(2)];
            ylp=[ap(1) bp(1)];
            if (abs(ap(2)- bp(2)))> (abs(ap(1)- bp(1)))
                if ap(2) < bp(2)
                    Xlp= ap(2): bp(2);
                else
                    Xlp= ap(2): -1: bp(2);
                end
                Xlp= Xlp.';
                Ylp=interp1(xlp,ylp,Xlp);
            else
                if ap(1) < bp(1)
                    Ylp= ap(1): bp(1);
                else
                    Ylp= ap(1): -1: bp(1);
                end
                
                Ylp= Ylp.';
                Xlp=interp1(ylp,xlp,Ylp);
            end
            kin=[];
            for k= 1: size(Xlp,1)
                ind1 = inpolygon(Ylp(k,1), Xlp(k,1),xy_r(:,2),xy_r(:,1));
                kin= [kin ; ind1];
                
            end
            cp=0;
            for h= 1: size(kin,1)
                if kin(h,1)==1
                    cp=cp+1;
                end
            end
            
            if ( (cp == (size(kin,1)))) && aw==2
                CORNER=[ CORNER ; lim(r+1,:); lim(end-1,:)];
            end
            aw=0;
            ab=0;
            
        end
        
        
        aw=0;
        ab=0;
        kin=[];
        bp=1;
        if size(CORNER,1) > 2 && CORNER(end,3)~= 0
            for k=1: (size(CORNER,1)-3)
                if CORNER(end,1) == CORNER(k,1) && CORNER(end,2) == CORNER(k,2)
                    if k~=1
                        bp= 2;
                        break
                    end
                end
            end
        end
        if size(CORNER,1) > 2 && CORNER(end,3)== 0
            for k=1: (size(CORNER,1)-2)
                if CORNER(end-1,1) == CORNER(k,1) && CORNER(end-1,2) == CORNER(k,2)
                    bp= 2;
                    break
                end
            end
        end
        if bp==2
            break
        end
        
    end
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 2)
%     hold on
%     plot(CORNER(:,1),CORNER(:,2),'--ys','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','r',...
%         'MarkerSize',4)
%     title('Potential Corner Regions')
    
    
    
    %% TO EXTRACT RIGHT CORNER FROM POTENTIAL CORNERS
    seg=[];
    bn=1;
    for i=1:size(CORNER,1)
        for j=bn: size(CORNER,1)
            % the non repeating adjacent values are stored in rep matrix and repeated values are stored in rep matrix
            if CORNER(j,3)~= CORNER(j+1,3)
                seg=[seg; CORNER(j,:) j];
                bn=bn+1;
                break
            end
            if CORNER(j,3)== CORNER(j+1,3)
                bn=bn+2;
                break
            end
        end
        if bn== size(CORNER,1)
            seg=[seg; CORNER(bn,:) bn];
        end
        if  bn >= size(CORNER,1)
            break
        end
    end
    CONVEX=seg(:,1:3);
    
    CONVEX=[];
    if size(seg,1)==2
        i=1;
        rep=[];
        % seg(i,4) stores the indicies of corner coordinates in CORNER matrix
        for j= seg(i,4): seg(i+1,4)
            if j== seg(i+1,4) && i == size(seg,1)-1 % to break the loop when 1 cycle is complete
                break
            end
            if CORNER (j,3)== CORNER (j+1,3)
                rep=[rep; j];
            end
        end
        rep=[rep; seg(i+1,4)];
        CONVEX= [ CONVEX; CORNER(seg(i,4),:)];
        for k=1:floor(size(rep,1)/10): size(rep,1)-1
            
            for l= k+1:floor(size(rep,1)/10): size(rep,1)
                if CORNER(rep(l-1,1),3)== CONVEX(end,3)
                    break
                end
                ap=CORNER(rep(k,1),1:2);
                bp=CORNER(rep(l,1),1:2);
                
                xlp=[ap(2) bp(2)];
                ylp=[ap(1) bp(1)];
                if (abs(ap(2)- bp(2)))> (abs(ap(1)- bp(1)))
                    if ap(2) < bp(2)
                        Xlp= ap(2): bp(2);
                    else
                        Xlp= ap(2): -1: bp(2);
                    end
                    Xlp= Xlp.';
                    Ylp=interp1(xlp,ylp,Xlp);
                else
                    if ap(1) < bp(1)
                        Ylp= ap(1): bp(1);
                    else
                        Ylp= ap(1): -1: bp(1);
                    end
                    
                    Ylp= Ylp.';
                    Xlp=interp1(ylp,xlp,Ylp);
                end
                kin=[];
                for m= 1: size(Xlp,1)
                    ind1 = inpolygon(Ylp(m,1), Xlp(m,1),xy_r(:,2),xy_r(:,1));
                    kin= [kin ; ind1];
                    
                end
                cp=0;
                for h= 1: size(kin,1)
                    if kin(h,1)==1
                        cp=cp+1;
                    end
                end
                if ( (cp == (size(kin,1))))
                    continue
                    
                    
                else
                    
                    break
                end
            end
            
            if ( (cp ~= (size(kin,1))))
                break
            end
        end
        CONVEX= [ CONVEX; CORNER(seg(i+1,4),:)];
        if ( (cp ~= (size(kin,1))))
            
            for i=1: size(seg,1)
                if mod (i,2) == 0 % because odd to even index is a corner pair
                    continue
                end
                rep=[];
                for j= seg(i,4): seg(i+1,4)
                    if j== seg(i+1,4) && i == size(seg,1)-1 % to break the loop when 1 cycle is complete
                        break
                    end
                    if CORNER (j,3)== CORNER (j+1,3)
                        rep=[rep; j];
                    end
                end
                rep=[rep; seg(i+1,4)];
                if size(rep,1)==1
                    CONVEX= [ CONVEX; CORNER(seg(i,4),:); CORNER(seg(i+1,4),:) ];
                else
                    CONVEX= [ CONVEX; CORNER(seg(i,4),:)];
                    for k=1: size(rep,1)-1
                        
                        for l= k+1: size(rep,1)
                            if CORNER(rep(l-1,1),3)== CONVEX(end,3)
                                break
                            end
                            ap=CORNER(rep(k,1),1:2);
                            bp=CORNER(rep(l,1),1:2);
                            
                            xlp=[ap(2) bp(2)];
                            ylp=[ap(1) bp(1)];
                            if (abs(ap(2)- bp(2)))> (abs(ap(1)- bp(1)))
                                if ap(2) < bp(2)
                                    Xlp= ap(2): bp(2);
                                else
                                    Xlp= ap(2): -1: bp(2);
                                end
                                Xlp= Xlp.';
                                Ylp=interp1(xlp,ylp,Xlp);
                            else
                                if ap(1) < bp(1)
                                    Ylp= ap(1): bp(1);
                                else
                                    Ylp= ap(1): -1: bp(1);
                                end
                                
                                Ylp= Ylp.';
                                Xlp=interp1(ylp,xlp,Ylp);
                            end
                            kin=[];
                            for m= 1: size(Xlp,1)
                                ind1 = inpolygon(Ylp(m,1), Xlp(m,1),xy_r(:,2),xy_r(:,1));
                                kin= [kin ; ind1];
                                
                            end
                            cp=0;
                            for h= 1: size(kin,1)
                                if kin(h,1)==1
                                    cp=cp+1;
                                end
                            end
                            if ( (cp == (size(kin,1))))
                                continue
                                
                                
                            else
                                CONVEX=[CONVEX; CORNER(rep(l-1,1),:);  CORNER(rep(l-1,1),:)];
                                break
                            end
                        end
                        
                        
                    end
                    CONVEX= [ CONVEX; CORNER(seg(i+1,4),:)];
                end
                
            end
        end
        
        
    else
        
        for i=1: size(seg,1)
            if mod (i,2) == 0 % because odd to even index is a corner pair
                continue
            end
            rep=[];
            for j= seg(i,4): seg(i+1,4)
                if j== seg(i+1,4) && i == size(seg,1)-1 % to break the loop when 1 cycle is complete
                    break
                end
                if CORNER (j,3)== CORNER (j+1,3)
                    rep=[rep; j];
                end
            end
            rep=[rep; seg(i+1,4)];
            if size(rep,1)==1
                CONVEX= [ CONVEX; CORNER(seg(i,4),:); CORNER(seg(i+1,4),:) ];
            else
                CONVEX= [ CONVEX; CORNER(seg(i,4),:)];
                for k=1: size(rep,1)-1
                    
                    for l= k+1: size(rep,1)
                        if CORNER(rep(l-1,1),3)== CONVEX(end,3)
                            break
                        end
                        ap=CORNER(rep(k,1),1:2);
                        bp=CORNER(rep(l,1),1:2);
                        
                        xlp=[ap(2) bp(2)];
                        ylp=[ap(1) bp(1)];
                        if (abs(ap(2)- bp(2)))> (abs(ap(1)- bp(1)))
                            if ap(2) < bp(2)
                                Xlp= ap(2): bp(2);
                            else
                                Xlp= ap(2): -1: bp(2);
                            end
                            Xlp= Xlp.';
                            Ylp=interp1(xlp,ylp,Xlp);
                        else
                            if ap(1) < bp(1)
                                Ylp= ap(1): bp(1);
                            else
                                Ylp= ap(1): -1: bp(1);
                            end
                            
                            Ylp= Ylp.';
                            Xlp=interp1(ylp,xlp,Ylp);
                        end
                        kin=[];
                        for m= 1: size(Xlp,1)
                            ind1 = inpolygon(Ylp(m,1), Xlp(m,1),xy_r(:,2),xy_r(:,1));
                            kin= [kin ; ind1];
                            
                        end
                        cp=0;
                        for h= 1: size(kin,1)
                            if kin(h,1)==1
                                cp=cp+1;
                            end
                        end
                        if ( (cp == (size(kin,1))))
                            continue
                            
                            
                        else
                            CONVEX=[CONVEX; CORNER(rep(l-1,1),:);  CORNER(rep(l-1,1),:)];
                            break
                        end
                    end
                    
                    
                end
                CONVEX= [ CONVEX; CORNER(seg(i+1,4),:)];
            end
            
        end
    end
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'r', 'LineWidth', 2)
%     hold on
%     plot(CONVEX(:,1),CONVEX(:,2),'--bs','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','y',...
%         'MarkerSize',4)
%     title('Corner Regions')
    %% to store index values of B matrix for corners and plot corner portions
    
    CORNER_PLOT=[];
    for i= 1: size(CONVEX,1)+1
        kk=0;
        if mod(i,2)~= 0 && i~=1
            if CONVEX(i-2,4)> CONVEX(i-1,4) && (CONVEX(i-1,3)-CONVEX(i-2,3))< (CONVEX(end,3)*0.5)
                CORNER_PLOT=[CORNER_PLOT; xy_r(CONVEX(i-2,4): size(xy_r), :); xy_r(1: CONVEX(i-1,4), :)];
            else if CONVEX(i-2,4)< CONVEX(i-1,4) && (CONVEX(i-1,3)-CONVEX(i-2,3))< (CONVEX(end,3)*0.5)
                    CORNER_PLOT=[CORNER_PLOT; xy_r(CONVEX(i-2,4): CONVEX(i-1,4) , :)];
                else if CONVEX(i-2,4)> CONVEX(i-1,4) && (CONVEX(i-1,3)-CONVEX(i-2,3))> (CONVEX(end,3)*0.5)
                        CORNER_PLOT=[CORNER_PLOT; xy_r(CONVEX(i-1,4): CONVEX(i-2,4) , :)];
                    end
                end
            end
        end
        if i <= size(CONVEX,1)
            for j= 1: size(xy_r,1)
                if (CONVEX(i,1)== xy_r(j,2) && CONVEX(i,2)== xy_r(j,1))
                    CONVEX(i,4)= j; % to store relevant index of xy_r matrix
                    kk=kk+1;
                    break
                end
            end
            dmin=[];
            if j== size(xy_r,1) && kk == 0
                for k= 1: size(xy_r,1)
                    dmin = [dmin; sqrt((CONVEX(i,1)-xy_r(k,2))^2 + (CONVEX(i,2)- xy_r(k,1))^2) k];
                end
                [idx,t]= min(dmin(:,1));
                CONVEX(i,4)= t;
            end
        else
            break
        end
        
    end
    FINAL_CORNER=[];
    if CONVEX(end,4)> CONVEX(1,4) && size(CONVEX,1)~=2
        FINAL_CORNER = CONVEX(3:end-1,:);
        FINAL_CORNER = [FINAL_CORNER; CONVEX(2,:)];
    else
        FINAL_CORNER= CONVEX;
    end
    
    FINAL_CORNER1=FINAL_CORNER;
    % % to remove straight line segments
    %
    % for i = 1: size(FINAL_CORNER,1)-1
    %     if mod(i,2)== 0
    %         continue
    %     end
    %     if FINAL_CORNER(i,3)< FINAL_CORNER(i+1,3)
    %         for j= FINAL_CORNER(i,3) : FINAL_CORNER(i+1,3) % to run for each segment
    %             pt = [COOR_DPF1(j,1),COOR_DPF1(j,2),0];
    %             s= [FINAL_CORNER(i,1:2),0];
    %             p= [FINAL_CORNER(i+1,1:2),0];
    %             a = s - p;%Vector
    %             b = pt-p;%Vector
    %             dist_p = [dist_p; (norm(cross(a,b)) / norm(a))];
    %         end
    %         if max(dist_p)> 0.01
    %             FINAL_CORNER1=[ FINAL_CORNER1; FINAL_CORNER(i,:); FINAL_CORNER(i+1,:)];
    %         end
    %         dist_p=[];
    %     else
    %         for j= FINAL_CORNER(i,3) : size(COOR_DPF1,1)/2 % to run for each segment
    %             pt = [COOR_DPF1(j,1),COOR_DPF1(j,2),0];
    %             s= [FINAL_CORNER(i,1:2),0];
    %             p= [FINAL_CORNER(i+1,1:2),0];
    %             a = s - p;%Vector
    %             b = pt-p;%Vector
    %             dist_p = [dist_p; (norm(cross(a,b)) / norm(a))];
    %         end
    %         for j= 1 : FINAL_CORNER(i+1,3) % to run for each segment
    %             pt = [COOR_DPF1(j,1),COOR_DPF1(j,2),0];
    %             s= [FINAL_CORNER(i,1:2),0];
    %             p= [FINAL_CORNER(i+1,1:2),0];
    %             a = s - p;%Vector
    %             b = pt-p;%Vector
    %             dist_p = [dist_p; (norm(cross(a,b)) / norm(a))];
    %         end
    %         if max(dist_p)> 0.01
    %             FINAL_CORNER1=[ FINAL_CORNER1; FINAL_CORNER(i,:); FINAL_CORNER(i+1,:)];
    %         end
    %         dist_p=[];
    %     end
    % end
    
    % to get correct corner plots
    
    CORNER_PLOT=[];
    for i=1: size(FINAL_CORNER1,1)
        kk=0;
        if mod(i,2)~= 0
            if FINAL_CORNER1(i,4)> FINAL_CORNER1(i+1,4)
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER1(i,4): size(xy_r), :); xy_r(1: FINAL_CORNER1(i+1,4), :)];
                
            else
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER1(i,4): FINAL_CORNER1(i+1,4) , :)];
                
            end
        end
    end
    toc
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'm', 'LineWidth', 2)
%     hold on
%     
%     plot(COOR_DPFK(:,1),COOR_DPFK(:,2),'--ys','LineWidth',2,...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',5)
%     hold on
%     
%     plot(FINAL_CORNER1(:,1),FINAL_CORNER1(:,2),'--gs','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
%     plot(xy_r(:,2),xy_r(:,1), 'm', 'LineWidth', 2)
%     hold on
%     plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'r*','LineWidth',2,...
%         'MarkerEdgeColor','r',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
    save('CORN','FINAL_CORNER1')
    save('DOMI','COOR_DPF1');
    save('CORN_PL','CORNER_PLOT');
    %% ROUNDNESS BY MODIFIED DOUBLE DERIVATIVE FORMULA
    clear all
    load('input_temp.mat');
    load('h')
    load('P')
    result=glob(h,filename);
    he = imread(result{P});% to read the current image
    level = graythresh(he);
    bw = im2bw(he,level);
    bw = bwareaopen(bw, 50);
    BW_filled = imfill(bw,'holes');
    im_final = bw;
    graindata = regionprops(im_final, 'Area','Perimeter','Centroid','MajorAxisLength','MinorAxisLength','ConvexHull','Orientation','Eccentricity'); % using regionprops all the information is extracted from the image and stored in graindata matrix.
    [B,L] = bwboundaries(im_final,'noholes');
    boundary =  B{1};
    
    %% to calculate slope by modified double derivative formula
    load('newb.mat')
    load('CORN.mat')
    load('CORN_PL.mat')
    load('DOMI.mat')
    load('open_fig.mat');
    kk=0;
    R_C=[];
    Ndd=[];
    Ddd=[];
    alpha_x=[];
    y1p=0;
    beta_y=[];
    slope1=[];
    slope_s=[];
    slope_diff=[];
    tic
    for j= 1:(size(FINAL_CORNER1,1)-1)
        if mod(j,2)== 0
            continue
        end
        Ndd=[];
        Ddd=[];
        y1p=0;
        
        if FINAL_CORNER1(j,4) < FINAL_CORNER1(j+1,4)
            k=1;
            p=1;
            xs=xy_r(FINAL_CORNER1(j,4),2);
            ys=xy_r(FINAL_CORNER1(j,4),1);
            ym=xy_r(k+FINAL_CORNER1(j,4),1);
            xm=xy_r(k+FINAL_CORNER1(j,4),2);
            for i=FINAL_CORNER1(j,4):FINAL_CORNER1(j+1,4)
                p = p+1;
                if p*k+1<=(FINAL_CORNER1(j+1,4)-FINAL_CORNER1(j,4)+1);
                    ye=xy_r(p*k+FINAL_CORNER1(j,4),1);
                    xe=xy_r(p*k+FINAL_CORNER1(j,4),2);
                    ys1=(ym-ys)/(xm-xs);
                    ys2=(ye-ym)/(xe-xm);
                    y1=(ys1+ys2)/2;
                    yy=atand(y1);
                    h = xe-xm;
                    Nd=xm-xs;
                    Nay=Nd*(ye-ym)+h*(ys-ym);
                    Ndd=[Ndd Nay];
                    Nax=(Nd*h*(Nd+h))/2;
                    Ddd=[Ddd Nax];
                    y2=Nay/Nax;
                    slope1=[slope1; y2 xs ys xm ym xe ye];
                    R=((1+(y1)^2)^1.5)/y2;
                    %Center of curvature
                    alpha=xs-((y1*(1+(y1)^2))/y2);
                    alpha_x=[alpha_x; alpha FINAL_CORNER1(j,3) FINAL_CORNER1(j+1,3) ];
                    beta=ys+((1+(y1)^2)/y2);
                    beta_y=[beta_y; beta FINAL_CORNER1(j,3) FINAL_CORNER1(j+1,3)];
                    xs=xm;
                    ys=ym;
                    xm=xe;
                    ym=ye;
                    R_C=[R_C; R FINAL_CORNER1(j,3) FINAL_CORNER1(j+1,3)];
                end
                y1p=y1;
            end
            
        else
            k=1;
            p=1;
            pp= xy_r(FINAL_CORNER1(j,4): size(xy_r,1),:);
            pp1=xy_r(2: FINAL_CORNER1(j+1,4),:);
            ph=[pp; pp1];
            xs=ph(1,2);
            ys=ph(1,1);
            ym=ph(k+1,1);
            xm=ph(k+1,2);
            for i=1:size(ph,1)
                p = p+1;
                if p*k+1<=size(ph,1);
                    ye=ph(p*k+1,1);
                    xe=ph(p*k+1,2);
                    ys1=(ym-ys)/(xm-xs);
                    ys2=(ye-ym)/(xe-xm);
                    y1=(ys1+ys2)/2;
                    yy=atand(y1);
                    h = xe-xm;
                    Nd=xm-xs;
                    Nay=Nd*(ye-ym)+h*(ys-ym);
                    Ndd=[Ndd Nay];
                    Nax=(Nd*h*(Nd+h))/2;
                    Ddd=[Ddd Nax];
                    y2=Nay/Nax;
                    slope1=[slope1; y2 xs ys xm  ym xe ye];
                    R=((1+(y1)^2)^1.5)/y2;
                    %Center of curvature
                    alpha=xs-((y1*(1+(y1)^2))/y2);
                    alpha_x=[alpha_x; alpha FINAL_CORNER1(j,3) FINAL_CORNER1(j+1,3)];
                    beta=ys+((1+(y1)^2)/y2);
                    beta_y=[beta_y; beta FINAL_CORNER1(j,3) FINAL_CORNER1(j+1,3)];
                    xs=xm;
                    ys=ym;
                    xm=xe;
                    ym=ye;
                    R_C=[R_C; R FINAL_CORNER1(j,3) FINAL_CORNER1(j+1,3)];
                end
                y1p=y1;
            end
        end
    end
    inf=[];
    for i = 1: size(slope1,1)-4
        if (abs(slope1(i,1)-slope1(i+1,1)) < 0.0010 && abs(slope1(i,1)) < 0.0010)|| (abs(slope1(i,1)-slope1(i+1,1)) > 8000 && abs(slope1(i,1)) > 8000 && abs(slope1(i+1,1)-slope1(i+2,1)) > 8000 && abs(slope1(i+1,1)) > 8000 && abs(slope1(i+2,1)-slope1(i+3,1)) > 8000 && abs(slope1(i+2,1)) > 8000 && abs(slope1(i+3,1)-slope1(i+4,1)) > 8000 && abs(slope1(i+3,1)) > 8000 )
            inf=[inf; slope1(i,:) ];
        end
    end
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%     hold on
%     plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*',...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',3)
%     hold on
    
    if size(inf,1)>0
%                 hold on
%                 plot(inf(:, 2),inf(:, 3), 'b*')
%                 hold on
%                 plot(inf(:, 4),inf(:, 5), 'b*')
%                 hold on
%                 plot(inf(:, 6),inf(:, 7), 'b*')
         infi1=[];
        for i= 1: size(inf,1)
            infi1= [ infi1; inf(i,2) inf(i,3); inf(i,4) inf(i,5); inf(i,6) inf(i,7)];
        end
        
        %to remove duplicate values from infi1
        infi= unique(infi1(:,1:2), 'rows', 'stable');
        infi = [infi; infi(1,:)];
    end
    %     title('Inflexion Points')
    %% inflexion points removed
    if size(inf,1)>0
        for i= 1: size(FINAL_CORNER1,1)
            FINAL_CORNER1(i,5) = i;
        end
        p=1;
        for i = 1: size(FINAL_CORNER1,1)
            if mod(i,2) == 0
                continue
            end
            if FINAL_CORNER1(i,4) < FINAL_CORNER1(i+1,4)
                for j= FINAL_CORNER1(i,4) : FINAL_CORNER1(i+1,4)
                    if j == FINAL_CORNER1(i+1,4) && (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2)) && (xy_r(j-1,2) == infi(p-1,1) && xy_r(j-1,1) == infi(p-1,2))
                        FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        p=p+1;
                        break
                    else
                        if j == FINAL_CORNER1(i+1,4)
                            break
                        end
                    end
                    
                    if j == FINAL_CORNER1(i,4) && (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2))
                        FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                    end
                    
                    if (xy_r(j,2) ~= infi(p,1) || xy_r(j,1) ~= infi(p,2))
                        if (xy_r(j+1,2) == infi(p,1) && xy_r(j+1,1) == infi(p,2) && xy_r(j+1,2) ~= FINAL_CORNER1(i+1,4))
                            FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        end
                    end
                    
                    if (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2)) && size(infi,1)>= (p+1)
                        if (xy_r(j+1,1) ~= infi(p+1,2)) || (xy_r(j+1,2) ~= infi(p+1,1))
                            FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        end
                        p=p+1;
                    end
                    nb=0;
                    for k= size(FINAL_CORNER1,1):-1:2 % to check if odd number of points are inserted in a corner
                        if FINAL_CORNER1(k,5)== FINAL_CORNER1(k-1,5)
                            nb= nb+1;
                        end
                    end
                    if j == FINAL_CORNER1(i+1,4) && mod(nb,2)~=0
                        FINAL_CORNER1=[FINAL_CORNER1; FINAL_CORNER1(i+1,1) FINAL_CORNER1(i+1,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                    end
                    
                    
                end
            else
                
                for j=  FINAL_CORNER1(i,4): size(xy_r,1)
                    if j == size(xy_r,1) && (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2)) && (xy_r(j-1,2) == infi(p-1,1) && xy_r(j-1,1) == infi(p-1,2))
                        FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        p=p+1;
                        break
                    else
                        if j== size(xy_r,1)
                            break
                        end
                    end
                    if j == FINAL_CORNER1(i,4) && (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2))
                        FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                    end
                    if (xy_r(j,2) ~= infi(p,1) || xy_r(j,1) ~= infi(p,2))
                        if (xy_r(j+1,2) == infi(p,1) && xy_r(j+1,1) == infi(p,2))
                            FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        end
                    end
                    
                    if (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2)) && size(infi,1)>= (p+1)
                        if (xy_r(j+1,1) ~= infi(p+1,2)) || (xy_r(j+1,2) ~= infi(p+1,1))
                            FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        end
                        p=p+1;
                    end
                    
                end
                for j=  1: FINAL_CORNER1(i+1,4)
                    if j == FINAL_CORNER1(i+1,4) && (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2)) && (xy_r(j-1,2) == infi(p-1,1) && xy_r(j-1,1) == infi(p-1,2))
                        FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        p=p+1;
                        break
                    else
                        if j == FINAL_CORNER1(i+1,4)
                            break
                        end
                    end
                    if (xy_r(j,2) ~= infi(p,1) || xy_r(j,1) ~= infi(p,2))
                        if (xy_r(j+1,2) == infi(p,1) && xy_r(j+1,1) == infi(p,2) && xy_r(j+1,2) ~= FINAL_CORNER1(i+1,4))
                            FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        end
                    end
                    
                    if (xy_r(j,2) == infi(p,1) && xy_r(j,1) == infi(p,2)) && size(infi,1)>= (p+1)
                        if (xy_r(j+1,1) ~= infi(p+1,2)) || (xy_r(j+1,2) ~= infi(p+1,1)) || ((xy_r(j+1,2) == infi(p,1) && xy_r(j+1,1) == infi(p,2)) && j == FINAL_CORNER1(i,4))
                            FINAL_CORNER1=[FINAL_CORNER1; infi(p,1) infi(p,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                        end
                        p=p+1;
                    end
                    nb=0;
                    for k= size(FINAL_CORNER1,1):-1:2 % to check if odd number of points are inserted in a corner
                        if FINAL_CORNER1(k,5)== FINAL_CORNER1(k-1,5)
                            nb= nb+1;
                        end
                    end
                    if j == FINAL_CORNER1(i+1,4) && mod(nb,2)~=0
                        FINAL_CORNER1=[FINAL_CORNER1; FINAL_CORNER1(i+1,1) FINAL_CORNER1(i+1,2) 0 j (FINAL_CORNER1(i,5)+FINAL_CORNER1(i+1,5))/2];
                    end
                    
                end
            end
        end
        
        FINAL_CORNER1= sortrows(FINAL_CORNER1,5);
        FINAL_CORNER=[];
        for i=1: size(FINAL_CORNER1,1)-1
            if mod(i,2) == 0
                continue
            end
            if abs(FINAL_CORNER1(i,4)- FINAL_CORNER1(i+1,4)) > 4
                FINAL_CORNER=[FINAL_CORNER;  FINAL_CORNER1(i,:); FINAL_CORNER1(i+1,:) ];
            end
        end
        
        % inserting new corner ends in the dominant points
        COOR_D=[];
        for i=1: size(COOR_DPF1,1)
            COOR_D=[COOR_D; COOR_DPF1(i,1:2)];
        end
        
        for i= 1: size(FINAL_CORNER,1)
            if floor(FINAL_CORNER(i,5))~= FINAL_CORNER(i,5)
                Y= [FINAL_CORNER(i,1) FINAL_CORNER(i,2) ];
                IDX = knnsearch(COOR_D,Y) ;
                if (COOR_DPF1(IDX,1)< COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)> COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)< COOR_DPF1(IDX+1,2))
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX+1,3))/2];
                end
                if (COOR_DPF1(IDX,1)>COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)< COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)> COOR_DPF1(IDX+1,2) )
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX+1,3))/2];
                end
                if (COOR_DPF1(IDX,1)< COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)> COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)> COOR_DPF1(IDX+1,2) )
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX+1,3))/2];
                end
                if (COOR_DPF1(IDX,1)>COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)< COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)<COOR_DPF1(IDX+1,2))
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX+1,3))/2];
                end
                if (COOR_DPF1(IDX,1)< COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)< COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)< COOR_DPF1(IDX+1,2))&& IDX > 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX-1,3))/2];
                end
                if (COOR_DPF1(IDX,1)< COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)< COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)>COOR_DPF1(IDX+1,2))&& IDX > 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX-1,3))/2];
                end
                if (COOR_DPF1(IDX,1)> COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)> COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)< COOR_DPF1(IDX+1,2)) && IDX > 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX-1,3))/2];
                end
                if (COOR_DPF1(IDX,1)> COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)> COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)> COOR_DPF1(IDX+1,2))&& IDX > 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) (COOR_DPF1(IDX,3)+COOR_DPF1(IDX-1,3))/2];
                end
                  if (COOR_DPF1(IDX,1)< COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)< COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)< COOR_DPF1(IDX+1,2))&& IDX == 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) 0.5];
                end
                if (COOR_DPF1(IDX,1)< COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)< COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)>COOR_DPF1(IDX+1,2))&& IDX == 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) 0.5];
                end
                if (COOR_DPF1(IDX,1)> COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)> COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)< COOR_DPF1(IDX+1,2)) && IDX == 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) 0.5];
                end
                if (COOR_DPF1(IDX,1)> COOR_DPF1(IDX+1,1) && FINAL_CORNER(i,1)> COOR_DPF1(IDX,1)) && (COOR_DPF1(IDX,2)> COOR_DPF1(IDX+1,2))&& IDX == 1
                    COOR_DPF1=[COOR_DPF1; FINAL_CORNER(i,1) FINAL_CORNER(i,2) 0.5];
                end
            end
        end
        COOR_DPF1= sortrows(COOR_DPF1,3);
        COOR_D=[];
        for i=1: size(COOR_DPF1,1)
            COOR_D=[COOR_D; COOR_DPF1(i,1:2)];
        end
        
        % to input FINAL_CORNER(:,3)
        
        for i =1: size(FINAL_CORNER,1)
            Y= [FINAL_CORNER(i,1) FINAL_CORNER(i,2) ];
            IDX = knnsearch(COOR_D,Y) ;
            FINAL_CORNER(i,3) = IDX;
        end
    else
        FINAL_CORNER=  FINAL_CORNER1;
    end
    %% circles at all points in corner by revised double derivative formula
    kk=0;
    R_C=[];
    Ndd=[];
    Ddd=[];
    alpha_x=[];
    y1p=0;
    beta_y=[];
    slope=[];
    slope1=[];
    for j= 1:(size(FINAL_CORNER,1)-1)
        if mod(j,2)== 0
            continue
        end
        Ndd=[];
        Ddd=[];
        y1p=0;
        
        if FINAL_CORNER(j,4) < FINAL_CORNER(j+1,4)
            nop = abs(FINAL_CORNER(j+1,4)-FINAL_CORNER(j,4))+1;
            k=1;
            p=1;
            xs=xy_r(FINAL_CORNER(j,4),2);
            ys=xy_r(FINAL_CORNER(j,4),1);
            ym=xy_r(k+FINAL_CORNER(j,4),1);
            xm=xy_r(k+FINAL_CORNER(j,4),2);
            slope1=[slope1; 0 0 0 0 0 0 0 ];
            for i=FINAL_CORNER(j,4):FINAL_CORNER(j+1,4)
                p = p+1;
                if p*k+1<=(FINAL_CORNER(j+1,4)-FINAL_CORNER(j,4)+1);
                    ye=xy_r(p*k+FINAL_CORNER(j,4),1);
                    xe=xy_r(p*k+FINAL_CORNER(j,4),2);
                    ys1=(ym-ys)/(xm-xs);
                    ys2=(ye-ym)/(xe-xm);
                    y1=(ys1+ys2)/2;
                    yy=atand(y1);
                    slope=[slope; yy];
                    h = xe-xm;
                    Nd=xm-xs;
                    Nay=Nd*(ye-ym)+h*(ys-ym);
                    Ndd=[Ndd Nay];
                    Nax=(Nd*h*(Nd+h))/2;
                    Ddd=[Ddd Nax];
                    y2=Nay/Nax;
                    slope1=[slope1; y2 xe ye xm ym xs ys];
                    R=((1+(y1)^2)^1.5)/y2;
                    %Center of curvature
                    alpha=xs-((y1*(1+(y1)^2))/y2);
                    alpha_x=[alpha_x; alpha FINAL_CORNER(j,3) FINAL_CORNER(j+1,3) ];
                    beta=ys+((1+(y1)^2)/y2);
                    beta_y=[beta_y; beta FINAL_CORNER(j,3) FINAL_CORNER(j+1,3)];
                    xs=xm;
                    ys=ym;
                    xm=xe;
                    ym=ye;
                    R_C=[R_C; R FINAL_CORNER(j,3) FINAL_CORNER(j+1,3)];
                end
                y1p=y1;
            end
            
        else
            k=1;
            p=1;
            pp= xy_r(FINAL_CORNER(j,4): size(xy_r,1),:);
            pp1=xy_r(2: FINAL_CORNER(j+1,4),:);
            ph=[pp; pp1];
            if size(ph,1) > 3
            xs=ph(1,2);
            ys=ph(1,1);
            ym=ph(k+1,1);
            xm=ph(k+1,2);
            slope1=[slope1; 0 0 0 0 0 0 0 ];
            for i=1:size(ph,1)
                p = p+1;
                if p*k+1<=size(ph,1);
                    ye=ph(p*k+1,1);
                    xe=ph(p*k+1,2);
                    ys1=(ym-ys)/(xm-xs);
                    ys2=(ye-ym)/(xe-xm);
                    y1=(ys1+ys2)/2;
                    yy=atand(y1);
                    slope=[slope; yy];
                    h = xe-xm;
                    Nd=xm-xs;
                    Nay=Nd*(ye-ym)+h*(ys-ym);
                    Ndd=[Ndd Nay];
                    Nax=(Nd*h*(Nd+h))/2;
                    Ddd=[Ddd Nax];
                    y2=Nay/Nax;
                    slope1=[slope1; y2 xe ye xm ym xs ys];
                    R=((1+(y1)^2)^1.5)/y2;
                    %Center of curvature
                    alpha=xs-((y1*(1+(y1)^2))/y2);
                    alpha_x=[alpha_x; alpha FINAL_CORNER(j,3) FINAL_CORNER(j+1,3)];
                    beta=ys+((1+(y1)^2)/y2);
                    beta_y=[beta_y; beta FINAL_CORNER(j,3) FINAL_CORNER(j+1,3)];
                    xs=xm;
                    ys=ym;
                    xm=xe;
                    ym=ye;
                    R_C=[R_C; R FINAL_CORNER(j,3) FINAL_CORNER(j+1,3)];
                end
                y1p=y1;
            end
        end
        end
    end
    
    %% to plot new corners
    
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%     hold on
%     
    CORNER_PLOT=[];
    for i=1: size(FINAL_CORNER,1)
        kk=0;
        if mod(i,2)~= 0
            if FINAL_CORNER(i,4)> FINAL_CORNER(i+1,4)
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER(i,4): size(xy_r), :); xy_r(1: FINAL_CORNER(i+1,4), :)];
                
            else
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER(i,4): FINAL_CORNER(i+1,4) , :)];
                
            end
        end
    end
%     plot(FINAL_CORNER(:,1),FINAL_CORNER(:,2),'g*','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
%     plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',3,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     title('Final Corner from Corner Detection Algorithm')
    %% threshold radius, plotting and roundness by inscribed circle
    % inscibed circle
    ds=[];
    error1=[];
    error2=[];
    p=[];
    error11=[];
    for b= 1: size(BW_filled,1)
        for j= 1: size(BW_filled,2)
            if BW_filled(b,j)== 1
                ds(b,j)=0;
            end
            if BW_filled(b,j)== 0
                ds(b,j)=1;
            end
        end
    end
    D= bwdist(ds);
    [radius_in,ind] = max(D(:));
    [m,nn] = ind2sub(size(D),ind);
    % figure,
    % imagesc(D);
    % C = colormap;
    % L = size(C,1);
    % Gs = round(interp1(linspace(min(D(:)),max(D(:)),L),1:L,D));
    % H = reshape(C(Gs,:),[size(Gs) 3]);
    % figure,
    % subimage(mat2gray(D)), title('Euclidean')
    % hold on, imcontour(D)
    % hold on
    % plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
    % hold on
    % viscircles([nn,m],radius_in,'EdgeColor','w', 'LineWidth',1.5)
    % h=impoint(gca, m,nn);
    % setColor(h,'k');
    % colorbar
    % title('Maximum Inscribing Circle')
    save('m', 'm')
    save('nn', 'nn')
    save('radius_in', 'radius_in')
    
    %% to remove circles completely outside figure and more than inscribed circle
    rd_final=R_C;
    ccxe= alpha_x;
    ccye=beta_y;
    ccxe11=[];
    ccye11=[];
    rd_final11=[];
    
    for i=1:1:size(ccxe,1)
        if (inpolygon((ccxe(i,1)),(ccye(i,1)),xy_r(:,2),xy_r(:,1))) == 1
            ccxe11=[ ccxe11; ccxe(i,:)];
            ccye11=[ ccye11; ccye(i,:)];
            rd_final11=[rd_final11; rd_final(i,:)];
        end
    end
%     figure, imshow(bw)
%         hold on
%         plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%         hold on
%         plot(FINAL_CORNER(:,1),FINAL_CORNER(:,2),'g*','LineWidth',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','k',...
%             'MarkerSize',2)
%         hold on
%         plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','k',...
%             'MarkerSize',2)
%         hold on
%     for i=1:size(ccxe11,1)
%         syms x y
%         f(x,y)=(x-ccxe11(i,1))^2 + (y-ccye11(i,1))^2 - rd_final11(i,1)^2;
%         h=ezplot(f, [-1000,4000,0,4000]);
%         set(h, 'Color', 'g');
%         hold on
%     end
    
    ccxe10=[];
    ccye10=[];
    rd_final10=[];
    
    for i=1:1:size(ccxe11,1)
        if  ((abs(rd_final11(i,1))) <= (radius_in))
            ccxe10=[ ccxe10; ccxe11(i,:)];
            ccye10=[ ccye10; ccye11(i,:)];
            rd_final10=[rd_final10; rd_final11(i,:)];
        end
    end
    % revised corner
    FINAL_CORNER2=[];
    p=1;
    for i= 1: size(ccxe10,1)
        if i == size(ccxe10,1)
            for j= 1: size(ccxe10,1)
                if FINAL_CORNER(p,3)== ccxe10(i,2) && FINAL_CORNER(p+1,3) == ccxe10(i,3)
                    FINAL_CORNER2=  [FINAL_CORNER2; FINAL_CORNER(p,:); FINAL_CORNER(p+1,:)];
                    
                else
                    p=p+2;
                    if p == size(FINAL_CORNER,1)-1
                        FINAL_CORNER2=  [FINAL_CORNER2; FINAL_CORNER(p,:); FINAL_CORNER(p+1,:)];
                    end
                end
                if p ~= size(FINAL_CORNER,1)-1
                    continue
                else
                    break
                end
            end
            break
        end
        if ccxe10(i,2)~= ccxe10(i+1,2)
            for j= 1: size(ccxe10,1)
                if FINAL_CORNER(p,3)== ccxe10(i,2) && FINAL_CORNER(p+1,3) == ccxe10(i,3)
                    FINAL_CORNER2=  [FINAL_CORNER2; FINAL_CORNER(p,:); FINAL_CORNER(p+1,:)];
                    p=p+2;
                    break
                else
                    p=p+2;
                    continue
                end
            end
        end
    end
    CORNER_PLOT=[];
    for i=1: size(FINAL_CORNER2,1)
        kk=0;
        if mod(i,2)~= 0
            if FINAL_CORNER2(i,4)> FINAL_CORNER2(i+1,4)
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER2(i,4): size(xy_r), :); xy_r(1: FINAL_CORNER2(i+1,4), :)];
                
            else
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER2(i,4): FINAL_CORNER2(i+1,4) , :)];
                
            end
        end
    end
%     figure, imshow(bw)
%         hold on
%         plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%         hold on
%         plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'g*','LineWidth',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','k',...
%             'MarkerSize',2)
%         hold on
%         plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','k',...
%             'MarkerSize',2)
%         hold on
%     for i=1 : size(ccxe10,1)
%         syms x y
%         f(x,y)=(x-ccxe10(i,1))^2 + (y-ccye10(i,1))^2 - rd_final10(i,1)^2;
%         h=ezplot(f, [-1000,4000,0,4000]);
%         set(h, 'Color', 'm');
%         hold on
%     end
%     
    save('FINAL_CORNER_AN', 'FINAL_CORNER2')
    save('CORNER_PLOT', 'CORNER_PLOT')
    %% finding maximum curvature point
    xy_r1=[];
    for i=1: size(xy_r,1)
        xy_r1= [xy_r1; xy_r(i,:) 0];
    end
    dist_p=[];
    xh=[];
    yh=[];
    indexh=[];
    for i=1: size(FINAL_CORNER2,1)
        if mod(i,2)==0
            continue
        end
        if FINAL_CORNER2(i,4)< FINAL_CORNER2(i+1,4)
            for j= FINAL_CORNER2(i,4) : FINAL_CORNER2(i+1,4) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsp,indh] = max(dist_p(:,3));
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            dist_p=[];
        else
            for j= FINAL_CORNER2(i,4) : size(xy_r1,1) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            for j= 1 : FINAL_CORNER2(i+1,4) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsp,indh] = max(dist_p(:,3));
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            dist_p=[];
        end
    end
    save('xh', 'xh')
    save('yh', 'yh')
    % hold on
    % plot(xh(:, 1),yh(:, 1), 'k*','MarkerSize',7)
    % title('corners with maximum curvature points')
    % hold on
    % for i =1: size(ccxe,1)
    %     syms x y
    %     f(x,y)=(x-ccxe(i,1))^2 + (y-ccye(i,1))^2 - rd_final(i,1)^2;
    %     h=ezplot(f, [-1000,4000,0,4000]);
    %     set(h, 'Color', 'g');
    %     hold on
    % end
    
    
    % to find circles tangential to stationary point
    
    ccxe13=[];
    ccye13=[];
    rd_final13=[];
    p=1;
    for i=1:1:size(ccxe10,1)
        
        if i == size(ccxe10,1)
            if abs(rd_final10(i,1) - sqrt((ccxe10(i,1)- xh(p,1))^2 + (ccye10(i,1)- yh(p,1))^2)) < 2
                ccxe13=[ ccxe13; ccxe10(i,:)];
                ccye13=[ ccye13; ccye10(i,:)];
                rd_final13=[rd_final13; rd_final10(i,:)];
            end
            break
        end
        if (ccxe10(i,3)== ccxe10(i+1,3))
            if abs(abs(rd_final10(i,1)) - sqrt((ccxe10(i,1)- xh(p,1))^2 + (ccye10(i,1)- yh(p,1))^2)) < 2
                ccxe13=[ ccxe13; ccxe10(i,:)];
                ccye13=[ ccye13; ccye10(i,:)];
                rd_final13=[rd_final13; rd_final10(i,:)];
            end
        else
            p=p+1;
        end
    end
    
   %revised corner
    FINAL_CORNER3=[];
    p=1;
    for i= 1: size(ccxe13,1)
        if i == size(ccxe13,1)
            for j= 1: size(ccxe13,1)
                if FINAL_CORNER2(p,3)== ccxe13(i,2) && FINAL_CORNER2(p+1,3) == ccxe13(i,3)
                    FINAL_CORNER3=  [FINAL_CORNER3; FINAL_CORNER2(p,:); FINAL_CORNER2(p+1,:)];
                    
                else
                    p=p+2;
                    if p == size(FINAL_CORNER2,1)-1
                        FINAL_CORNER3=  [FINAL_CORNER3; FINAL_CORNER2(p,:); FINAL_CORNER2(p+1,:)];
                    end
                end
                if p ~= size(FINAL_CORNER2,1)-1
                    continue
                else
                    break
                end
            end
            break
        end
        if ccxe13(i,2)~= ccxe13(i+1,2)
            for j= 1: size(ccxe13,1)
                if FINAL_CORNER2(p,3)== ccxe13(i,2) && FINAL_CORNER2(p+1,3) == ccxe13(i,3)
                    FINAL_CORNER3=  [FINAL_CORNER3; FINAL_CORNER2(p,:); FINAL_CORNER2(p+1,:)];
                    p=p+2;
                    break
                else
                    p=p+2;
                    continue
                end
            end
        end
    end
    CORNER_PLOT=[];
    for i=1: size(FINAL_CORNER3,1)
        kk=0;
        if mod(i,2)~= 0
            if FINAL_CORNER3(i,4)> FINAL_CORNER3(i+1,4)
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER3(i,4): size(xy_r), :); xy_r(1: FINAL_CORNER3(i+1,4), :)];
                
            else
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER3(i,4): FINAL_CORNER3(i+1,4) , :)];
                
            end
        end
    end
    FINAL_CORNER2= FINAL_CORNER3;
    % to update stationary point
    dist_p=[];
    xh=[];
    yh=[];
    indexh=[];
    for i=1: size(FINAL_CORNER2,1)
        if mod(i,2)==0
            continue
        end
        if FINAL_CORNER2(i,4)< FINAL_CORNER2(i+1,4)
            for j= FINAL_CORNER2(i,4) : FINAL_CORNER2(i+1,4) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsp,indh] = max(dist_p(:,3));
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            dist_p=[];
        else
            for j= FINAL_CORNER2(i,4) : size(xy_r1,1) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            for j= 1 : FINAL_CORNER2(i+1,4) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsp,indh] = max(dist_p(:,3));
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            dist_p=[];
        end
    end
    save('xh', 'xh')
    save('yh', 'yh')
    
    
    
    
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%     hold on
%     plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'g*','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
%     plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
%     for i=1:size(ccxe13,1)
%         syms x y
%         f(x,y)=(x-ccxe13(i,1))^2 + (y-ccye13(i,1))^2 - rd_final13(i,1)^2;
%         h=ezplot(f, [-1000,4000,0,4000]);
%         set(h, 'Color', 'g');
%         hold on
%     end
    
    % to filter out circles wich are going out from ends
    ccxei=[];
    ccyei=[];
    rd_finali=[];
    ccxei1=[];
    ccyei1=[];
    rd_finali1=[];
    p=1;
    p1=1;
    d1=[];
    d2=[];
    for i=1:1:size(ccxe13,1)
        if i == size(ccxe13,1) || (ccxe13(i,3)~= ccxe13(i+1,3))
            if size(ccxei1,1)~= 0
                [l,ind]= min(d1(:,1));
                [l1,ind1]= min(d2(:,1));
                ccxei=[ ccxei; ccxei1(ind,:); ccxei1(ind1,:)];
                ccyei=[ ccyei; ccyei1(ind,:); ccyei1(ind1,:)];
                rd_finali=[rd_finali; rd_finali1(ind,:); rd_finali1(ind1,:)];
            end
            if i ~= size(ccxe13,1)
                p=p+1;
                p1=p1+2;
            end
            d1=[];
            d2=[];
            ccxei1=[];
            ccyei1=[];
            rd_finali1=[];
        end
        
        if i == size(ccxe13,1)
            if abs(rd_final13(i,1)) <=  sqrt((ccxe13(i,1)- xh(p,1))^2 + (ccye13(i,1)- yh(p,1))^2) && abs(rd_final13(i,1)) <=  sqrt((ccxe13(i,1)- FINAL_CORNER2(p1,1))^2 + (ccye13(i,1)- FINAL_CORNER2(p1,2))^2) && abs(rd_final13(i,1)) <=  sqrt((ccxe13(i,1)- FINAL_CORNER2(p1+1,1))^2 + (ccye13(i,1)- FINAL_CORNER2(p1+1,2))^2)
                ccxei=[ ccxei; ccxe13(i,:)];
                ccyei=[ ccyei; ccye13(i,:)];
                rd_finali=[rd_finali; rd_final13(i,:)];
                
            end
            break
        end
        
        if (ccxe13(i,3)== ccxe13(i+1,3))
            
            if abs(rd_final13(i,1)) <=  sqrt((ccxe13(i,1)- xh(p,1))^2 + (ccye13(i,1)- yh(p,1))^2) && abs(rd_final13(i,1)) <=  sqrt((ccxe13(i,1)- FINAL_CORNER2(p1,1))^2 + (ccye13(i,1)- FINAL_CORNER2(p1,2))^2) && abs(rd_final13(i,1)) <=  sqrt((ccxe13(i,1)- FINAL_CORNER2(p1+1,1))^2 + (ccye13(i,1)- FINAL_CORNER2(p1+1,2))^2)
                
                ccxei=[ ccxei; ccxe13(i,:)];
                ccyei=[ ccyei; ccye13(i,:)];
                rd_finali=[rd_finali; rd_final13(i,:)];
            else
                d1=[d1; abs(rd_final13(i,1))- sqrt((ccxe13(i,1)- FINAL_CORNER2(p1,1))^2 + (ccye13(i,1)- FINAL_CORNER2(p1,2))^2)];
                d2=[d2; abs(rd_final13(i,1)) -  sqrt((ccxe13(i,1)- FINAL_CORNER2(p1+1,1))^2 + (ccye13(i,1)- FINAL_CORNER2(p1+1,2))^2)];
                ccxei1= [ ccxei1; ccxe13(i,:)];
                ccyei1=[ ccyei1; ccye13(i,:)];
                rd_finali1=[rd_finali1; rd_final13(i,:)];
            end
            
        end
    end
    
    %revised corner
    FINAL_CORNER3=[];
    p=1;
    for i= 1: size(ccxei,1)
        if i == size(ccxei,1)
            for j= 1: size(ccxei,1)
                if FINAL_CORNER2(p,3)== ccxei(i,2) && FINAL_CORNER2(p+1,3) == ccxei(i,3)
                    FINAL_CORNER3=  [FINAL_CORNER3; FINAL_CORNER2(p,:); FINAL_CORNER2(p+1,:)];
                    
                else
                    p=p+2;
                    if p == size(FINAL_CORNER2,1)-1
                        FINAL_CORNER3=  [FINAL_CORNER3; FINAL_CORNER2(p,:); FINAL_CORNER2(p+1,:)];
                    end
                end
                if p ~= size(FINAL_CORNER2,1)-1
                    continue
                else
                    break
                end
            end
            break
        end
        if ccxei(i,2)~= ccxei(i+1,2)
            for j= 1: size(ccxei,1)
                if FINAL_CORNER2(p,3)== ccxei(i,2) && FINAL_CORNER2(p+1,3) == ccxei(i,3)
                    FINAL_CORNER3=  [FINAL_CORNER3; FINAL_CORNER2(p,:); FINAL_CORNER2(p+1,:)];
                    p=p+2;
                    break
                else
                    p=p+2;
                    continue
                end
            end
        end
    end
    CORNER_PLOT=[];
    for i=1: size(FINAL_CORNER3,1)
        kk=0;
        if mod(i,2)~= 0
            if FINAL_CORNER3(i,4)> FINAL_CORNER3(i+1,4)
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER3(i,4): size(xy_r), :); xy_r(1: FINAL_CORNER3(i+1,4), :)];
                
            else
                CORNER_PLOT=[CORNER_PLOT; xy_r(FINAL_CORNER3(i,4): FINAL_CORNER3(i+1,4) , :)];
                
            end
        end
    end
    FINAL_CORNER2= FINAL_CORNER3;
    % to update stationary point
    dist_p=[];
    xh=[];
    yh=[];
    indexh=[];
    for i=1: size(FINAL_CORNER2,1)
        if mod(i,2)==0
            continue
        end
        if FINAL_CORNER2(i,4)< FINAL_CORNER2(i+1,4)
            for j= FINAL_CORNER2(i,4) : FINAL_CORNER2(i+1,4) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsp,indh] = max(dist_p(:,3));
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            dist_p=[];
        else
            for j= FINAL_CORNER2(i,4) : size(xy_r1,1) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            for j= 1 : FINAL_CORNER2(i+1,4) % to run for each segment
                pt = [xy_r1(j,1),xy_r1(j,2),0];
                a = xy_r1(FINAL_CORNER2(i,4),:)-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                b = pt-xy_r1(FINAL_CORNER2(i+1,4),:);%Vector
                dist_p = [dist_p; xy_r1(j,2) xy_r1(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsp,indh] = max(dist_p(:,3));
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            dist_p=[];
        end
    end
    save('xh', 'xh')
    save('yh', 'yh')
    
    
    
    %
%     figure, imshow(bw)
%     hold on
%     plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%     hold on
%     plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'g*','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
%     plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',4,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
%     for i=1:size(ccxei,1)
%         syms x y
%         f(x,y)=(x-ccxei(i,1))^2 + (y-ccyei(i,1))^2 - rd_finali(i,1)^2;
%         h=ezplot(f, [-1000,4000,0,4000]);
%         set(h, 'Color', 'b');
%         hold on
%     end
%     
%     
    %% circles tangential to maximum number of dominant points minimum error
    % to find distance of each circle centre in the corner from its corner
    % dominant point and finding the circle with maximum number of tangency
    % points and minimum error among the ones which have the same number of
    % tangency points
    count=[];
    p=1;
    err=[];
    dp=[];
    for i=1:size(ccxei,1)
        c=0;
        dp=[];
        if i == size(ccxei,1) % for last circle
            if ccxei(i,2)< ccxei(i,3)
                for j = ccxei(i,2) : ccxei(i,3)
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                
                if size(dp,1)==0
                    err=[err; 0];
                else
                    err=[err; mean(dp,1)];
                end
                count=[count; c ccxei(i,1) ccyei(i,1) rd_finali(i,1) ccxei(i,2) ccxei(i,3)];
            else
                for j = ccxei(i,2) : size(COOR_DPF1,1)/2
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                for j = 1 : ccxei(i,3)
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                
                count=[count; c ccxei(i,1) ccyei(i,1) rd_finali(i,1) ccxei(i,2) ccxei(i,3)];
                if size(dp,1)==0
                    err=[err; 0];
                else
                    err=[err; mean(dp,1)];
                end
            end
            break
        end
        if ccxei(i,3)== ccxei(i+1,3) % for same corner
            if ccxei(i,2)< ccxei(i,3)
                for j = ccxei(i,2) : ccxei(i,3)
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                
                count=[count; c ccxei(i,1) ccyei(i,1) rd_finali(i,1) ccxei(i,2) ccxei(i,3)];
                if size(dp,1)==0
                    err=[err; 0];
                else
                    err=[err; mean(dp,1)];
                end
            else % for re entrant corner
                for j = ccxei(i,2) : size(COOR_DPF1,1)/2
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                for j = 1 : ccxei(i,3)
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                
                count=[count; c ccxei(i,1) ccyei(i,1) rd_finali(i,1) ccxei(i,2) ccxei(i,3)];
                if size(dp,1)==0
                    err=[err; 0];
                else
                    err=[err; mean(dp,1)];
                end
            end
        else % for change of corner
            if ccxei(i,2)< ccxei(i,3)
                for j = ccxei(i,2) : ccxei(i,3)
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                count=[count; c ccxei(i,1) ccyei(i,1) rd_finali(i,1) ccxei(i,2) ccxei(i,3)];
                count=[count; 0 0 0 0 0 0];
                if size(dp,1)==0
                    err=[err; 0];
                else
                    err=[err; mean(dp,1)];
                end
                err=[err;1000];
                
            else
                for j = ccxei(i,2) : size(COOR_DPF1,1)/2
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                for j = 1 : ccxei(i,3)
                    dpo= ((sqrt((COOR_DPF1(j,1)- ccxei(i,1))^2 + (COOR_DPF1(j,2)- ccyei(i,1))^2 ))- abs(rd_finali(i,1)))^2;
                    if dpo < p
                        c=c+1;
                        dp=[dp; dpo];
                    end
                end
                count=[count; c ccxei(i,1) ccyei(i,1) rd_finali(i,1) ccxei(i,2) ccxei(i,3)];
                count=[count; 0 0 0 0 0 0];
                if size(dp,1)==0
                    err=[err; 0];
                else
                    err=[err; mean(dp,1)];
                end
                err=[err;1000];
                
            end
        end
    end
    err=[err;1000];
    
    % to find circles in corners with maximum no. of tangency points
    p=1;
    for i= 1: size(count,1)
        if count(i,2)==0
            p=[p;i];
        end
    end
    p=[p; size(count,1)];
    ccxe2=[];
    ccye2=[];
    rd_final2=[];
    ccxe3=[0 0 0];
    ccye3=0;
    rd_final3=0;
    err2=[];
    err3=10;
    for i= 1: size(p,1)-1
        [idx,k]= max(count(p(i,1):p(i+1,1),1));
        [r c] = find(count(p(i,1):p(i+1,1)) == idx);
        ccxe2=[ccxe2; count((p(i,1)+c-1),2)];
        ccye2=[ccye2; count((p(i,1)+c-1),3)];
        rd_final2=[rd_final2; count((p(i,1)+c-1),4)];
        err2=[err2;err((p(i,1)+c-1),1)];
        ccxe3=[ccxe3; count((p(i,1)+c-1),2) count((p(i,1)+c-1),5) count((p(i,1)+c-1),6) ;0 0 0];
        ccye3=[ccye3; count((p(i,1)+c-1),3);0];
        rd_final3=[rd_final3; count((p(i,1)+c-1),4);0];
        err3=[err3;err((p(i,1)+c-1),1);10];
    end
    % to plot circles with maximum tangency
    %     figure, imshow(bw)
    %     hold on
    %     plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
    %     hold on
    %     plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'--gs','LineWidth',2,...
    %         'MarkerEdgeColor','g',...
    %         'MarkerFaceColor','k',...
    %         'MarkerSize',2)
    %     hold on
    %     plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'--gs','LineWidth',2,...
    %         'MarkerEdgeColor','g',...
    %         'MarkerFaceColor','k',...
    %         'MarkerSize',2)
    %     hold on
    %     for i =1: size(ccxe2,1)
    %         syms x y
    %         f(x,y)=(x-ccxe2(i,1))^2 + (y-ccye2(i,1))^2 - rd_final2(i,1)^2;
    %         h=ezplot(f, [-1000,4000,0,4000]);
    %         set(h, 'Color', 'b');
    %         hold on
    %     end
    
    % to get circles with minimum error and maximum no of tangent points
    po=[];
    for i= 1: size(err3,1)
        if err3(i,1)==10
            po=[po;i];
        end
    end
    rd_final22=[];
    ccye22=[];
    ccxe22=[];
    for i=1: size(po,1)-1
        [idx,k]= min(err3(po(i,1):po(i+1,1),1));
        ccxe22=[ccxe22; ccxe3((po(i,1)+k-1),:)];
        ccye22=[ccye22; ccye3((po(i,1)+k-1),1)];
        rd_final22=[rd_final22; rd_final3((po(i,1)+k-1),1)];
    end
    % to plot circles with minimum error
%         figure, imshow(bw)
%         hold on
%         plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
%         hold on
%         plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'--gs','LineWidth',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','k',...
%             'MarkerSize',2)
%         hold on
%         plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'--gs','LineWidth',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','k',...
%             'MarkerSize',2)
%         hold on
%         for i =1: size(ccxe22,1)
%             syms x y
%             f(x,y)=(x-ccxe22(i,1))^2 + (y-ccye22(i,1))^2 - rd_final22(i,1)^2;
%             h=ezplot(f, [-1000,4000,0,4000]);
%             set(h, 'Color', 'b');
%             hold on
%         end
    
    
    %% circles inside figure from outside the extreme points of a corner
    centrex=[];
    centrey=[];
    radiusxy=[];
    p=1;
    g1=1;
    NZ=[];
    NZ1=[];
    TRIM_CORNER=[];
    for i= 1:1:size(ccxe22,1)
        x1 = ccxe22(i,1);
        y1 = ccye22(i,1);
        c= [x1+rd_final22(i,1)*cos(0), y1+rd_final22(i,1)*sin(0)];
        DirVector1=[ccxe22(i,1),ccye22(i,1)]-c;
        DirVector2=[ccxe22(i,1),ccye22(i,1)]-[COOR_DPF1(ccxe22(i,2),1),COOR_DPF1(ccxe22(i,2),2)];
        DirVector3=[ccxe22(i,1),ccye22(i,1)]-[COOR_DPF1(ccxe22(i,3),1),COOR_DPF1(ccxe22(i,3),2)];
        Angle1=acos( dot(DirVector1,DirVector2)/norm(DirVector1)/norm(DirVector2) );
        Angle2=acos( dot(DirVector1,DirVector3)/norm(DirVector1)/norm(DirVector3) );
        tmp=[];
        for angle = 0:0.017:2*pi
            x1 = ccxe22(i,1);
            y1 = ccye22(i,1);
            tmp=[tmp; x1+rd_final22(i,1)*cos(angle)  y1+rd_final22(i,1)*sin(angle)];
        end
        angle= 0:0.017:2*pi;
        Y= [COOR_DPF1(ccxe22(i,2),1),COOR_DPF1(ccxe22(i,2),2)];
        IDX = knnsearch(tmp,Y) ;
        if IDX > 185
            a1= 2*pi-Angle1;
        else
            a1= Angle1;
        end
        Y= [COOR_DPF1(ccxe22(i,3),1),COOR_DPF1(ccxe22(i,3),2)];
        IDX = knnsearch(tmp,Y) ;
        if IDX > 185
            a2= 2*pi-Angle2;
        else
            a2= Angle2;
        end
        Y= [xh(p,1),yh(p,1)];
        IDX = knnsearch(tmp,Y) ;
        if a1<a2
            p1= a1 :0.017: a2;
            p1=p1.';
            if abs(p1(knnsearch(p1,angle(1,IDX)),1)- angle(1,IDX))>0.01
                tmp=[];
                tmp_rev=[];
                for a = a1 :0.017: a2
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp=[tmp; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                for a = a2:0.017:2*pi
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp_rev=[tmp_rev; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                for a = 0:0.017:a1
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp_rev=[tmp_rev; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
            else
                tmp=[];
                tmp_rev=[];
                for a = a2:0.017:2*pi
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp=[tmp; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                
                for a = 0:0.017:a1
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp=[tmp; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                
                for a = a1 :0.017: a2
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp_rev=[tmp_rev; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
            end
        else
            p1=a1 :0.017: 2*pi;
            p1=[p1 0:0.017:a2];
            p1=p1.';
            if abs(p1(knnsearch(p1,angle(1,IDX)),1)- angle(1,IDX))< 0.01
                tmp=[];
                tmp_rev=[];
                for a = a2:0.017: a1
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp=[tmp; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                for a = a1:0.017:2*pi
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp_rev=[tmp_rev; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                
                for a = 0:0.017:a2
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp_rev=[tmp_rev; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
            else
                tmp=[];
                tmp_rev=[];
                for a = a1:0.017:2*pi
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp=[tmp; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                
                for a = 0:0.017:a2
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp=[tmp; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
                for a = a2:0.017: a1
                    x1 = ccxe22(i,1);
                    y1 = ccye22(i,1);
                    tmp_rev=[tmp_rev; x1+rd_final22(i,1)*cos(a)  y1+rd_final22(i,1)*sin(a)];
                end
            end
        end
        % trimmed corner
        sx=[];
        for k=1: size(tmp_rev,1)
          Y= [tmp_rev(k,2) tmp_rev(k,1) ];
          IDX = knnsearch(xy_r,Y) ;
          dp= sqrt((xy_r(IDX,2)- tmp_rev(k,1))^2 + (xy_r(IDX,1)- tmp_rev(k,2))^2);
          if dp <= 5
           TRIM_CORNER=[ TRIM_CORNER ; xy_r(IDX,2) xy_r(IDX,1) g1 ]; 
           sx= [sx; xy_r(IDX,2) xy_r(IDX,1) g1 ];
           break
          end
        end
        for k =size(tmp_rev,1):-1: 1
          Y= [tmp_rev(k,2) tmp_rev(k,1) ];
          IDX = knnsearch(xy_r,Y) ;
          dp= sqrt((xy_r(IDX,2)- tmp_rev(k,1))^2 + (xy_r(IDX,1)- tmp_rev(k,2))^2);
          if dp <= 5
           TRIM_CORNER=[ TRIM_CORNER ; xy_r(IDX,2) xy_r(IDX,1) g1 ]; 
           sx= [sx; xy_r(IDX,2) xy_r(IDX,1) g1 ];
           break
          end
        end
        if size(sx,1) == 0
         TRIM_CORNER=[ TRIM_CORNER ; FINAL_CORNER2(g1,1:2) g1 ; FINAL_CORNER2(g1+1,1:2) g1 ];    
        end
                        
        dist_p=[];
        for l= 1 : size(tmp_rev,1)
            if inpolygon(tmp_rev(l,1),tmp_rev(l,2),xy_r(:,2),xy_r(:,1)) == 0
                Y= [tmp_rev(l,2) tmp_rev(l,1) ];
                IDX = knnsearch(xy_r,Y) ;
                dist_p=[dist_p; tmp_rev(l,1) tmp_rev(l,2) sqrt((xy_r(IDX,1)-tmp_rev(l,2))^2 + (xy_r(IDX,2)-tmp_rev(l,1))^2 )];
            end
        end
        if size(dist_p,1) == 0
            dsp = 0;
        else
            [dsp,indh] = max(dist_p(:,3));
            x_rev= dist_p(indh,1);
            y_rev=  dist_p(indh,2);
        end
        dist_p=[];
        ind=[];
        for j= 1: size(tmp,1)
            ind = [ind; inpolygon(tmp(j,1),tmp(j,2),xy_r(:,2),xy_r(:,1))];
        end
        nz=0;
        for k = 1: size(ind,1)
            if ind(k,1)== 0
                nz=nz+1;
            end
        end
        NZ=[NZ;nz];
        NZ1=[NZ1; nz i size(ind,1)];
        if min(ind(:,1) == 1) && dsp <= 5
            centrex=[ centrex; ccxe22(i,:)];
            centrey=[ centrey; ccye22(i,:)];
            radiusxy=[radiusxy; rd_final22(i,:)];
        end
        if min(NZ) ~=0 || (min(NZ) == 0 && dsp>5)
            [k,h]= min(NZ);
            if ((k <= 0.2*NZ1(h,3) && rd_final22(i,1)<=60) || (k <= 0.1*NZ1(h,3))) && dsp <= 5
                centrex=[ centrex; ccxe22(NZ1(h,2),:)];
                centrey=[centrey; ccye22(NZ1(h,2),:)];
                radiusxy=[radiusxy; rd_final22(NZ1(h,2),:)];
            else % if no circle finally satisfied in a corner
                ccxe_st=[];
                ccye_st=[];
                rd_st=[];
                % finding all the circles tangent to the stationary point in that corner
                for j= 1: size(ccxe13,1)
                    if ccxe13(j,2) == ccxe22(i,2) && ccxe13(j,3) == ccxe22(i,3)
                        ccxe_st = [ccxe_st; ccxe13(j,:)];
                        ccye_st = [ccye_st; ccye13(j,:)];
                        rd_st = [rd_st; rd_final13(j,:)];
                    end
                end
%                                                 figure, imshow(bw)
%                                                 hold on
%                                                 plot(xy_r(:, 2),xy_r(:, 1), 'm', 'LineWidth', 2)
%                                                 hold on
%                                                 plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'g*','LineWidth',2,...
%                                                     'MarkerEdgeColor','g',...
%                                                     'MarkerFaceColor','k',...
%                                                     'MarkerSize',2)
%                                                 hold on
%                                                 plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',4,...
%                                                     'MarkerEdgeColor','r',...
%                                                     'MarkerFaceColor','k',...
%                                                     'MarkerSize',2)
%                                                 hold on
%                                                 for j=1:size(ccxe_st,1)
%                                                     syms x y
%                                                     f(x,y)=(x-ccxe_st(j,1))^2 + (y-ccye_st(j,1))^2 - rd_st(j,1)^2;
%                                                     h=ezplot(f, [-1000,4000,0,4000]);
%                                                     set(h, 'Color', 'g');
%                                                     hold on
%                                                 end
%                 %
                % finding the circles inside the particle from extreme points in that corner
                p2=p;
                NZ1i=[];
                NZi=[];
                ccxe_ex=[];
                ccye_ex=[];
                rd_ex=[];
                for j= 1:1:size(ccxe_st,1)
                    x1 = ccxe_st(j,1);
                    y1 = ccye_st(j,1);
                    c= [x1+rd_st(j,1)*cos(0), y1+rd_st(j,1)*sin(0)];
                    DirVector1=[ccxe_st(j,1),ccye_st(j,1)]-c;
                    DirVector2=[ccxe_st(j,1),ccye_st(j,1)]-[COOR_DPF1(ccxe_st(j,2),1),COOR_DPF1(ccxe_st(j,2),2)];
                    DirVector3=[ccxe_st(j,1),ccye_st(j,1)]-[COOR_DPF1(ccxe_st(j,3),1),COOR_DPF1(ccxe_st(j,3),2)];
                    Angle1=acos( dot(DirVector1,DirVector2)/norm(DirVector1)/norm(DirVector2) );
                    Angle2=acos( dot(DirVector1,DirVector3)/norm(DirVector1)/norm(DirVector3) );
                    tmp=[];
                    for angle = 0:0.017:2*pi
                        x1 = ccxe_st(j,1);
                        y1 = ccye_st(j,1);
                        tmp=[tmp; x1+rd_st(j,1)*cos(angle)  y1+rd_st(j,1)*sin(angle)];
                    end
                    angle= 0:0.017:2*pi;
                    Y= [COOR_DPF1(ccxe_st(j,2),1),COOR_DPF1(ccxe_st(j,2),2)];
                    IDX = knnsearch(tmp,Y) ;
                    if IDX > 185
                        a1= 2*pi-Angle1;
                    else
                        a1= Angle1;
                    end
                    Y= [COOR_DPF1(ccxe_st(j,3),1),COOR_DPF1(ccxe_st(j,3),2)];
                    IDX = knnsearch(tmp,Y) ;
                    if IDX > 185
                        a2= 2*pi-Angle2;
                    else
                        a2= Angle2;
                    end
                    Y= [xh(p2,1),yh(p2,1)];
                    IDX = knnsearch(tmp,Y) ;
                    if a1<a2
                        p1= a1 :0.017: a2;
                        p1=p1.';
                        if abs(p1(knnsearch(p1,angle(1,IDX)),1)- angle(1,IDX))>0.01
                            tmp=[];
                            tmp_rev=[];
                            for a = a1 :0.017: a2
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp=[tmp; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            for a = a2:0.017:2*pi
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp_rev=[tmp_rev; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            
                            for a = 0:0.017:a1
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp_rev=[tmp_rev; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                        else
                            tmp=[];
                            tmp_rev=[];
                            for a = a2:0.017:2*pi
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp=[tmp; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            
                            for a = 0:0.017:a1
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp=[tmp; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            for a = a1 :0.017: a2
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp_rev=[tmp_rev; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                        end
                        
                    else
                        p1=a1 :0.017: 2*pi;
                        p1=[p1 0:0.017:a2];
                        p1=p1.';
                        if abs(p1(knnsearch(p1,angle(1,IDX)),1)- angle(1,IDX))< 0.01
                            tmp=[];
                            tmp_rev=[];
                            for a = a2:0.017: a1
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp=[tmp; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            for a = a1:0.017:2*pi
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp_rev=[tmp_rev; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            
                            for a = 0:0.017:a2
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp_rev=[tmp_rev; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                        else
                            tmp=[];
                            tmp_rev=[];
                            for a = a1:0.017:2*pi
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp=[tmp; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            
                            for a = 0:0.017:a2
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp=[tmp; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                            for a = a2:0.017: a1
                                x1 = ccxe_st(j,1);
                                y1 = ccye_st(j,1);
                                tmp_rev=[tmp_rev; x1+rd_st(j,1)*cos(a)  y1+rd_st(j,1)*sin(a)];
                            end
                        end
                    end
                    ind=[];
                    for k= 1: size(tmp,1)
                        ind = [ind; inpolygon(tmp(k,1),tmp(k,2),xy_r(:,2),xy_r(:,1))];
                    end
                    nz=0;
                    for k = 1: size(ind,1)
                        if ind(k,1)== 0
                            nz=nz+1;
                        end
                    end
                    NZi=[NZi; nz];
                    NZ1i=[NZ1i; nz j size(ind,1)];
                    dist_p=[];
                    for l= 1 : size(tmp_rev,1)
                        if inpolygon(tmp_rev(l,1),tmp_rev(l,2),xy_r(:,2),xy_r(:,1)) == 0
                            Y= [tmp_rev(l,2) tmp_rev(l,1) ];
                            IDX = knnsearch(xy_r,Y) ;
                            dist_p=[dist_p; tmp_rev(l,1) tmp_rev(l,2) sqrt((xy_r(IDX,1)-tmp_rev(l,2))^2 + (xy_r(IDX,2)-tmp_rev(l,1))^2 )];
                        end
                    end
                    if size(dist_p,1) == 0
                        dsp = 0;
                    else
                        [dsp,indh] = max(dist_p(:,3));
                        x_rev= dist_p(indh,1);
                        y_rev=  dist_p(indh,2);
                    end
                    dist_p=[];
                    if min(ind(:,1) == 1) && dsp <= 15
                        ccxe_ex=[ ccxe_ex; ccxe_st(j,:)];
                        ccye_ex=[ ccye_ex; ccye_st(j,:)];
                        rd_ex=[rd_ex; rd_st(j,:)];
                    end
                    if  j== size(ccxe_st,1)
                        if min(NZi) ~=0 || (min(NZi) == 0 && dsp > 15)
                            [k,h]= min(NZi);
                            if ((k <= 0.2*NZ1i(h,3) && rd_st(j,1)<=60) || (k <= 0.1*NZ1i(h,3))) && dsp <= 15
                                ccxe_ex=[ ccxe_ex; ccxe_st(NZ1i(h,2),:)];
                                ccye_ex=[ccye_ex; ccye_st(NZ1i(h,2),:)];
                                rd_ex=[rd_ex; rd_st(NZ1i(h,2),:)];
                            end
                        end
                    end
                    NZi=[];
                    NZ1i=[];
                end
                if size(ccxe_ex,1) == 0
                    p=p+1;
                    g1=g1+2;
                    NZ=[];
                    NZ1=[];
                    continue
                end
%                                 figure, imshow(bw)
%                                 hold on
%                                 plot(xy_r(:, 2),xy_r(:, 1), 'm', 'LineWidth', 2)
%                                 hold on
%                                 plot(FINAL_CORNER2(:,1),FINAL_CORNER2(:,2),'g*','LineWidth',2,...
%                                     'MarkerEdgeColor','g',...
%                                     'MarkerFaceColor','k',...
%                                     'MarkerSize',2)
%                                 hold on
%                                 plot(CORNER_PLOT(:,2),CORNER_PLOT(:,1),'g*','LineWidth',4,...
%                                     'MarkerEdgeColor','r',...
%                                     'MarkerFaceColor','k',...
%                                     'MarkerSize',2)
%                                 hold on
%                                 for j=1:size(ccxe_ex,1)
%                                     syms x y
%                                     f(x,y)=(x-ccxe_ex(j,1))^2 + (y-ccye_ex(j,1))^2 - rd_ex(j,1)^2;
%                                     h=ezplot(f, [-1000,4000,0,4000]);
%                                     set(h, 'Color', 'g');
%                                     hold on
%                                 end
                %
                % finding the maximum tangential and minimum error circle in that corner finally
                
                count=[];
                p3=1;
                err=[];
                dp=[];
                for k=1:size(ccxe_ex,1)
                    c=0;
                    dp=[];
                    if k == size(ccxe_ex,1) % for last circle
                        if ccxe_ex(k,2)< ccxe_ex(k,3)
                            for j = ccxe_ex(k,2) : ccxe_ex(k,3)
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            
                            if size(dp,1)==0
                                err=[err; 0];
                            else
                                err=[err; mean(dp,1)];
                            end
                            count=[count; c ccxe_ex(k,1) ccye_ex(k,1) rd_ex(k,1) ccxe_ex(k,2) ccxe_ex(k,3)];
                        else
                            for j = ccxe_ex(k,2) : size(COOR_DPF1,1)/2
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            for j = 1 : ccxe_ex(k,3)
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            
                            count=[count; c ccxe_ex(k,1) ccye_ex(k,1) rd_ex(k,1) ccxe_ex(k,2) ccxe_ex(k,3)];
                            if size(dp,1)==0
                                err=[err; 0];
                            else
                                err=[err; mean(dp,1)];
                            end
                        end
                        break
                    end
                    if ccxe_ex(k,3)== ccxe_ex(k+1,3) % for same corner
                        if ccxe_ex(k,2)< ccxe_ex(k,3)
                            for j = ccxe_ex(k,2) : ccxe_ex(k,3)
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            
                            count=[count; c ccxe_ex(k,1) ccye_ex(k,1) rd_ex(k,1) ccxe_ex(k,2) ccxe_ex(k,3)];
                            if size(dp,1)==0
                                err=[err; 0];
                            else
                                err=[err; mean(dp,1)];
                            end
                        else % for re entrant corner
                            for j = ccxe_ex(k,2) : size(COOR_DPF1,1)/2
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            for j = 1 : ccxe_ex(k,3)
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            
                            count=[count; c ccxe_ex(k,1) ccye_ex(k,1) rd_ex(k,1) ccxe_ex(k,2) ccxe_ex(k,3)];
                            if size(dp,1)==0
                                err=[err; 0];
                            else
                                err=[err; mean(dp,1)];
                            end
                        end
                    else % for change of corner
                        if ccxe_ex(k,2)< ccxe_ex(k,3)
                            for j = ccxe_ex(k,2) : ccxe_ex(k,3)
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            count=[count; c ccxe_ex(k,1) ccye_ex(k,1) rd_ex(k,1) ccxe_ex(k,2) ccxe_ex(k,3)];
                            count=[count; 0 0 0 0 0 0];
                            if size(dp,1)==0
                                err=[err; 0];
                            else
                                err=[err; mean(dp,1)];
                            end
                            err=[err;1000];
                            
                        else
                            for j = ccxe_ex(k,2) : size(COOR_DPF1,1)/2
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            for j = 1 : ccxe_ex(k,3)
                                dpo= ((sqrt((COOR_DPF1(j,1)- ccxe_ex(k,1))^2 + (COOR_DPF1(j,2)- ccye_ex(k,1))^2 ))- abs(rd_ex(k,1)))^2;
                                if dpo < p3
                                    c=c+1;
                                    dp=[dp; dpo];
                                end
                            end
                            count=[count; c ccxe_ex(k,1) ccye_ex(k,1) rd_ex(k,1) ccxe_ex(k,2) ccxe_ex(k,3)];
                            count=[count; 0 0 0 0 0 0];
                            if size(dp,1)==0
                                err=[err; 0];
                            else
                                err=[err; mean(dp,1)];
                            end
                            err=[err;1000];
                            
                        end
                    end
                end
                err=[err;1000];
                
                % to find circles in corners with maximum no. of tangency points
                p3=1;
                for h= 1: size(count,1)
                    if count(h,2)==0
                        p3=[p3;h];
                    end
                end
                p3=[p3; size(count,1)];
                ccxe_max=[];
                ccye_max=[];
                rd_max=[];
                ccxe_m=[0 0 0];
                ccye_m=0;
                rd_m=0;
                err2=[];
                err3=10;
                for h= 1: size(p3,1)-1
                    [idx,k]= max(count(p3(h,1):p3(h+1,1),1));
                    [r c] = find(count(p3(h,1):p3(h+1,1)) == idx);
                    ccxe_max=[ccxe_max; count((p3(h,1)+c-1),2)];
                    ccye_max=[ccye_max; count((p3(h,1)+c-1),3)];
                    rd_max=[rd_max; count((p3(h,1)+c-1),4)];
                    err2=[err2;err((p3(h,1)+c-1),1)];
                    ccxe_m=[ccxe_m; count((p3(h,1)+c-1),2) count((p3(h,1)+c-1),5) count((p3(h,1)+c-1),6) ;0 0 0];
                    ccye_m=[ccye_m; count((p3(h,1)+c-1),3);0];
                    rd_m=[rd_m; count((p3(h,1)+c-1),4);0];
                    err3=[err3;err((p3(h,1)+c-1),1);10];
                end
                % to get circles with minimum error and maximum no of tangent points
                po=[];
                for h= 1: size(err3,1)
                    if err3(h,1)==10
                        po=[po;h];
                    end
                end
                for h=1: size(po,1)-1
                    [idx,k]= min(err3(po(h,1):po(h+1,1),1));
                    centrex=[centrex; ccxe_m((po(h,1)+k-1),:)];
                    centrey=[centrey; ccye_m((po(h,1)+k-1),1)];
                    radiusxy=[radiusxy; rd_m((po(h,1)+k-1),1)];
                end
            end
        end
        NZ=[];
        p=p+1;
        g1=g1+2;
        NZ1=[];
    end
    toc
    % revised corner and corner plot
    CORNER=[];
    CORNER_PLOT=[];
    p=1;
    for i= 1: size(centrex,1)
        if i == size(centrex,1)
            if FINAL_CORNER2(p,3)== centrex(i,2) && FINAL_CORNER2(p+1,3) == centrex(i,3)
                CORNER=  [CORNER; FINAL_CORNER2(p,:) p; FINAL_CORNER2(p+1,:) p];
                p=p+2;
            else
                p=p+2;
                if FINAL_CORNER2(p,3)== centrex(i,2) && FINAL_CORNER2(p+1,3) == centrex(i,3)
                    CORNER=  [CORNER; FINAL_CORNER2(p,:) p; FINAL_CORNER2(p+1,:) p];
                    p=p+2;
                end
            end
            break
        end
        if centrex(i,2)~= centrex(i+1,2)
            while (i < size(centrex,1))
                if FINAL_CORNER2(p,3)== centrex(i,2) && FINAL_CORNER2(p+1,3) == centrex(i,3)
                    CORNER=  [CORNER; FINAL_CORNER2(p,:) p; FINAL_CORNER2(p+1,:) p];
                    p=p+2;
                    break
                else
                    p=p+2;
                    continue
                end
            end
        end
    end
    for i=1: size(CORNER,1)
        kk=0;
        if mod(i,2)~= 0
            if CORNER(i,4)> CORNER(i+1,4)
                CORNER_PLOT=[CORNER_PLOT; xy_r(CORNER(i,4): size(xy_r), :); xy_r(1:CORNER(i+1,4), :)];
                
            else
                CORNER_PLOT=[CORNER_PLOT; xy_r(CORNER(i,4): CORNER(i+1,4) , :)];
                
            end
        end
    end
    
    % FINAL TRIMMED CORNERS
    % endpoints
    TRIM_CORN=[];
    TRIM_CORNER_PLOT=[];
    q=1;
    for i= 1: size(TRIM_CORNER,1)
        if mod(i,2) == 0
            continue
        end
        if q<= size(CORNER,1) && TRIM_CORNER(i,3) == CORNER(q,end)
            TRIM_CORN=[ TRIM_CORN; TRIM_CORNER(i,:); TRIM_CORNER(i+1,:)];
            q=q+2;
        else
           continue
        end
    end
    % corresponding endpoints in xy_r AND CORNER PLOT
    
    TRIM_CORNER_PLOT=[];
    for i= 1: size(TRIM_CORN,1)+1
        kk=0;
        if mod(i,2)~= 0 && i~=1
            if TRIM_CORN(i-2,4)> TRIM_CORN(i-1,4) && (TRIM_CORN(i-1,3)-TRIM_CORN(i-2,3))< (TRIM_CORN(end,3)*0.5)
                TRIM_CORNER_PLOT=[TRIM_CORNER_PLOT; xy_r(TRIM_CORN(i-2,4): size(xy_r), :); xy_r(1: TRIM_CORN(i-1,4), :)];
            else if TRIM_CORN(i-2,4)< TRIM_CORN(i-1,4) && (TRIM_CORN(i-1,3)-TRIM_CORN(i-2,3))< (TRIM_CORN(end,3)*0.5)
                    TRIM_CORNER_PLOT=[ TRIM_CORNER_PLOT; xy_r(TRIM_CORN(i-2,4): TRIM_CORN(i-1,4) , :)];
                else if TRIM_CORN(i-2,4)> TRIM_CORN(i-1,4) && (TRIM_CORN(i-1,3)-TRIM_CORN(i-2,3))> (TRIM_CORN(end,3)*0.5)
                        TRIM_CORNER_PLOT=[TRIM_CORNER_PLOT; xy_r(TRIM_CORN(i-1,4): TRIM_CORN(i-2,4) , :)];
                    end
                end
            end
        end
        if i <= size(TRIM_CORN,1)
            for j= 1: size(xy_r,1)
                if (TRIM_CORN(i,1)== xy_r(j,2) && TRIM_CORN(i,2)== xy_r(j,1))
                    TRIM_CORN(i,4)= j; % to store relevant index of xy_r matrix
                    kk=kk+1;
                    break
                end
            end
            dmin=[];
            if j== size(xy_r,1) && kk == 0
                for k= 1: size(xy_r,1)
                    dmin = [dmin; sqrt((TRIM_CORN(i,1)-xy_r(k,2))^2 + (TRIM_CORN(i,2)- xy_r(k,1))^2) k];
                end
                [idx,t]= min(dmin(:,1));
                TRIM_CORN(i,4)= t;
            end
        else
            break
        end
        
    end
    
    
    
    %to plot final circles
    h= figure,imshow(bw)

    load('ho')
    load('P')
%     result=glob(ho,filename);
%     he = imread(result{P});% to read the current image
%     figure, imshow(he)
    hold on
    plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
    hold on
%     plot(TRIM_CORNER(:,1),TRIM_CORNER(:,2),'--gs','LineWidth',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','k',...
%         'MarkerSize',2)
%     hold on
    plot(TRIM_CORNER_PLOT(:,2),TRIM_CORNER_PLOT(:,1),'g*','LineWidth',3,...
        'MarkerEdgeColor','g',...
        'MarkerFaceColor','k',...
        'MarkerSize',2)
    hold on
    viscircles([nn,m],[radius_in],'EdgeColor','g', 'LineWidth',2 )
    hold on
    for i =1: size(centrex,1)
        tmp=[];
        for angle = 0:0.017:2*pi
            x1 = centrex(i,1);
            y1 = centrey(i,1);
            tmp=[tmp; x1+radiusxy(i,1)*cos(angle)  y1+radiusxy(i,1)*sin(angle)];
        end
        plot(tmp(:,1), tmp(:,2),'k',  'LineWidth',2.5);
        hold on
    end
    hold on
    num_corn= size(TRIM_CORN,1)/2;
   % to find perimeter ratio
    d=0;
    h=[];
    for i= 1: size(TRIM_CORN,1)
        if mod(i,2)~= 0
            
            for j= TRIM_CORN(i,4): TRIM_CORN(i+1,4)-1
                d= d+ sqrt((xy_r(j,1)-xy_r(j+1,1))^2 + (xy_r(j,2)-xy_r(j+1,2))^2);
            end
        end
    end
    load('nm')
    Number_of_pixels= nm;
    perimeter_ratio=  d/graindata.Perimeter;
    % roundness
    roundness = mean(abs(radiusxy))/(radius_in);
    num_corn= size(CORNER,1)/2;
    load('Rq');
    NRq = Rq/nm * 100;
   
%% minimum circumscribing circle
% finding the maximum length of chord within the grain to find the maximum distance of boundary point from centroid
% 
x = xy_r(:,2); % x component of coordinate
y = xy_r(:,1); % y component of coordinate
% figure, imshow(bw)
% hold on
% plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
g= pdist2(xy_r,xy_r,'euclidean');
[nm,q] = max(g(:));
[ro co] = find(g== nm);
x3= xy_r(ro(1,1),2);
y3= xy_r(ro(1,1),1);
x4= xy_r(ro(2,1),2);
y4= xy_r(ro(2,1),1);
% h=impoint(gca, x3,y3)
% setColor(h,'g')
% h=impoint(gca, x4,y4)
% setColor(h,'g')
x1 = [x3 x4];
y1 = [y3 y4];
% line(x1,y1, 'Color','b','LineWidth',3)
xc= (x3+x4)/2;
yc= (y3+y4)/2;
% h=impoint(gca, xc,yc)
% setColor(h,'k')
rc= sqrt((x3 - xc)^2 + (y4 - yc)^2);
% viscircles([xc,yc],rc)
xnc= [x3;x4];
ync=[y3;y4];
xt=[];
yt=[];
dnc=[];
% loop to check if boundary points are outside the circle and including
% the farthest point in forming a circle until all boundary points are within the circle

for i=1: size(graindata.ConvexHull,1)
    for j=1: size(graindata.ConvexHull,1)
        dfc= sqrt((xc - graindata.ConvexHull(j,1))^2 + (yc - graindata.ConvexHull(j,2))^2);
        if dfc > rc
            xt= [xt; graindata.ConvexHull(j,1)];
            yt= [yt; graindata.ConvexHull(j,2)];
            dnc=[dnc; dfc];
        end
    end
    if size(xt,1) == 0
        break
    end
    [c2,ind2] = max(dnc(:,1));
    xnc= [xnc; xt(ind2,1)];
    ync= [ync; yt(ind2,1)];
    [xc1,yc1,Re,a] = circfit(xnc,ync); % function to fit a circle to a given set of points
    xc=xc1;
    yc=yc1;
    rc=Re;
    xt=[];
    yt=[];
    dnc=[];
end

syms xi yi
f(xi,yi)=(xi-xc)^2 + (yi-yc)^2 - rc^2;
h=ezplot(f,[-1000,4000,0,4000]);
coor_cir = get(h,'contourMatrix');
set(h, 'Color', 'g', 'LineWidth',2);
radius_cir= rc;
%title('final corners with circles, inscribed and circumscribed circles and length and width')

%% SPHERICITY CALCULATION by roman(2015)

%Area sphericity
Sa = (graindata.Area)/(pi*radius_cir*radius_cir);

%width to length ratio sphericity
tic
% figure, imshow(bw)
% hold on
% plot(xy_r(:, 2),xy_r(:, 1), 'r', 'LineWidth', 2)
% hold on
% syms xi yi
% f(xi,yi)=(xi-xc)^2 + (yi-yc)^2 - rc^2;
% h=ezplot(f,[-1000,4000,0,4000]);
% coor_cir = get(h,'contourMatrix');
% set(h, 'Color', 'r');
% radius_cir= rc;
hold on
g= pdist2(xy_r,xy_r,'euclidean');
[nm,q] = max(g(:));
[ro co] = find(g== nm);
x3= xy_r(ro(1,1),2);
y3= xy_r(ro(1,1),1);
x4= xy_r(ro(2,1),2);
y4= xy_r(ro(2,1),1);
% h=impoint(gca, x3,y3)
% setColor(h,'g')
% h=impoint(gca, x4,y4)
% setColor(h,'g')
x1 = [x3 x4];
y1 = [y3 y4];
%line(x1,y1, 'Color','b','LineWidth',3)
len = sqrt((x3-x4)^2 + (y3-y4)^2);

m1 = (diff(y1)/diff(x1));
minv= -1/m1;
index = ro(1,1);
index1 = ro(2,1);

%to find breadth of the parle

xb=[];
yb=[];
v=[ (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0; (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0 ]; % v stores coordinates of the two extreme points
im=[ min(index,index1); max(index,index1)]; % index of extreme points are stored in im
dist_p=[];
for j= im(1) : im(2) % to run for each segment
    pt = [x(j),y(j),0];
    a = v(1,:)-v(2,:);%Vector
    b = pt-v(2,:);%Vector
    dist_p = [dist_p; x(j) y(j) (norm(cross(a,b)) / norm(a)) j];
end
[dsp,indh] = max(dist_p(:,3));
xb= [xb; dist_p(indh,1)];
yb= [yb; dist_p(indh,2)];
dist_p=[];
vn=[ (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0; (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0  ];
imn=[ max(index,index1); length(x)];
for j= imn(1) : imn(2) % to run for each segment
    pt = [x(j),y(j),0];
    a = vn(1,:)-vn(2,:);%Vector
    b = pt-vn(2,:);%Vector
    dist_p = [dist_p; x(j) y(j) (norm(cross(a,b)) / norm(a)) j];
end
vn=[ (xy_r(max(index,index1),2)),(xy_r(max(index,index1),1)),0; (xy_r(min(index,index1),2)),(xy_r(min(index,index1),1)),0  ];
imn=[ 1; min(index,index1)];
for j= imn(1) : imn(2) % to run for each segment
    pt = [x(j),y(j),0];
    a = vn(1,:)-vn(2,:);%Vector
    b = pt-vn(2,:);%Vector
    dist_p = [dist_p; x(j) y(j) (norm(cross(a,b)) / norm(a)) j];
end
[dsp,indh] = max(dist_p(:,3));
xb= [xb; dist_p(indh,1)];
yb= [yb; dist_p(indh,2)];
% h=impoint(gca, xb(1,1),yb(1,1))
% setColor(h,'g')
% h=impoint(gca, xb(2,1),yb(2,1))
% setColor(h,'g')
%equation of parallel line through the 2nd point

if abs(x3- x4) < abs(y3-y4)
    syms yo
    xo = (yo - yb(2,1) + m1* xb(2,1))/m1;
    wm_y= 0:1:min(size(BW_filled,1),size(BW_filled,2));
    wm_x= subs(xo, wm_y);
else
    syms xo
    yo = m1*(xo - xb(2,1)) + yb(2,1);
    wm_x= 0:1:max(size(BW_filled,1),size(BW_filled,2));
    wm_y= subs(yo, wm_x);
end
% plot(wm_x, wm_y, 'g')
% plot(wm_x, wm_y, 'g')

%to draw parallel lines at max distance points; 1st point

if abs(x3- x4) < abs(y3-y4)
    syms yo
    xo = (yo - yb(1,1) + m1* xb(1,1))/m1;
    wm_y= 0:1:min(size(BW_filled,1),size(BW_filled,2));
    wm_x= subs(xo, wm_y);
else
    syms xo
    yo = m1*(xo - xb(1,1)) + yb(1,1);
    wm_x= 0:1:max(size(BW_filled,1),size(BW_filled,2));
    wm_y= subs(yo, wm_x);
end
% plot(wm_x, wm_y, 'g')
% plot(wm_x, wm_y, 'g')
wm_x= wm_x.';
wm_y= wm_y.';

%to find distance of 2nd point from the vector line plotted

c=[];
c(1,:)= [xb(1,1) yb(1,1) 0];
c(2,:)= [wm_x(end-5,1) wm_y(end-5,1) 0]; % any two points on a parallel line
pt = [xb(2,1),yb(2,1),0];
a = c(1,:)-c(2,:);%Vector
b = pt-c(2,:);%Vector
bre= (norm(cross(a,b)) / norm(a));

%to plot the perpendicular distance of breadth

A = [-m1 1; -minv 1];
B = [ -m1* xb(2,1) + yb(2,1) ; -minv* xb(1,1)+ yb(1,1)];
Xi = linsolve(A,B);
%impoint(gca,Xi(1,1), Xi(2,1))
X2=[xb(1,1) Xi(1,1)];
Y2=[yb(1,1) Xi(2,1)];


%equation of perpendicular line through the 2nd point

if abs(x3- x4) > abs(y3-y4)
    syms yo
    xo = (yo - y3(1,1) + minv* x3(1,1))/minv;
    wm_y= 0:1:min(size(BW_filled,1),size(BW_filled,2));
    wm_x= subs(xo, wm_y);
else
    syms xo
    yo = minv*(xo - x3(1,1)) + y3(1,1);
    wm_x= 0:1:max(size(BW_filled,1),size(BW_filled,2));
    wm_y= subs(yo, wm_x);
end
% plot(wm_x, wm_y, 'g')
% plot(wm_x, wm_y, 'g')

%equation of perpendicular line through the 1st point

if abs(x3- x4) > abs(y3-y4)
    syms yo
    xo = (yo - y4(1,1) + minv* x4(1,1))/minv;
    wm_y= 0:1:min(size(BW_filled,1),size(BW_filled,2));
    wm_x= subs(xo, wm_y);
else
    syms xo
    yo = minv*(xo - x4(1,1)) + y4(1,1);
    wm_x= 0:1:max(size(BW_filled,1),size(BW_filled,2));
    wm_y= subs(yo, wm_x);
end
% plot(wm_x, wm_y, 'g')
% plot(wm_x, wm_y, 'g')
% line(X2,Y2, 'Color','b','LineWidth',3)
Swl = min(bre,len)/max(bre,len);
%title('width and length of particle for sphericity calculation')
toc
 load('Round_spheri_rough');
    Round_spheri_rough =[Round_spheri_rough; Swl nm Rq NRq num_corn roundness perimeter_ratio];
    save('Round_spheri_rough','Round_spheri_rough')

fid = fopen(strcat('results_', filename, '.txt'),'a');
if fid<0
    error('Couldn''t open results.txt for writing!');
end

fprintf(fid, '%g\n', roundness);
fprintf(fid, '%g\n', Swl);
fprintf(fid, '%g\n', NRq);

disp(roundness);
disp(Swl);
disp(NRq);    
% xbar = size(BW_filled,2)-260;   % to call the centroids stored  in graindata matrix
% ybar = 50;
% info_string = sprintf('%2.2f',roundness);
% text(xbar + 180 ,ybar,info_string,'Color','g',...
%     'FontSize',10,'FontWeight','bold');
% text(xbar,ybar,'Roundness = ','Color','g',...
%     'FontSize',10,'FontWeight','bold');
% ybar = 80;
% info_string = sprintf('%2.2f',Swl);
% text(xbar + 180,ybar,info_string,'Color','g',...
%     'FontSize',10,'FontWeight','bold');
% text(xbar,ybar,'Sphericity = ','Color','g',...
%     'FontSize',10,'FontWeight','bold');
% ybar = 120;
% info_string = sprintf('%2.2f',NRq);
% text(xbar + 180 ,ybar,info_string,'Color','g',...
%     'FontSize',10,'FontWeight','bold');
% text(xbar,ybar,'NRq = ','Color','g',...
%     'FontSize',10,'FontWeight','bold');
% ybar = 900;
% info_string = sprintf('%2.2f',Angularity);
% text(xbar + 400,ybar,info_string,'Color','r',...
%     'FontSize',14,'FontWeight','bold');
% text(xbar,ybar,'Angularity = ','Color','r',...
%     'FontSize',14,'FontWeight','bold');
export_fig final.png
load('P')
%saveas(h,sprintf('FIG1 %d.fig',P))
saveas(h,sprintf(strcat('output2_',filename, '.jpg'),P))
end


