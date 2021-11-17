
function [registered,delta_x,delta_y]=register(ref_rgb,target_rgb,eta,lr,lt,n,m,nuclei)
%Input:
    % ref_rgb: reference image RGB
    % target_rgb: target image RGB
    % eta: normalized cross-correlation threshold
    % lr: patch size of reference image
    % lt: patch size of target image, lt>lr
    % n: the number of nearest landmarks for partitioning strong and weak landmarks.
    % m: the number of nearest landmarks for selecting landmarks from strong landmarks.
    % nuclei: 1, nuclei as landmarks; 0, Harris points as landmarks.
 % Output:
    % registered: registered iamge RGB
    % delta_x: transformation of x
    % delta_y: transformation of y
if nargin < 2
        error(message('MATLAB:plregister:NotEnoughInputs'));
end
switch nargin
    case 2
        eta=0.25;lr=200;lt=500;n=4; m=2;nuclei=1;
    case 3
        lr=200;lt=500;n=4; m=2;nuclei=1;
    case 4
        lt=500;n=4; m=2;nuclei=1;
    case 5
        n=4; m=2;nuclei=1;
    case 6
        m=2;nuclei=1;
    case 7
        nuclei=1;
end

ref= normalizeStaining(ref_rgb);
ref=mat2gray(ref);
ref=rgb2gray(ref);
if nuclei==1
    centroid=extract_centroid(ref_rgb,20);
    [x,y]=find(centroid>=1);
    position=[x,y];
else
    corners = detectHarrisFeatures(ref);
    position=corners.selectStrongest(2*round(sqrt(size_n(1)*size_n(2))));
    position=position.Location;
    position=round(position);
end
N=size(position,1);

disp('Pre-aligned image')
target = normalizeStaining(target_rgb);
target=mat2gray(target);
target=rgb2gray(target);
optimizer = registration.optimizer.OnePlusOneEvolutionary;
optimizer.InitialRadius=1e-04;
metric = registration.metric.MattesMutualInformation;
tform = imregtform(target, ref, 'similarity', optimizer, metric);
target = imwarp(target,tform,'OutputView',imref2d(size(ref)));
target_rgb(:,:,1)=imwarp(target_rgb(:,:,1),tform,'OutputView',imref2d(size(ref)));
target_rgb(:,:,2)=imwarp(target_rgb(:,:,2),tform,'OutputView',imref2d(size(ref)));
target_rgb(:,:,3)=imwarp(target_rgb(:,:,3),tform,'OutputView',imref2d(size(ref)));
size_n=size(ref_rgb);
registered=zeros(size_n);
size_n=size_n(1:2);
delta_x=NaN(N);
delta_y=NaN(N);
position2=NaN(N,2);
delta1=round(lr/2);
delta2=round(lt/2);

disp('Matching landmarks by NCC')
for i=1:N
    try
        pos=position(i,1:2);
        adj_x=0;
        adj_y=0;
        adj_x1=0;
        adj_y1=0;
        x11=pos(1)-delta1;
        x12=pos(1)+delta1;
        y11=pos(2)-delta1;
        y12=pos(2)+delta1;
        x21=pos(1)-delta2;
        x22=pos(1)+delta2;
        y21=pos(2)-delta2;
        y22=pos(2)+delta2;
        if x11<1
            adj_x1=delta2-delta1;
            x11=1;
            x12=x11+2*delta1;
        end
        if x12>size(ref,1)
            adj_x1=delta1-delta2;
            x12=size(ref,1);
            x11=x12-2*delta1;
            
        end
        if y11<1
            adj_y1=delta2-delta1;
            y11=1;
            y12=y11+2*delta1;
            
        end
        if y12>size(ref,2)
           adj_y1=delta1-delta2;
           y12=size(ref,2);
           y11=y12-2*delta1;
           
        end
        if x21<1
            adj_x=delta2-pos(1);
            adj_x1=0;
            x21=1;
            x22=x21+2*delta2;
            
        end
        if x22>size(ref,1)
            adj_x=pos(1)-delta2;
            adj_x1=0;
            x22=size(ref,1);
            x21=x22-2*delta2;
            
        end
        if y21<1
            adj_y=delta2-pos(2);
            adj_y1=0;
            y21=1;
            y22=y21+2*delta2;
            
        end
        if y22>size(ref,2)
            adj_y=pos(2)-delta2;
            adj_y1=0;
            y22=size(ref,2);
            y21=y22-2*delta2;
            
        end
        x1=x11:x12;
        y1=y11:y12;
        x2=x21:x22;
        y2=y21:y22;
        ref_cent=ref(x1,y1);
        target_cent=target(x2,y2);
        cor=normxcorr2(ref_cent,target_cent);
        cor=cor(2*delta1:2*delta2,2*delta1:2*delta2);
        max_cor=max(cor(:));
        if max_cor>=eta
            [xpeak,ypeak] = find(cor==max(cor(:)));
            xpeak=xpeak(1);
            ypeak=ypeak(1);
            delta_x(i)=xpeak-delta1-delta2+2*delta1-1+adj_x1+adj_x;
            delta_y(i)=ypeak-delta1-delta2+2*delta1-1+adj_y1+adj_y;
            position2(i,:)=[pos(1)+delta_x(i),pos(2)+delta_y(i)];
        end
    catch ErrorInfo
      disp(ErrorInfo)
    end
end

disp('Partitioning landmarks into strong and weak')
z=find(~isnan(delta_x));
delta_x=delta_x(z);
delta_y=delta_y(z);
position=position(z,:);
position2=position2(z,:);
d1=pdist(position);
d1=squareform(d1);
min_n=n+1;
diff_l=0.1;
max_diff=1+diff_l;
min_diff=1-diff_l;
k=1;
clear z;
for i=1:size(d1,1)
    t=sort(d1(i,:));
    ind=find(d1(i,:)<=t(min_n));
    ind=setdiff(ind,i);
    xy1=position(i,:);
    xy2=position2(i,:);
    xy_d1=position(ind,:);
    xy_d2=position2(ind,:);
    flag=1;
    for j=1:(min_n-1)
        diff1=xy1(1)-xy_d1(j,1);
        diff2=xy2(1)-xy_d2(j,1);
        diff3=xy1(2)-xy_d1(j,2);
        diff4=xy2(2)-xy_d2(j,2);
        if diff1==0
            diff1=sign(diff2);
        end
        if diff3==0
            diff3=sign(diff4);
        end
        if sign(diff1)~=sign(diff2) || sign(diff3)~=sign(diff4)
            flag=0;
            break;
        end
        diff1=abs(diff1);
        diff2=abs(diff2);
        diff3=abs(diff3);
        diff4=abs(diff4);
        if diff2<min_diff*diff1 || diff2>max_diff*diff1 || diff4<min_diff*diff3 || diff4>max_diff*diff3
            flag=0;
           break;
        end
    end
    if flag==1
        z(k)=i;
        k=k+1;
    end
end
position_opt=position(z,:);
position2_opt=position2(z,:);
d1=pdist(position_opt);
d1=squareform(d1);
min_n=m+1;
diff_l=0.1;
max_diff=1+diff_l;
min_diff=1-diff_l;
k=1;
clear z;
for i=1:size(d1,1)
    t=sort(d1(i,:));
    ind=find(d1(i,:)<=t(min_n));
    ind=setdiff(ind,i);
    xy1=position_opt(i,:);
    xy2=position2_opt(i,:);
    xy_d1=position_opt(ind,:);
    xy_d2=position2_opt(ind,:);
    flag=1;
    for j=1:(min_n-1)
        diff1=xy1(1)-xy_d1(j,1);
        diff2=xy2(1)-xy_d2(j,1);
        diff3=xy1(2)-xy_d1(j,2);
        diff4=xy2(2)-xy_d2(j,2);
        if diff1==0
            diff1=sign(diff2);
        end
        if diff3==0
            diff3=sign(diff4);
        end
        if sign(diff1)~=sign(diff2) || sign(diff3)~=sign(diff4)
            flag=0;
            break;
        end
        diff1=abs(diff1);
        diff2=abs(diff2);
        diff3=abs(diff3);
        diff4=abs(diff4);
        if diff2<min_diff*diff1 || diff2>max_diff*diff1 || diff4<min_diff*diff3 || diff4>max_diff*diff3
            flag=0;
           break;
        end
    end
    if flag==1
        z(k)=i;
        k=k+1;
    end
end
position_opt=position_opt(z,:);
position2_opt=position2_opt(z,:);

disp('Spatial shift modification')
max_dist=0;
min_dist=sqrt((size_n(1)-1)^2+(size_n(2)-1)^2);
n=size(delta_x);
for i=1:n
    pq=position(i,:);
    [~,dist]=dsearchn(position_opt,pq);
    dist=dist(1);
    if dist==0
        continue
    end
    if dist>max_dist
        max_dist=dist;
    end
    if dist<min_dist
        min_dist=dist;
    end
end
z1=[];
j=1;
for i=1:n
    pq=position(i,:);
    [k,dist]=dsearchn(position_opt,pq);
    k=k(1);
    dist=dist(1);
    if dist==0
        continue;
    end
    w=(dist-min_dist)/(max_dist-min_dist);
    w=max(0,exp(-w)-0.5);
    if dist>0
        p1=position_opt(k,:);
        p2=position2_opt(k,:);
        delta_x(i)=w*(p2(1)-p1(1))+(1-w)*delta_x(i);
        delta_y(i)=w*(p2(2)-p1(2))+(1-w)*delta_y(i);
        if pq(1)+delta_x(i)<1
            delta_x(i)=1-pq(1);
        end
        if pq(1)+delta_x(i)>size_n(1)
            delta_x(i)=size_n(1)-pq(1);
        end
       if pq(2)+delta_y(i)<1
            delta_y(i)=1-pq(2);
        end
        if pq(2)+delta_y(i)>size_n(2)
            delta_y(i)=size_n(2)-pq(2);
        end
    else
        z1(j)=i;
        j=j+1;
    end
end
z1=setdiff((1:n),z1);
delta_x=delta_x(z1);
delta_y=delta_y(z1);
position=position(z1,:);
position2=position2(z1,:);
diff_l=0.1;
max_diff=1+diff_l;
min_diff=1-diff_l;
j=1;
for i=1:size(delta_x)
    pq=position(i,:);
    pq2=position2(i,:);
    [k,dist]=dsearchn(position_opt,pq);
    k=k(1);
    dist=dist(1);
    if dist==0
        z2(j)=i;
        j=j+1;
        continue;
    end
    p1=position_opt(k,:);
    p2=position2_opt(k,:);
    diff1=pq(1)-p1(1);
    diff2=pq2(1)-p2(1);
    diff3=pq(2)-p1(2);
    diff4=pq2(2)-p2(2);
    if diff1==0
       diff1=sign(diff2);
    end
    if diff3==0
        diff3=sign(diff4);
    end
    if sign(diff1)~=sign(diff2) || sign(diff3)~=sign(diff4)
        continue;
    end
    diff1=abs(diff1);
    diff2=abs(diff2);
    diff3=abs(diff3);
    diff4=abs(diff4);
    if diff2<min_diff*diff1 || diff2>max_diff*diff1 || diff4<min_diff*diff3 || diff4>max_diff*diff3
       continue;
    end
    z2(j)=i;
    j=j+1;
end
delta_x=delta_x(z2);
delta_y=delta_y(z2);
position=position(z2,:);
position2=position2(z2,:);

disp('Selection by texture and spatial proximity')
n=size(delta_x);
min1=10000;
max1=0;
for i=1:size(delta_x)
    dd=position;
    dd(i,:)=[];
    pq=position(i,:);
    [~,dist]=dsearchn(dd,pq);
    if dist>max1
        max1=dist;
    end
    if dist<min1
        min1=dist;
    end
end
j=1;
for i=1:n
    flag=1;
    pq=position(i,:);
    dd=position;
    dd(i,:)=[];
    [~,dist]=dsearchn(dd,pq);
    thresh2=0.1*exp(-1/dist);
    pq2=position2(i,:);
    x11=pq(1)-100;
    x12=pq(1)+100;
    y11=pq(2)-100;
    y12=pq(2)+100;
    x21=pq(1)-100;
    x22=pq(1)+100;
    y21=pq(2)-100;
    y22=pq(2)+100;
    if x11<1 || x12>size_n(1) || y11<1 || y12>size_n(2) || x21<1 || x22>size_n(1) || y21<1 || y22>size_n(2)
        z3(j)=i;
        j=j+1;
        continue
    end
    I1=ref(x11:x12,y11:y12);
    I2=target(x21:x22,y21:y22);
    points1 = detectSURFFeatures(I1);
    points2 = detectSURFFeatures(I2);
    f1 = extractFeatures(I1,points1);
    f2 = extractFeatures(I2,points2);
    [~,matchmetric] = matchFeatures(f1,f2);
    if size(matchmetric,1)==0
        z3(j)=i;
        j=j+1;
        continue
    end
    for ni=matchmetric'
        if ni>thresh2
            flag=0;
            break
        end
    end
    if flag==0
        continue
    end
    z3(j)=i;
    j=j+1;
end
delta_x=delta_x(z3);
delta_y=delta_y(z3);
position=position(z3,:);
x = position(:,1);
y = position(:,2);
n1=size_n(1);
n2=size_n(2);
maxn=max(n1,n2);

disp('generate registered image');
f = scatteredInterpolant(x,y,delta_x,'natural');
[xq,yq] = meshgrid(1:maxn,1:maxn);
delta_x = f(xq,yq);
f = scatteredInterpolant(x,y,delta_y,'natural');
delta_y = f(xq,yq);
delta_y1=yq+delta_y;
delta_x1=xq+delta_x;
delta_y1=round(delta_y1);
delta_x1=round(delta_x1);
for i=1:maxn
    for j=1:maxn
        if delta_x1(i,j)<=n1 && delta_x1(i,j)>=1 && delta_y1(i,j)>=1 && delta_y1(i,j)<=n2
            registered(j,i,:)=target_rgb(delta_x1(i,j),delta_y1(i,j),:);
        end
    end
end
registered=registered(1:n1,1:n2,:);
registered=uint8(registered);
fprintf('finished\n')
end
