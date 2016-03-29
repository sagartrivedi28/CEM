function [xoffset,yoffset,c] =ImgRegister(img,template,maxShift)
% template should be smaller than img
% find location of template in given img

% % img1 = imread('D:\Shared Drive Data\OPC DATA security-B\SEMV\ADI1\Sem\Archive\000001_RT_01_SemAdcDef1_EDF3C_45_RID_002463.tif');
% % img2 = imread('D:\Shared Drive Data\OPC DATA security-B\SEMV\ADI1\Sem\Archive\000002_RT_01_SemAdcDef1_EDF3C_45_RID_002489.tif');
% 
% % img = img1;
% % template = img2;

% img = imfilter(img,fspecial('gaussian',5,5),'symmetric');
% template = imfilter(template,fspecial('gaussian',5,5),'symmetric');

% imshow(templatecrop)
templatecrop = template;
% templatewin = [3*size(template,2)/8 3*size(template,1)/8 size(template,2)/8 size(template,1)/8];
% templatewin = round(templatewin);
% templatecrop = imcrop(template, templatewin);

% templatecrop = imfilter(templatecrop,fspecial('gaussian',3,0.5));
 
c = normxcorr2_general(templatecrop(:,:,1),img(:,:,1));

% cost modification
[h1,w1]=size(img(:,:,1));
% [h2,w2]=size(templatecrop(:,:,1));
% lh = 1:(h1+h2-1);                   lw = 1:(w1+w2-1);
[h2,w2] = size(template(:,:,1));
mnh = min(h1,h2);                    mnw=min(w1,w2);


[m,n]=size(c);
c2 = zeros(m,n);
for k1=1:m
    for k2=1:n
%         d = sqrt((k1-481)^2+(k2-481)^2)+481;
       d = sqrt((k1-m/2)^2+(k2-n/2)^2)*0.002*1.41+1;
       d=max(d,1+0.002*1.41*(mnh+mnw)/3)-0.002*1.41*(mnh+mnw)/3;
       c2(k1,k2)=d;%costmap(481,round(d));
    end
end
c = c./c2;%.*costmap2;
% cwin = floor([size(templatecrop,2)/2 size(templatecrop,1)/2 size(img,1) size(img,2)]);
cwin = floor([size(c,2)*(1-maxShift)/2 size(c,1)*(1-maxShift)/2 size(c,2)*maxShift size(c,1)*maxShift]);
c = imcrop(c,cwin);

c(abs(c(:))<=0.95*max(abs(c(:))))=0;
[Lm,grp]=bwlabeln(c~=0);
dist2cntr=zeros(grp,1);

x_cntr=(1+size(c,2))/2;
y_cntr=(1+size(c,1))/2;

for i=1:grp
    [sub_y,sub_x]=ind2sub(size(Lm),find(Lm==i));
    dist2cntr(i)=min(sqrt((sub_x-x_cntr).^2+(sub_y-y_cntr).^2));
end

[~,mingrp]=min(dist2cntr);
[Cx,listy]=max((Lm==mingrp).*c);
[~,listx]=max(Cx);
xpeak=listx;
ypeak=listy(listx);

% disp(['cost =' num2str(cmax)]);
% corr_offset = [(xpeak-ceil(size(templatecrop,2)/2))
%                (ypeak-ceil(size(templatecrop,1)/2))];
corr_offset = [(xpeak + floor(size(c2,2)*(1-maxShift)/2)-ceil(size(c2,2)/2)+1)
               (ypeak+ floor(size(c2,1)*(1-maxShift)/2)-ceil(size(c2,1)/2)+1)];
%   imshow(c,[]);
% offset = -round(templatewin(1:2)); 
offset =[-1 -1];
xoffset = corr_offset(1) + offset(1);
yoffset = corr_offset(2) + offset(2);
% disp(['x= ' num2str(xoffset) ' y= ' num2str(yoffset)]);
% [offset] = [yoffset xoffset];
% y =yoffset;
% x =xoffset;
% imagwin =  [xoffset+offset yoffset+offset size(template,2)-1 size(template,1)-1];
% imgAPI = imcrop(img,imagwin);

% [h,w]=size(img1);
% x1 = max(1,2+xoffset);
% y1 = max(1,2+yoffset);
% x2 = min(w,w+xoffset+1);
% y2 = min(h,h+yoffset+1);
% w1 = [x1 y1 x2-x1 y2-y1];
% 
% x1 = max(1,-xoffset);
% y1 = max(1,-yoffset);
% x2 = min(w,w-xoffset-1);
% y2 = min(h,h-yoffset-1);
% w2 = [x1 y1 x2-x1 y2-y1];
% 
% img1 = imcrop(img,w1);
% img2 = imcrop(template,w2);
% subplot(1,2,1);imshow(img1);
% subplot(1,2,2);imshow(img2);
end
