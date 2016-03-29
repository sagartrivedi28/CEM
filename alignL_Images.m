function [new_img,xoffset2,yoffset2] = alignL_Images(imgs,factx,facty,range1)

% Last image align by previous image

imgno = length(imgs)-1;
base_img = imgs{1};
am = size(base_img);
xoffset2 = zeros(imgno,1);
yoffset2 = zeros(imgno,1);
for i=2:imgno
    [yoffset2(i),xoffset2(i),~] = ImgRegister(base_img(1:round(end/factx),1:round(end/facty)),imgs{i}(1:round(end/factx),1:round(end/facty)),range1);
end;

xoffset2 = (xoffset2-1); xoffset2(1)=0; 
yoffset2 = (yoffset2-1); yoffset2(1)=0;

Xmax = max(xoffset2); Ymax = max(yoffset2);
Xmin = min(xoffset2); Ymin = min(yoffset2);

for i=1:imgno
    new_img{i} = imgs{i}(Xmax-xoffset2(i)+1:end-(xoffset2(i)-Xmin),Ymax-yoffset2(i)+1:end-(yoffset2(i)-Ymin));
end;

new_img{imgno+1} = imgs{imgno+1}(Xmax-xoffset2(imgno)+1:end-(xoffset2(imgno)-Xmin),Ymax-yoffset2(imgno)+1:end-(yoffset2(imgno)-Ymin));

% gds1 = handles.gds(sem1_xmax+1:N-abs(sem1_xmin),sem1_ymax+1:N-abs(sem1_ymin));
% sem1_align(:,:,i-1) = sem1_bw(abs(sem1_xoffset(i)-sem1_xmax)+1:N-abs(sem1_xoffset(i)-sem1_xmin), abs(sem1_yoffset(i)-sem1_ymax)+1:N-abs(sem1_yoffset(i)-sem1_ymin),i-1);

end

