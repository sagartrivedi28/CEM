function [map_in,map_on,xx,yy,xoffset,yoffset] = DieOffset(W_die,H_die,wsize,ax)
%DIEOFFSET Summary of this function goes here
%   Detailed explanation goes here

r=wsize/2;
rs=3;
notch=9;

%assign map size and die count
W_die_count=ceil((r-W_die/2)/W_die)*2+3;
H_die_count=ceil((r-H_die/2)/H_die)*2+3;
map_in=false(H_die_count,W_die_count);
map_on=false(H_die_count,W_die_count);

%x=[xl,xr,xr,xl]; y=[yd,yd,yu,yu]
xind=-(W_die_count-1)/2:(W_die_count-1)/2;
yind=-(H_die_count-1)/2:(H_die_count-1)/2;
[xx,yy]=meshgrid(xind,yind);
x=ones(length(yind),1)*(xind*W_die); 
x=x(:);
x=x*ones(1,4)+ones(length(x),1)*[-1,1,1,-1]*W_die/2;
y=(yind*H_die).'*ones(1,length(xind)); 
y=y(:);
y=y*ones(1,4)+ones(length(y),1)*[-1,-1,1,1]*H_die/2;

%Condition
%x.^2+y.^2<=(r-rs).^2 & y>=-r+notch
%in_v=sum(((x.^2+y.^2<=(r-rs)^2) & (y>=-r+notch))*1,2);
%Find Offset
x_step=100;
y_step=100;
xoffset=0;
yoffset=0;
count=0;
dist=0;

for i=-y_step:y_step
    for j=0:x_step
        xs=j*W_die/2/x_step;
        ys=i*H_die/2/y_step;
        v=sum((((x+xs).^2+(y+ys).^2<=(r-rs)^2) & ((y+ys)>=-r+notch))*1,2);
        in_v_tmp=v==4;
        on_v_tmp=v>=1 & v<=3;
        count_tmp=sum(in_v_tmp*1);
        dist_tmp=min(min((r-rs)-sqrt((x(repmat(in_v_tmp(:),4,1))+xs).^2+(y(repmat(in_v_tmp(:),4,1))+ys).^2)),...
                     min(y(repmat(in_v_tmp(:),4,1))+ys+r-notch));
        if count_tmp>count || (count_tmp==count && dist_tmp>dist)
            
           xoffset=xs; yoffset=ys;
           count=count_tmp;
           in_v=in_v_tmp;
           on_v=on_v_tmp;
           dist=dist_tmp;
           
        end
    end
end

%Draw Map
map_in(in_v)=true;
map_on(on_v)=true;
x_aug=[x,x(:,1),NaN*ones(size(x,1),1)];
y_aug=[y,y(:,1),NaN*ones(size(x,1),1)];
x_in=x_aug(in_v,:).'; x_in=x_in(:)+xoffset;
y_in=y_aug(in_v,:).'; y_in=y_in(:)+yoffset;
x_on=x_aug(on_v,:).'; x_on=x_on(:)+xoffset;
y_on=y_aug(on_v,:).'; y_on=y_on(:)+yoffset;

axes(ax);
plot(x_in,y_in,'b','Linewidth',2);
hold on;
plot(x_on,y_on,'r--','Linewidth',1);
theta=0:2*pi/360:2*pi;
xc=(r-rs)*cos(theta);
yc=(r-rs)*sin(theta);
plot(xc,yc,'k');
xyt=axis(gca); 
text(xyt(1),xyt(3),[' Die Count = ' num2str(sum(map_in(:)),'%.0f')],'HorizontalAlignment','left','VerticalAlignment','bottom');
text(xyt(1),xyt(4),[' OFFSET = (' num2str(xoffset,'%.1f') 'mm,' num2str(yoffset,'%.1f') 'mm)'],'HorizontalAlignment','left','VerticalAlignment','top');
text(xyt(2),xyt(4),[' CNR2EDGE = ' num2str(dist,'%.1f') 'mm'],'HorizontalAlignment','right','VerticalAlignment','top');
plot(xoffset,yoffset,'+');
end

