function [data1,data2] = Im_align(xoffset2,yoffset2,inmat1,inmat2)
%--------------INPUT----------------
% xoffset2 = row offset
% yoffset2 = column offset
% inmat1 = reference matrix
% inmat2 = aligned matrix
%--------------OUTPUT---------------
% data1 = aligned reference matrix
% data2 = aligned test matrix
 
if(xoffset2>=0 && yoffset2>=0)
    xoffset1 = abs(xoffset2);
    yoffset1 = abs(yoffset2);
    data1 = inmat1(xoffset1+1:end,yoffset1+1:end);
    data2 = inmat2(1:end-xoffset1,1:end-yoffset1);
    
end;
 
if(xoffset2<=0 && yoffset2<=0)
    xoffset1 = abs(xoffset2);
    yoffset1 = abs(yoffset2);
    data1 = inmat1(1:end-xoffset1,1:end-yoffset1);
    data2 = inmat2(xoffset1+1:end,yoffset1+1:end);
    
end;
 
if(xoffset2<=0 && yoffset2>=0)
    xoffset1 = abs(xoffset2);
    yoffset1 = abs(yoffset2);
    data1 = inmat1(1:end-xoffset1,yoffset1+1:end);
    data2 = inmat2(xoffset1+1:end,1:end-yoffset1);
    
end;
 
if(xoffset2>=0 && yoffset2<=0)
    xoffset1 = abs(xoffset2);
    yoffset1 = abs(yoffset2);
    data1 = inmat1(xoffset1+1:end,1:end-yoffset1);
    data2 = inmat2(1:end-xoffset1,yoffset1+1:end);
    
end;
 
end
 

