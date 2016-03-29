% clc; clear all; close all

%% Paramerters
stepZ = 0.01;                       % Pupil Resolution
Sann = false;                        % Annular Light-Source Flag
Dout = 0.86;                        % Annular outer diameter
Din = 0.86*0.85;                    % Annular inner diameter

NA=0.35;                            % For Wafer NA = 1.4 or 1.35;
Sigma = 1;                          % Source Bandwidth
lembda=193;                         % unit is "nm"
n=1.44;                             % Refractive Index Mask=1, Wafer=1.44;
p = 0:24;                           % Zernike Polynomial Length;

ams_pz = 17.488;                    % unit is "nm"
ams_grid = [568,568];               % AIMS image Size
ams_fov = ams_pz*ams_grid(1);       % AIMS FOV
sem_fov = 4500;                     % SEM FOV unit is "nm"
Ngrid = 21;                         % TCC Grid
base_ord = 100;                     % Eigen Vector order
Sampleup = 1;                       % Eigen Vector Up-Sampling

Mypara.stepZ = stepZ;
Mypara.Sann = Sann;
Mypara.Dout = Dout;
Mypara.Din = Din;
Mypara.NA = NA;
Mypara.Sigma = Sigma;
Mypara.lembda = lembda;
Mypara.n = n;
Mypara.p = p;
Mypara.ams_pz = ams_pz;
Mypara.ams_grid = ams_grid;
Mypara.ams_fov = ams_fov;
Mypara.sem_fov = sem_fov;
Mypara.Ngrid = Ngrid;
Mypara.base_ord =  base_ord;
Mypara.Sampleup = Sampleup;

% Not used parameters
% M=0.25;                             % Magnification AIMS = 450 , Scanner = 1/4;


%% Image Path

sr_path1 = 'D:\SAGAR\WORK SPACE\SAGAR\AIMS\Main Project\AIMS simulation project\Production cases(N16)\Source File';                   % Light Source File
sr_name1 = 'HA1128_R0 to ZEISS.src';



sr_path2 = 'D:\SAGAR\WORK SPACE\SAGAR\AIMS\Main Project\AIMS simulation project\Production cases(N16)\GDS\GDS_img';  % OPC GDS
sr_name2 = '10.png';
tlist = dir(fullfile(sr_path2,'*.png'));
gdspath = sr_path2;
gdslist = {tlist.name};

sr_path3 = 'D:\SAGAR\WORK SPACE\SAGAR\AIMS\Main Project\AIMS simulation project\Production cases(N16)\SEM Images\SEM_B';           % SEM Image
sr_name3 = '10-ST.tif'; 
tlist = dir(fullfile(sr_path3,'*.tif'));
sempath = sr_path3;
semlist = {tlist.name};

sr_path = 'D:\SAGAR\WORK SPACE\SAGAR\AIMS\Main Project\AIMS simulation project\Production cases(N16)\AIMS_mode\AIMS_B'; % AIMS Mode image
sr_name = '10_ST.csv';
tlist = dir(fullfile(sr_path,'*.csv'));
Aamspath = sr_path;
Aamslist = {tlist.name};

sr_pathS = 'D:\SAGAR\WORK SPACE\SAGAR\AIMS\Main Project\AIMS simulation project\Production cases(N16)\Scanner_mode\Scanner_B'; % Scanner Mode image
sr_nameS = '10_ST.csv';
tlist = dir(fullfile(sr_path,'*.csv'));
Samspath = sr_pathS;
Samslist = {tlist.name};

%% Generate Pupil & Light Source
[Paims, Pscanner, P0_norm] = Vector_Pupil(Mypara);
am = size(P0_norm);
if ~Sann
    %%%%%%%%%%%%%% Light Source File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mysrf = fopen(fullfile(sr_path1,sr_name1),'r');
    my_source = textscan(mysrf,'%s','delimiter','\n');
    my_source=my_source{1,1};
    ind1 = find(~cellfun(@isempty,strfind(my_source,'[DATA]')))+1;
    Sdata = cell2mat(cellfun(@str2num,my_source(ind1:end), 'UniformOutput', false));
    tempZ = abs(Sdata(1,1)-Sdata(2,1));
    Smat = reshape(Sdata(:,3),sqrt(length(Sdata(:,3))),sqrt(length(Sdata(:,3))))';
    % Smat = imresize(Smat,[21 21]);
%     am = size(Smat);
    Smat = imresize(Smat,am);
    figure(); imshow(Smat); colormap(jet);
else
    %%%%%%%%%%%%%% Light Source Annular %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Smat = imerode(P0_norm,strel('disk',floor(am(1)*(1-Dout)/2),0)) & (~imerode(P0_norm,strel('disk',floor(am(1)*(1-Din)/2),0)));
    figure(); imshow(Smat,[]);
end;

% am = size(X); Smat = imresize(Smat,am,'nearest')
Ypol = xor(triu(true(size(Smat))),rot90(triu(true(size(Smat))))); 
Exs = Smat.*double(~Ypol);
Eys = Smat.*double(Ypol);

figure(); subplot 121; imshow(Exs,[]); subplot 122; imshow(Eys,[]); colormap jet;
figure(); imshow(Exs | Eys,[]); colormap jet;
Jsource{1,1} = Exs;
Jsource{1,2} = 0*Exs;
Jsource{2,1} = 0*Eys;
Jsource{2,2} = Eys;
close all;
%% Grid Alignment & TCC Generation
% [TCC4_ams, lm_ams, TCC2_ams,slm_ams] = Vector_TCCfreq(Paims,Jsource,Mypara,false);
ATCC_freqxx = {TCC2_ams{1,:}}; ATCC_freqxy = {TCC2_ams{2,:}}; 
% ATCC_freqyx = {TCC4_ams{3,:}}; ATCC_freqyy = {TCC4_ams{4,:}}; 
Alm = slm_ams;
% save('NFreqMxTCC_ams.mat','TCC4_ams','lm_ams','TCC2_ams','slm_ams');

% [TCC4_scanner, lm_scanner, TCC2_scanner,slm_scanner] = Vector_TCCfreq(Pscanner,Jsource,Mypara,false);
STCC_freqxx = {TCC2_scanner{1,:}}; STCC_freqxy = {TCC2_scanner{2,:}}; 
% STCC_freqyx = {TCC4_scanner{3,:}}; STCC_freqyy = {TCC4_scanner{4,:}}; 
Slm = slm_scanner;
% save('NFreqMxTCC_scaneer.mat','TCC4_scanner','lm_scanner','TCC2_scanner','slm_scanner');

%% Mask Transmission Function
mask_ts = imread(fullfile(gdspath,gdslist{1}));
mask_ts = ~(mask_ts(:,:,1)~=0); %figure; imshow(mask_ts,[]);

mask_sem = imread(fullfile(sempath,semlist{1})); 
mask_sem = im2double(mask_sem(:,:,1)); figure; imshow(mask_sem,[]);

sz = size(mask_ts);
if sz(1)>sz(2)
    mask_ts = mask_ts(1:sz(2),:);
else
    mask_ts = mask_ts(:,1:sz(1));
end;
figure(); imshow(mask_ts);

%%%%%%%%%%%%%%%%%%% Mask Modeling: Bias + Cr Rounding %%%%%%%%
[ctr_out,gds_bias,close_rounding,open_rounding] = ctr_match(mask_ts,mask_sem,ams_fov,sem_fov,false);
figure(); imshow(ctr_out);

%%%%%%%%% Anti-Aliasing Filtering %%%%%%%%%%%%%%%%%
nz = size(ctr_out);
mask_pz = ams_fov/nz(1);
FGsigma =  2*pi*((4*(4*NA)/lembda)+(1/ams_pz))/2; 
TGsigma = 1/(mask_pz*FGsigma);
AA = fspecial('gaussian',round(TGsigma*10),TGsigma);
close all;
%% TCC Normalization
% 
% for tf=1:100
%     ATCC_freqfinal{tf}=ifftshift(ATCC_freqxx{tf} + ATCC_freqxy{tf});% + ATCC_freqyx{tf} + ATCC_freqyy{tf});    
%     STCC_freqfinal{tf}=ifftshift(STCC_freqxx{tf} + STCC_freqxy{tf});% + STCC_freqyx{tf} + STCC_freqyy{tf});
% end;

Sord = 42;
Alm = slm_ams(1:Sord)/1500;
Slm = slm_scanner(1:Sord)/1500;
% ATCCweight=0;
% for k=1:Sord
%     ATCCweight = ATCCweight + Alm(k)*sum(sum(abs(ATCC_freqfinal{k}))).^2;
% end;
% Alm = (Alm(1:Sord)/sum(Alm(1:Sord)));

% STCCweight=0;
% for k=1:(2*Sord)
%     STCCweight = STCCweight + Slm(k)*sum(sum(abs(STCC_freqfinal{k}))).^2;
% end;
% Slm = (Slm(1:(2*Sord))/sum(Slm(1:2*Sord))); % Scanner Gain

%% Reading Files

for No=1
    close all;
%%%%%%%%%%%%%%%% GDS Read & Mask Model %%%%%%%%%%%%%%%%%
mask_ts = imread(fullfile(gdspath,gdslist{No}));
mask_ts = ~(mask_ts(:,:,1)~=0); %figure; imshow(mask_ts,[]);
sz = size(mask_ts);
if sz(1)>sz(2)
    mask_ts = mask_ts(1:sz(2),:);
else
    mask_ts = mask_ts(:,1:sz(1));
end;
figure(1); imshow(mask_ts);
Smask_ts = imresize(mask_ts,nz,'nearest');
if(gds_bias<0)
    Smask_ts = imerode(Smask_ts,strel('square',-(2*gds_bias-1)));        
end;    
if(gds_bias>=0)
    Smask_ts = imdilate(Smask_ts,strel('square',2*gds_bias+1));        
end;
ctr_out = imclose(imopen(Smask_ts,strel('disk',open_rounding,0)),strel('disk',close_rounding,0));

Nctr_out = conv2(double(ctr_out),AA,'same');
mask_ts1 = imresize(Nctr_out,ams_grid,'nearest');
% mask_ts1 = (mask_ts1+0.06);
figure(1); imshow(mask_ts1,[]); colormap jet;


%%%%%%%%%%%%%%%%%%% AIMS Mode %%%%%%%%%%%%%%%%%%%%%%%%

mycsv = fopen(fullfile(Aamspath,Aamslist{No}),'r');
mydata1 = textscan(mycsv,'%f','delimiter',{',',';'});
mydata1=mydata1{1,1};
mydata1(isnan(mydata1))=[];
mydata = mydata1(7:end);
ams = reshape(mydata,Mypara.ams_grid(1),Mypara.ams_grid(2));
ams=flipud(rot90(ams,1));
AIMS_ams = ams;%mat2gray(ams);
figure(2); imshow(AIMS_ams,[]); colormap(jet);


%%%%%%%%%%%%%%%%%%% Scanner Mode %%%%%%%%%%%%%%%%%%%%%%%%
mycsv = fopen(fullfile(Samspath,Samslist{No}),'r');
mydata = textscan(mycsv,'%f','delimiter',{',',';'});
mydata=mydata{1,1};
mydata(isnan(mydata))=[];
mydata = mydata(7:end);
ams = reshape(mydata,Mypara.ams_grid(1),Mypara.ams_grid(2));
ams=flipud(rot90(ams,1));
Scanner_ams = ams;%mat2gray(ams);
figure(3); imshow((Scanner_ams),[]); colormap(jet); colorbar;

%% Initial Condition

G = fspecial('gaussian',21,5);
AMSscale = [min(AIMS_ams(:)) max(AIMS_ams(:))];
mask_ts2 = (conv2((double(mask_ts1)*AMSscale(2)+0.06),G,'same'));
Fmask_ts2 = (fft2(mask_ts2));
my_aimst = zeros(size(AIMS_ams));
for k=1:Sord
    my_aimst = my_aimst + Alm(k)*(ifft2(Fmask_ts2.*ifftshift(ATCC_freqxx{k} + ATCC_freqxy{k}))).^2;
end;

figure(4); imshow(abs(my_aimst),[]); colormap jet; colorbar;

AMS_temp = abs(my_aimst);
AMS_temp1 = AMS_temp;
figure(5); imshow(abs(AMS_temp1),[]); colormap jet; colorbar;


% [yoffset2,xoffset2,~] = ImgRegister(abs(AMS_temp1),AIMS_ams,0.5);
% xoffset2 = xoffset2-1; yoffset2 = yoffset2-1;
% [AMS_sim,AMS_org] = Im_align(xoffset2,yoffset2,AMS_temp1,AIMS_ams);
% [mask_ts3,~] = Im_align(xoffset2,yoffset2,mask_ts2,AIMS_ams);

imgs{1}=AIMS_ams; imgs{2}=Scanner_ams; imgs{3}=AMS_temp1; imgs{4}=mask_ts2;
[new_img,xoffset2,yoffset2] = alignL_Images(imgs,1,1,0.5);
AMS_org = new_img{1}; Scanner_org = new_img{2}; AMS_sim = new_img{3}; mask_ts3 = new_img{4};
Fmask_ts3 = fft2(mask_ts3);

Fsize = size(mask_ts3);
Fcor = floor((ams_grid-size(mask_ts3))/2+1);
for tf=1:base_ord
    ATCC_freqfinal{tf}=ifftshift(ATCC_freqxx{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1) + ...
                    ATCC_freqxy{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1));% + ...
                    %ATCC_freqyx{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1) + ...
                    %ATCC_freqyy{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1));
    STCC_freqfinal{tf}=ifftshift(STCC_freqxx{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1) + ...
                    STCC_freqxy{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1)); %+ ...
                    %STCC_freqyx{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1) + ...
                    %STCC_freqyy{tf}(Fcor(1):Fcor(1)+Fsize(1)-1,Fcor(2):Fcor(2)+Fsize(2)-1));
end;


figure(6); imshow(abs(AMS_sim(41:end-40,41:end-40)-AMS_org(41:end-40,41:end-40)),[]); colormap(jet); colorbar;


%% MNF Calibration

int=0;
old_sim = zeros(size(AMS_sim));
Region1 = zeros(size(AMS_sim));
Region1(11:end-10,11:end-10) = 1;
old_diff = 1;
% diff_gif = zeros(size(AMS_sim));
% AMS_gif = zeros(size(AMS_sim));

while((max(max(abs(AMS_sim(21:end-20,21:end-20)-AMS_org(21:end-20,21:end-20)))))>0.008)
    Mweight=zeros(size(mask_ts3));
    Icost = (AMS_sim-AMS_org);
%     diff_gif(:,:,int+1)=Icost;
    for k=1:Sord
        Tmfr = Icost.*(ifft2((Fmask_ts3).*(ATCC_freqfinal{k})));
        Mweight = Mweight + Alm(k)*((fft2(Tmfr).*conj(ATCC_freqfinal{k})));
    end;
    Mweight = ifft2(Mweight);
    fact = (0.05*max(abs(Icost(:))))/(max(abs(Mweight(:))));
    old_mask = mask_ts3;
    mask_ts3 = mask_ts3 - (fact*(Mweight));    
    Fmask_ts3 = fft2(mask_ts3);
    figure(5); imshow(mask_ts3,[]); colormap jet; colorbar;   
    
    AMS_sim=zeros(size(mask_ts3));
    for k=1:Sord
        AMS_sim = AMS_sim + (Alm(k))*(ifft2(Fmask_ts3.*ATCC_freqfinal{k})).^2;
    end;
    
%     if ~int
%         Inorm = max(max(abs(AMS_sim(11:end-10,11:end-10))));
%     end;  

    AMS_sim = abs(AMS_sim);%/Inorm;
%     AMS_gif(:,:,int+1) = (AMS_sim);
    
    figure(10); imshow((AMS_sim(21:end-20,21:end-20)-AMS_org(21:end-20,21:end-20)),[]); colormap jet; colorbar;
    int = int+1;
    Merr = max(max(abs(AMS_sim(21:end-20,21:end-20)-AMS_org(21:end-20,21:end-20))));
    figure(20); hold on;  plot(int,Merr,'*r');
    if ~mod(int,20)
        int=int;
    end;
    
end;

%%  Store Simulation Process
% for i=1:40
% if i == 1;
% imwrite(AMS_gif(:,:,i),colormap(jet),'SimulatedAMS','gif', 'Loopcount',Inf);
% else
% imwrite(AMS_gif(:,:,i),colormap(jet),'SimulatedAMS','gif','WriteMode','append');
% end
% end;
%% Vector Effect Embedding

Scanner_sim=zeros(size(mask_ts3));
for k=1:(Sord)
    Scanner_sim = Scanner_sim + (Slm(k))*(ifft2(Fmask_ts3.*STCC_freqfinal{k})).^2;
end;
Sgain = 2.2846;
Scanner_sim = (Sgain*abs(Scanner_sim));
Nsize = size(Scanner_sim(11:end-10,11:end-10));
figure(30); imshow(Scanner_sim,[]); colormap jet; colorbar;
Tempdiff = [abs(AMS_org(11:end-10,11:end-10)-Scanner_org(11:end-10,11:end-10)) abs(Scanner_sim(11:end-10,11:end-10)-Scanner_org(11:end-10,11:end-10))];
figure(40); imshow(Tempdiff,[]); colormap jet; colorbar;

%% Data Storing
mkdir(gdspath,'VectorImage');
save(fullfile(gdspath,'VectorImage',['MNR_' gdslist{No} '_.mat']),'mask_ts3');
MasterImg = [AMS_sim(11:end-10,11:end-10) AMS_org(11:end-10,11:end-10) Scanner_sim(11:end-10,11:end-10) Scanner_org(11:end-10,11:end-10)];
imwrite(uint8(round(255*abs(MasterImg))),colormap(jet(255)),fullfile(gdspath,'VectorImage',['MasterI_' gdslist{No} '_.bmp']));

MasterDiff = [abs(AMS_sim(11:end-10,11:end-10)-AMS_org(11:end-10,11:end-10)) abs(Scanner_sim(11:end-10,11:end-10)-Scanner_org(11:end-10,11:end-10)) ...
    abs(AMS_sim(11:end-10,11:end-10)-Scanner_sim(11:end-10,11:end-10)) abs(AMS_org(11:end-10,11:end-10)-Scanner_org(11:end-10,11:end-10))];
SV_OV = max(max(abs(Scanner_sim(11:end-10,11:end-10)-Scanner_org(11:end-10,11:end-10)))); SS_SV = max(max(abs(AMS_sim(11:end-10,11:end-10)-Scanner_sim(11:end-10,11:end-10))));
OS_OV = max(max(abs(AMS_org(11:end-10,11:end-10)-Scanner_org(11:end-10,11:end-10))));
imwrite(uint8(round(255*mat2gray(MasterDiff))),colormap(jet(255)),fullfile(gdspath,'VectorImage',['MasterD_' gdslist{No} num2str([SV_OV,SS_SV,OS_OV]) '_.bmp']));

mkdir(gdspath,'VectorCSV');
tempA = 0.5*ones(size(AIMS_ams)); Aind = round((size(AIMS_ams)-Nsize)/2);
tempA(Aind(1)+1:Aind(1)+Nsize(1),Aind(2)+1:Aind(2)+Nsize(2)) = AMS_sim(11:end-10,11:end-10);
dlmwrite(fullfile(gdspath,'VectorCSV',['AIMS_' Aamslist{No}]),mydata1(1:6)','delimiter',',');
dlmwrite(fullfile(gdspath,'VectorCSV',['AIMS_' Aamslist{No}]),tempA,'-append','delimiter',';');

tempA(Aind(1)+1:Aind(1)+Nsize(1),Aind(2)+1:Aind(2)+Nsize(2)) = Scanner_sim(11:end-10,11:end-10);
dlmwrite(fullfile(gdspath,'VectorCSV',['Scanner_' Samslist{No}]),mydata1(1:6)','delimiter',',');
dlmwrite(fullfile(gdspath,'VectorCSV',['Scanner_' Samslist{No}]),tempA,'-append','delimiter',';');

end;
