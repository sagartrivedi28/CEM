function [TCC4_ams, lm, TCC2_ams,slm] = Vector_TCCfreq(Paims,Jsource,Mypara,Flagfast)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mypara.Ngrid = TCC Grid Size
% Mypara.base_ord =  TCC basis order
% Mypara.Sampleup = Up-sampling for TCC
% Mypara.lambda = 193nm Optical Source
% Mypara.ams_pz = AIMS Pixel-Size
% Mypara.NA = AIMS tool NA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ngrid = 21;
Ngrid = Mypara.Ngrid;
Tam = [Ngrid Ngrid];
if Flagfast
    Paims{1,1} = imresize(Paims{1,1},Tam); Paims{1,2} = imresize(Paims{1,2},Tam);
    Paims{2,1} = imresize(Paims{2,1},Tam); Paims{2,2} = imresize(Paims{2,2},Tam);
    Paims{3,1} = imresize(Paims{3,1},Tam); Paims{3,2} = imresize(Paims{3,2},Tam);
    
    Jsource{1,1} = imresize(Jsource{1,1},Tam); Jsource{1,2} = imresize(Jsource{1,2},Tam);
    Jsource{2,1} = imresize(Jsource{2,1},Tam); Jsource{2,2} = imresize(Jsource{2,2},Tam);
end;



TCC_master = zeros(4*[(2*Tam(1)-1)^2 (2*Tam(2)-1)^2]);
tcc_size = size(TCC_master);
Nam = 2*Tam-1;
Pupil = Paims;
am = size(Pupil{1,1});
for ro=1:2
    for co=1:2
        TCCxx = zeros([(2*Tam-1) (2*Tam-1)]); TCCxy = zeros([(2*Tam-1) (2*Tam-1)]);
        TCCyx = zeros([(2*Tam-1) (2*Tam-1)]); TCCyy = zeros([(2*Tam-1) (2*Tam-1)]);
        
        Pt{1} = zeros(3*am-1); Pt{1}(am(1):2*am(1)-1,am(2):2*am(2)-1) = Pupil{1,ro};
        Pt{2} = zeros(3*am-1); Pt{2}(am(1):2*am(1)-1,am(2):2*am(2)-1) = Pupil{2,ro};
        Pt{3} = zeros(3*am-1); Pt{3}(am(1):2*am(1)-1,am(2):2*am(2)-1) = Pupil{3,ro};

        Pstar{1} = rot90(conj(Pupil{1,co}),2);
        Pstar{2} = rot90(conj(Pupil{2,co}),2);
        Pstar{3} = rot90(conj(Pupil{3,co}),2);
        

        for i=round(1:((2*am(1)-2)/(Nam(1)-1)):(2*am(1)-1))
            for j=round(1:((2*am(2)-2)/(Nam(2)-1)):(2*am(2)-1))
                
                n1 = round((i-1)./((2*am(1)-2)/(Nam(1)-1)))+1;
                n2 = round((j-1)./((2*am(2)-2)/(Nam(2)-1)))+1;
                TCCxx(n1,n2,:,:) = imresize(conv2(Pt{1}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{1,1}, Pstar{1})...
                    +conv2(Pt{2}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{1,1}, Pstar{2}) + conv2(Pt{3}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{1,1}, Pstar{3}),Nam,'nearest');
                
                TCCxy(n1,n2,:,:) = imresize(conv2(Pt{1}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{1,2}, Pstar{1})...
                    +conv2(Pt{2}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{1,2}, Pstar{2}) + conv2(Pt{3}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{1,2}, Pstar{3}),Nam,'nearest');
                
                TCCyx(n1,n2,:,:) = imresize(conv2(Pt{1}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{2,1}, Pstar{1})...
                    +conv2(Pt{2}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{2,1}, Pstar{2}) + conv2(Pt{3}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{2,1}, Pstar{3}),Nam,'nearest');
                
                TCCyy(n1,n2,:,:) = imresize(conv2(Pt{1}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{2,2}, Pstar{1})...
                    +conv2(Pt{2}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{2,2}, Pstar{2}) + conv2(Pt{3}(i:i+am(1)-1,j:j+am(2)-1).*Jsource{2,2}, Pstar{3}),Nam,'nearest');
            end;
        end;
        indx1 = ((ro-1)*(tcc_size(1)/2))+1; indx2 = ro*(tcc_size(1)/2);
        indy1 = ((co-1)*(tcc_size(2)/2))+1; indy2 = co*(tcc_size(2)/2);
        TCC_master(indx1:indx2,indy1:indy2) = [reshape(TCCxx,[(2*Tam(1)-1)^2 (2*Tam(2)-1)^2]) reshape(TCCxy,[(2*Tam(1)-1)^2 (2*Tam(2)-1)^2]); ...
                                reshape(TCCyx,[(2*Tam(1)-1)^2 (2*Tam(2)-1)^2]) reshape(TCCyy,[(2*Tam(1)-1)^2 (2*Tam(2)-1)^2])];  
        co
    end;
end;

% figure(); imshow(abs(TCC_master),[]); colormap jet;
% base_ord = 100;
base_ord = Mypara.base_ord;
[U1,S1,~] = svds(TCC_master,base_ord);
lm = diag(S1);

%%
% Sampleup = 1;
Sampleup = Mypara.Sampleup;
lembda = Mypara.lembda;
ams_pz = Mypara.ams_pz;
ams_grid = Mypara.ams_grid;
NA = Mypara.NA;

fact1 = (lembda/(4*ams_pz*(4*NA)));
Freqgrid = round(ams_grid/fact1);
original_sz = (2*Tam-1);
new_sz = round(Sampleup*original_sz*fact1);
% TCC_basis(:,:,k)=ifftshift(ifft2(fftshift(mytemp)));
% TCC_basisxx=[]; TCC_basisxy=[]; TCC_basisyx=[]; TCC_basisyy=[];
Ftemp = zeros(ams_grid);
Fcor = floor((ams_grid-Freqgrid)/2+1);
for k=1:base_ord
    TCC_freqxx{k} = imresize(reshape(U1(1:original_sz(1)^2,k),original_sz),Freqgrid);
    Ftemp(Fcor(1):Fcor(1)+Freqgrid(1)-1,Fcor(2):Fcor(2)+Freqgrid(2)-1) = TCC_freqxx{k};
    TCC_freqxx{k} = Ftemp;
    
    TCC_freqxy{k} = imresize(reshape(U1(1+original_sz(1)^2:2*original_sz(1)^2,k),original_sz),Freqgrid);
    Ftemp(Fcor(1):Fcor(1)+Freqgrid(1)-1,Fcor(2):Fcor(2)+Freqgrid(2)-1) = TCC_freqxy{k};
    TCC_freqxy{k} = Ftemp;
    
    TCC_freqyx{k} = imresize(reshape(U1(1+2*original_sz(1)^2:3*original_sz(1)^2,k),original_sz),Freqgrid);
    Ftemp(Fcor(1):Fcor(1)+Freqgrid(1)-1,Fcor(2):Fcor(2)+Freqgrid(2)-1) = TCC_freqyx{k};
    TCC_freqyx{k} = Ftemp;
    
    TCC_freqyy{k} = imresize(reshape(U1(1+3*original_sz(1)^2:4*original_sz(1)^2,k),original_sz),Freqgrid);
    Ftemp(Fcor(1):Fcor(1)+Freqgrid(1)-1,Fcor(2):Fcor(2)+Freqgrid(2)-1) = TCC_freqyy{k};
    TCC_freqyy{k} = Ftemp;
end;
TCC4_ams = {TCC_freqxx{:}; TCC_freqxy{:}; TCC_freqyx{:}; TCC_freqyy{:}};

%%
% Sampleup = 1;

rdim = [(2*Tam(1)-1)^2 (2*Tam(1)-1)^2 (2*Tam(1)-1)^2 (2*Tam(1)-1)^2];
cdim = [(2*Tam(2)-1)^2 (2*Tam(2)-1)^2 (2*Tam(2)-1)^2 (2*Tam(2)-1)^2];
STCC_master = mat2cell(TCC_master,rdim,cdim);
STCC_master = [STCC_master{1,1} STCC_master{1,4}; STCC_master{4,1} STCC_master{4,4}];

[SU1,SS1,~] = svds(STCC_master,base_ord);
slm = diag(SS1);

for k=1:base_ord
    STCC_freqxx{k} = imresize(reshape(SU1(1:original_sz(1)^2,k),original_sz),Freqgrid);
    Ftemp(Fcor(1):Fcor(1)+Freqgrid(1)-1,Fcor(2):Fcor(2)+Freqgrid(2)-1) = STCC_freqxx{k};
    STCC_freqxx{k} = Ftemp;
    
    STCC_freqxy{k} = imresize(reshape(SU1(1+original_sz(1)^2:2*original_sz(1)^2,k),original_sz),Freqgrid);
    Ftemp(Fcor(1):Fcor(1)+Freqgrid(1)-1,Fcor(2):Fcor(2)+Freqgrid(2)-1) = STCC_freqxy{k};
    STCC_freqxy{k} = Ftemp;
end;
TCC2_ams = {STCC_freqxx{:}; STCC_freqxy{:}};


end

