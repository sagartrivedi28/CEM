function [Paims, Pscanner, P0_norm] = Vector_RPupil(Mypara)


% stepZ = 0.01;
stepZ = Mypara.stepZ;
p = Mypara.p;
NA = Mypara.NA;
n = Mypara.n;
nr = Mypara.nr;
kr = Mypara.kr;
nr = nr + kr*sqrt(-1);

lembda = Mypara.lembda;


[X,Y] = meshgrid(-1:stepZ:1,-1:stepZ:1);
X(X==0) = eps^0.25; Y(Y==0) = eps^0.25;
am = size(X);
%%%%%%%%%%%%%% Pupil Cut-off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P0 = strel('disk',floor(am(1)/2),0);
P0_norm = (P0.getnhood);%/sum(sum(P0.getnhood));
P = P0_norm; %.*Vrc_AIMS;
% figure(); imshow(P,[]); colormap jet;

%%%%%%%%%%%%%% Aberration Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p =0:24;

zpoly = zeros(length(p),1); zpoly(1) = 1;
[theta,r] = cart2pol(X,Y);
idx = r<=1;
zer_basis = zeros([size(X,1).^2 length(p)]);
zer_temp = zernfun2(p,r(idx),theta(idx),'norm');
for i=2:length(p)
    zer_temp(:,i) = (zer_temp(:,i)-min(zer_temp(:,i)))/(max(zer_temp(:,i))-min(zer_temp(:,i)));
end;
zer_basis(idx(:),:)=zer_temp;
zer_basis = reshape(zer_basis,[size(X) length(p)]);
figure();
for i=1:length(p)
    subplot(5,5,i); imshow(zer_basis(:,:,i),[]); colormap jet;
end;

Phi=zpoly(1)*zer_basis(:,:,1);
for i=2:length(p)
    Phi = Phi + zpoly(i)*zer_basis(:,:,i);
end;
Abb = mat2gray(exp(i*2*pi*Phi));
% Abb = (Abb-min(Abb(:)))./(max(Abb(:))-min(Abb(:)));

%%%%%%%%%%%%%% Jone's Pupil (Pauli Decomposition) %%%%%%%%%%%%%%%%%%%
Ar = [1 0 0 0]; Ai = [0 0 0 0];
Pc = Ar + i*Ai;
J{1,1} = (Pc(1)+Pc(2))*ones(size(X));
J{1,2} = (Pc(3)-i*Pc(4))*ones(size(X));
J{2,1} = (Pc(3)+i*Pc(4))*ones(size(X));
J{2,2} = (Pc(1)-Pc(2))*ones(size(X));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% AIMS Mode Pupil%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Directional Unit Vector k=(A,B,R) and k'=(Ai,Bi,Ri)%%%%%%%%%%%
M=450;
A = X.*NA; B = Y.*NA; R = (1-A.^2-B.^2).^0.5;
Ai = A/(M*n); Bi = B/(M*n); Ri = (1-Ai.^2-Bi.^2).^0.5;

%%%%%%%%%%%% Rotation & Projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Normal = (A.^2+B.^2).^0.5;
Rmat{1,1} = A./Normal;          
Rmat{1,2} = B./Normal;         
Rmat{2,1} = Rmat{1,2};
Rmat{2,2} = -Rmat{1,1};

Normali = (1-Ri.^2).^0.5;
PSmat{1,1} = (Ai.*Ri)./Normali;    
PSmat{1,2} = Bi./Normali;          
PSmat{2,1} = (Bi.*Ri)./Normali;    
PSmat{2,2} = -Ai./Normali;         
PSmat{3,1} = -Normali;
PSmat{3,2} = 0*Normali;

PRmat{1,1} = PSmat{1,1}.*Rmat{1,1} + PSmat{1,2}.*Rmat{2,1};
PRmat{1,2} = PSmat{1,1}.*Rmat{1,2} + PSmat{1,2}.*Rmat{2,2};
PRmat{2,1} = PSmat{2,1}.*Rmat{1,1} + PSmat{2,2}.*Rmat{2,1};
PRmat{2,2} = PSmat{2,1}.*Rmat{1,2} + PSmat{2,2}.*Rmat{2,2};
PRmat{3,1} = PSmat{3,1}.*Rmat{1,1} + PSmat{3,2}.*Rmat{2,1};
PRmat{3,2} = PSmat{3,1}.*Rmat{1,2} + PSmat{3,2}.*Rmat{2,2};


%%%%%%%%%%%%%% Radio-Metric Correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Radio_crr = sqrt(R./Ri);
figure(); imshow(P.*Radio_crr,[]); colormap jet; colorbar;
%%%%%%%%%%%%%% Defocus Wafer Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fstep = 0;
Dfocus = exp(i*2*pi*n*Ri*Fstep/lembda);

%%%%%%%%%%%%%%%%%%% Pupil Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Slr_term = Radio_crr.*Abb.*P; %%%%%%%% Correct %%%%%%%%%

Pupil{1,1} = double(P~=0).*((PRmat{1,1}.*J{1,1} + PRmat{1,2}.*J{2,1}).*Slr_term);  %Pupil{1,1} = double(P~=0).*(Pupil{1,1}-min(Pupil{1,1}(:)))./(max(Pupil{1,1}(:))-min(Pupil{1,1}(:)));
Pupil{1,2} = double(P~=0).*((PRmat{1,1}.*J{1,2} + PRmat{1,2}.*J{2,2}).*Slr_term);  %Pupil{1,2} = double(P~=0).*(Pupil{1,2}-min(Pupil{1,2}(:)))./(max(Pupil{1,2}(:))-min(Pupil{1,2}(:)));
Pupil{2,1} = double(P~=0).*((PRmat{2,1}.*J{1,1} + PRmat{2,2}.*J{2,1}).*Slr_term);  %Pupil{2,1} = double(P~=0).*(Pupil{2,1}-min(Pupil{2,1}(:)))./(max(Pupil{2,1}(:))-min(Pupil{2,1}(:)));
Pupil{2,2} = double(P~=0).*((PRmat{2,1}.*J{1,2} + PRmat{2,2}.*J{2,2}).*Slr_term);  %Pupil{2,2} = double(P~=0).*(Pupil{2,2}-min(Pupil{2,2}(:)))./(max(Pupil{2,2}(:))-min(Pupil{2,2}(:)));
Pupil{3,1} = double(P~=0).*((PRmat{3,1}.*J{1,1} + PRmat{3,2}.*J{2,1}).*Slr_term);  %Pupil{3,1} = double(P~=0).*(Pupil{3,1}-min(Pupil{3,1}(:)))./(max(Pupil{3,1}(:))-min(Pupil{3,1}(:)));
Pupil{3,2} = double(P~=0).*((PRmat{3,1}.*J{2,1} + PRmat{3,2}.*J{2,2}).*Slr_term);  %Pupil{3,2} = double(P~=0).*(Pupil{3,2}-min(Pupil{3,2}(:)))./(max(Pupil{3,2}(:))-min(Pupil{3,2}(:)));

%%%%%%%%%%%%%%%%%%% Pupil Normalization%%%%%%%%%%%%%%%%%%%%%%%

Pupil_mat = (cell2mat(Pupil));
Paims = Pupil;

figure(); 
imshow(Pupil_mat,[]); colormap jet; colorbar;

% figure(); 
% subplot 321; imshow((Pupil{1,1}),[]); colormap jet; colorbar;
% subplot 322; imshow((Pupil{1,2}),[]); colormap jet; colorbar;
% subplot 323; imshow((Pupil{2,1}),[]); colormap jet; colorbar;
% subplot 324; imshow((Pupil{2,2}),[]); colormap jet; colorbar;
% subplot 325; imshow((Pupil{3,1}),[]); colormap jet; colorbar;
% subplot 326; imshow((Pupil{3,2}),[]); colormap jet; colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Scanner Mode Pupil%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Directional Unit Vector k=(A,B,R), k'=(Ai,Bi,Ri) & kr = (Ar,Br,Rr)%%%%%%%%%%%
M=0.25;
A = X.*NA; B = Y.*NA; R = (1-A.^2-B.^2).^0.5;
Ai = A/(M*n); Bi = B/(M*n); Ri = (1-Ai.^2-Bi.^2).^0.5;

Ar = A/(M*nr); Br = B/(M*nr); Rr = (1-Ar.^2-Br.^2).^0.5;
Ptm = (n*Rr)./(nr*Ri);  % Ex Component: TM Mode
Pte = (n*Ri)./(nr*Rr);  % Ey component: TE Mode

TMx = (1+Ptm)/2 - ((1-Ptm).^2)./(2*(1+Ptm));
TEy = (1+Pte)/2 - ((1-Pte).^2)./(2*(1+Pte));
TMz = TMx.*((n*Ri)./(nr*Rr));

%%%%%%%%%%%% In Rotation & Projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Normal = (A.^2+B.^2).^0.5;
Rmat{1,1} = A./Normal;          
Rmat{1,2} = B./Normal;         
Rmat{2,1} = Rmat{1,2};
Rmat{2,2} = -Rmat{1,1};

Normali = (1-Ri.^2).^0.5;
PSmat{1,1} = TMx.*((Ai.*Ri)./Normali);    
PSmat{1,2} = TEy.*(Bi./Normali);          
PSmat{2,1} = TMx.*((Bi.*Ri)./Normali);    
PSmat{2,2} = TEy.*(-Ai./Normali);         
PSmat{3,1} = TMz.*(-Normali);
PSmat{3,2} = 0*Normali;

PRmat{1,1} = PSmat{1,1}.*Rmat{1,1} + PSmat{1,2}.*Rmat{2,1};
PRmat{1,2} = PSmat{1,1}.*Rmat{1,2} + PSmat{1,2}.*Rmat{2,2};
PRmat{2,1} = PSmat{2,1}.*Rmat{1,1} + PSmat{2,2}.*Rmat{2,1};
PRmat{2,2} = PSmat{2,1}.*Rmat{1,2} + PSmat{2,2}.*Rmat{2,2};
PRmat{3,1} = PSmat{3,1}.*Rmat{1,1} + PSmat{3,2}.*Rmat{2,1};
PRmat{3,2} = PSmat{3,1}.*Rmat{1,2} + PSmat{3,2}.*Rmat{2,2};


%%%%%%%%%%%%%% Radio-Metric Correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Radio_crr = sqrt(R./Rr);
figure(); imshow(P.*Radio_crr,[]); colormap jet; colorbar;
%%%%%%%%%%%%%% Defocus Wafer Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fstep = 0;
Dfocus = exp(sqrt(-1)*2*pi*nr*Rr*Fstep/lembda);

%%%%%%%%%%%%%%%%%%% Pupil Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Slr_term = (Radio_crr).*Abb.*P; % Radio_crr/2 for Normalization

Pupil{1,1} = double(P~=0).*((PRmat{1,1}.*J{1,1} + PRmat{1,2}.*J{2,1}).*Slr_term);  %Pupil{1,1} = double(P~=0).*(Pupil{1,1}-min(Pupil{1,1}(:)))./(max(Pupil{1,1}(:))-min(Pupil{1,1}(:)));
Pupil{1,2} = double(P~=0).*((PRmat{1,1}.*J{1,2} + PRmat{1,2}.*J{2,2}).*Slr_term);  %Pupil{1,2} = double(P~=0).*(Pupil{1,2}-min(Pupil{1,2}(:)))./(max(Pupil{1,2}(:))-min(Pupil{1,2}(:)));
Pupil{2,1} = double(P~=0).*((PRmat{2,1}.*J{1,1} + PRmat{2,2}.*J{2,1}).*Slr_term);  %Pupil{2,1} = double(P~=0).*(Pupil{2,1}-min(Pupil{2,1}(:)))./(max(Pupil{2,1}(:))-min(Pupil{2,1}(:)));
Pupil{2,2} = double(P~=0).*((PRmat{2,1}.*J{1,2} + PRmat{2,2}.*J{2,2}).*Slr_term);  %Pupil{2,2} = double(P~=0).*(Pupil{2,2}-min(Pupil{2,2}(:)))./(max(Pupil{2,2}(:))-min(Pupil{2,2}(:)));
Pupil{3,1} = double(P~=0).*((PRmat{3,1}.*J{1,1} + PRmat{3,2}.*J{2,1}).*Slr_term);  %Pupil{3,1} = double(P~=0).*(Pupil{3,1}-min(Pupil{3,1}(:)))./(max(Pupil{3,1}(:))-min(Pupil{3,1}(:)));
Pupil{3,2} = double(P~=0).*((PRmat{3,1}.*J{2,1} + PRmat{3,2}.*J{2,2}).*Slr_term);  %Pupil{3,2} = double(P~=0).*(Pupil{3,2}-min(Pupil{3,2}(:)))./(max(Pupil{3,2}(:))-min(Pupil{3,2}(:)));

%%%%%%%%%%%%%%%%%%% Pupil Normalization%%%%%%%%%%%%%%%%%%%%%%%

Pupil_mat = (cell2mat(Pupil));
Pscanner = Pupil;

figure(); 
imshow(Pupil_mat,[]); colormap jet; colorbar;

end

