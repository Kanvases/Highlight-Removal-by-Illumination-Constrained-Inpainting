clear;
addpath('../L-BFGS-B-C-master/Matlab')
dataPath='../data/1/';
I0=double(imread([dataPath,'o.png']))./255;
trimap=double(rgb2gray(imread([dataPath,'trimap.png'])));

% specular color in RGB
Cs=[0.3,0.3,0.3];
% stopping function threshold
t=0.1;
% smooth 
gamma=1e-4;

ROI=zeros(size(trimap));
ROI(trimap>0)=1;

% inital alpha
It=double(imread('../data/1/truth.png'));
I=double(imread('../data/1/o.png'));
ta=sum(I-It,3)./(3*255);
aph=ta(ROI==1);

% aph=rand([sum(ROI(:)),1]);
% aph=zeros([sum(ROI(:)),1]);

% ---- energy function ----- %
% --- N, edge stopping coff. ---%
% llumination constraints line in each pixels
% in RGB space
% Cs=(n,p,m)
% I0=(x0,y0,z0)
% x=nt+x0; y=pt+y0; z=mt+z0;
LCL=I0;
% for convinece, make the point t=1 as the LCL
LCL(:,:,1)=I0(:,:,1)+Cs(1)*1;
LCL(:,:,2)=I0(:,:,2)+Cs(2)*1;
LCL(:,:,3)=I0(:,:,3)+Cs(3)*1;

% 3 points for a plane.
% (0,0,0); (x0,y0,z0) (nt+x0,pt+y0,mt+z0)
N=fa_vector(zeros(size(I0)),I0,LCL);
[N1x,N1y]=(gradient(N(:,:,1)));
[N2x,N2y]=(gradient(N(:,:,2)));
[N3x,N3y]=(gradient(N(:,:,3)));

del_Nx_L1=sum(cat(3,abs(N1x),abs(N2x),abs(N3x)),3);
del_Ny_L1=sum(cat(3,abs(N1y),abs(N2y),abs(N3y)),3);
del_Nx_L1(isnan(del_Nx_L1))=0;
del_Ny_L1(isnan(del_Ny_L1))=0;
clear N1x N1y N2x N2y N3x N3y




Ep.sx=edgestop(del_Nx_L1,t);
Ep.sy=edgestop(del_Ny_L1,t);
Ep.I0=I0;
Ep.Cs=Cs;
Ep.gamma=gamma;

E=energy(aph,Ep,ROI)
