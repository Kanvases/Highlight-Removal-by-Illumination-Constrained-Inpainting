clear;
addpath('../L-BFGS-B-C-master/Matlab')
addpath('../tvinpaint')
dataPath='../data/f1/';
I0=double(imread([dataPath,'o.png']))./255;
trimap=double(rgb2gray(imread([dataPath,'t.png'])));

% specular color in RGB
Cs=[0.25,0.375,0.375];
% stopping function threshold
t=0.2;
% r g smooth 
gamma=100;

ROI=zeros(size(trimap));
ROI(trimap>0)=1;


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

% inital alpha
aph1=initAlpha(Ep,ROI);
aph=aph1(ROI==1);
imwrite(aph1,[dataPath,'init_alpha.png']);

% ---- optimization ----%
l=zeros(size(aph));
u=ones(size(aph));
fcn=@(a)energy(a,Ep,ROI);
grad=@(a)grandientEnergy(a,Ep,ROI);
fun=@(a)fminunc_wrapper(a,fcn,grad);
opts=struct('x0',aph,'maxIts',200);
[fa,~,~]=lbfgsb(fun,l,u,opts);

final_alpha=zeros(size(ROI));
final_alpha(ROI==1)=fa;
% figure;
imwrite(final_alpha,[dataPath,'alpha.png']);

Res(:,:,1)=I0(:,:,1)-final_alpha*Cs(1);
Res(:,:,2)=I0(:,:,2)-final_alpha*Cs(2);
Res(:,:,3)=I0(:,:,3)-final_alpha*Cs(3);
imwrite(uint8(Res.*255),[dataPath,'optRes.png']);

mask=ones(size(ROI));
mask(final_alpha>0.7)=0;
disp('tv inpainting...')
lam=1;
tao=0.1;
inP(:,:,1)=Inpainting_TV(Res(:,:,1),mask,tao,lam);
inP(:,:,2)=Inpainting_TV(Res(:,:,2),mask,tao,lam);
inP(:,:,3)=Inpainting_TV(Res(:,:,3),mask,tao,lam);
imwrite(uint8(inP.*255),[dataPath,'inp.png']);
mask3=repmat(mask,1,1,3);
Ifin=Res;
Ifin(mask3==0)=inP(mask3==0);
imwrite(uint8(Ifin.*255),[dataPath,'finalRes.png']);
