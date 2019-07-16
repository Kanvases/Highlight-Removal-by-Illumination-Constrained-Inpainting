function E=energy(aph,Ep,ROI)
% energy function calculation
alpha=zeros(size(ROI));
alpha(ROI==1)=aph;

% specular color * alpha = highlight effect in image(CS) 
CS=repmat(alpha,1,1,3);
CS(:,:,1)=CS(:,:,1)*Ep.Cs(1);
CS(:,:,2)=CS(:,:,2)*Ep.Cs(2);
CS(:,:,3)=CS(:,:,3)*Ep.Cs(3);

% real color - highlight effect = diffuse color
Id=Ep.I0-CS;

% real color in chromaticity
IdS=sum(Id,3);
Idr=Id(:,:,1)./IdS;
Idg=Id(:,:,2)./IdS;

% grandient of r,g,I(IdS) in x y 
[del_rx,del_ry]=gradient(Idr);
[del_gx,del_gy]=gradient(Idg);
[del_Ix,del_Iy]=gradient(IdS);

% homo-texture part
E1x=Ep.sx.*(Ep.gamma.*(del_rx.^2+del_gx.^2)+del_Ix.^2);
E1y=Ep.sy.*(Ep.gamma.*(del_ry.^2+del_gy.^2)+del_Iy.^2);
E1=E1x+E1y;

% complex texture part
ICs=sum(CS,3);
[del_ICsx,del_ICsy]=gradient(ICs);
E2x=(1-Ep.sx).*abs(del_ICsx);
E2y=(1-Ep.sy).*abs(del_ICsy);
E2=E2x+E2y;

E_mat=ROI.*(E1+E2);
E_mat(isnan(E_mat))=0;
E=sum(E_mat(:));
