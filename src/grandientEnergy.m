function G= grandientEnergy(aph,Ep,ROI)
alpha=zeros(size(ROI));
alpha(ROI==1)=aph;

G_mat=zeros(size(ROI));

% useful in granient
sumC=sum(Ep.Cs);
sumI0=sum(Ep.I0,3);

% specular color * alpha = highlight effect in image(CS) 
CS=repmat(alpha,1,1,3);
CS(:,:,1)=CS(:,:,1)*Ep.Cs(1);
CS(:,:,2)=CS(:,:,2)*Ep.Cs(2);
CS(:,:,3)=CS(:,:,3)*Ep.Cs(3);

% real color - highlight effect = diffuse color
Id=Ep.I0-CS;

% real color in chromaticity
IdS=sum(Id,3);
Idr=Id(:,:,1)./(IdS+0.001);
Idg=Id(:,:,2)./(IdS+0.001);

for x=1:size(ROI,1)
    for y=1:size(ROI,2)
        if(ROI(x,y)==1)
        
            % --- (x,y-1) ---%
            del_rx_2_g=-0.5*(Idr(x,y)-Idr(x,y-2))*(Ep.Cs(1)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,1)-Ep.Cs(1)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_gx_2_g=-0.5*(Idg(x,y)-Idg(x,y-2))*(Ep.Cs(2)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,2)-Ep.Cs(2)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_Ix_2_g=-0.5*(IdS(x,y)-IdS(x,y-2))*(sumC);
            E1x_g=Ep.sx(x,y-1)*(Ep.gamma*(del_rx_2_g+del_gx_2_g)+del_Ix_2_g);
            if alpha(x,y)>alpha(x,y-2)
                E2x_g=0.5*(1-Ep.sx(x,y-1))*sumC;
            else
                E2x_g=-0.5*(1-Ep.sx(x,y-1))*sumC;
            end
            G_l=E1x_g+E2x_g;
            % --- in (x,y+1) ---%
            del_rx_2_g=0.5*(Idr(x,y+2)-Idr(x,y))*(Ep.Cs(1)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,1)-Ep.Cs(1)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_gx_2_g=0.5*(Idg(x,y+2)-Idg(x,y))*(Ep.Cs(2)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,2)-Ep.Cs(2)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_Ix_2_g=0.5*(IdS(x,y+2)-IdS(x,y))*(sumC);
            E1x_g=Ep.sx(x,y+1)*(Ep.gamma*(del_rx_2_g+del_gx_2_g)+del_Ix_2_g);
            if alpha(x,y+2)>alpha(x,y)
                E2x_g=-0.5*(1-Ep.sx(x,y+1))*sumC;
            else
                E2x_g=0.5*(1-Ep.sx(x,y+1))*sumC;
            end
            G_r=E1x_g+E2x_g;
            % --- in (x-1,y) ---%
            del_ry_2_g=-0.5*(Idr(x,y)-Idr(x-2,y))*(Ep.Cs(1)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,1)-Ep.Cs(1)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_gy_2_g=-0.5*(Idg(x,y)-Idg(x-2,y))*(Ep.Cs(2)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,2)-Ep.Cs(2)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_Iy_2_g=-0.5*(IdS(x,y)-IdS(x-2,y))*(sumC);
            E1y_g=Ep.sy(x-1,y)*(Ep.gamma*(del_ry_2_g+del_gy_2_g)+del_Iy_2_g);
            if alpha(x,y)>alpha(x-2,y)
                E2y_g=0.5*(1-Ep.sy(x-1,y))*sumC;
            else
                E2y_g=-0.5*(1-Ep.sy(x-1,y))*sumC;
            end
            G_u=E1y_g+E2y_g;
            % --- (x+1,y) ---%
            del_ry_2_g=0.5*(Idr(x+2,y)-Idr(x,y))*(Ep.Cs(1)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,1)-Ep.Cs(1)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_gy_2_g=0.5*(Idg(x+2,y)-Idg(x,y))*(Ep.Cs(2)*(sumI0(x,y)-sumC*alpha(x,y))-sumC*(Ep.I0(x,y,2)-Ep.Cs(2)*alpha(x,y)))/(sumI0(x,y)-sumC*alpha(x,y)+0.01)^2;
            del_Iy_2_g=0.5*(IdS(x+2,y)-IdS(x,y))*(sumC);
            E1y_g=Ep.sy(x+1,y)*(Ep.gamma*(del_ry_2_g+del_gy_2_g)+del_Iy_2_g);
            if alpha(x+2,y)>alpha(x,y)
                E2y_g=-0.5*(1-Ep.sy(x+1,y))*sumC;
            else
                E2y_g=0.5*(1-Ep.sy(x+1,y))*sumC;
            end
            G_d=E1y_g+E2y_g;
            G_mat(x,y)=G_l+G_r+G_u+G_d;
        end
    end
end

G=G_mat(ROI==1);

