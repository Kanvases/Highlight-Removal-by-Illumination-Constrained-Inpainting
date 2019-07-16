function alpha = initAlpha(Ep,aROI)
alpha=zeros(size(aROI));
ROI_lab=bwlabel(aROI);

sN=min(Ep.sx,Ep.sy);
sNo=bwmorph(sN,'open',1);
sNlab=bwlabel(sNo);
sNlab(sN==0)=-1;

for roil=1:max(ROI_lab(:))
ROI=zeros(size(aROI));
ROI(ROI_lab==roil)=1;
% when (x,y) in homotexture, sN=1

% find 'processed' colors
hlbound=edge(ROI);
hlbound=imdilate(hlbound,ones(5));
processedPix=zeros(size(ROI));
processedPix(hlbound&sN==1)=1;

Ir=Ep.I0(:,:,1);
pr=Ir(processedPix==1);
Ig=Ep.I0(:,:,2);
pg=Ig(processedPix==1);
Ib=Ep.I0(:,:,3);
pb=Ib(processedPix==1);
[px,py]=find(processedPix==1);
plab=zeros(size(px));
for i=1:length(px)
    plab(i)=sNlab(px(i),py(i));
end
% nDis=(max(px)-min(px))^2+(max(py)-min(py))^2;

% mark of pixels whether found alpha or not
findPix=ones(size(ROI));
findPix(ROI==1&processedPix==0&sN==1)=0;

% sN==1 homo-texture inpaint alpha by colors
while(sum(findPix(:))~=size(ROI,1)*size(ROI,2))
    % list ROI location
    [sx,sy]=find(findPix==0);
    for i=1:length(sx)
        x=sx(i);y=sy(i);
        % calculate alpha that have min distance to processed color in same
        % label
        minDis=10;
        tlab=sNlab(x,y);
        for ta=0:0.05:1
            tr=Ep.I0(x,y,1)-ta*Ep.Cs(1);
            tg=Ep.I0(x,y,2)-ta*Ep.Cs(2);
            tb=Ep.I0(x,y,3)-ta*Ep.Cs(3);
            prl=pr(plab==tlab);
            pgl=pg(plab==tlab);
            pbl=pb(plab==tlab);
            allDis=(prl-tr).^2+(pgl-tg).^2+(pbl-tb).^2;
            dis=min(allDis);
            if(dis<minDis)
                minDis=dis;
                minA=ta;
            end
        end
        
        alpha(x,y)=minA;
        findPix(x,y)=1;
        processedPix(x,y)=1;
    end
end

% sN==0 differnt texture inpaint alpha by smooth
findPix=ones(size(ROI));
findPix(sN==0&ROI==1)=0;
[sx,sy]=find(findPix==0);
for i=1:length(sx)
    x=sx(i);y=sy(i);
    meanAph=0;ia=0;
    for dx=x-5:x+5
        for dy=y-5:y+5
            if(sN(dx,dy)==1)
                meanAph=meanAph+alpha(dx,dy);
                ia=ia+1;
            end
        end
    end
    if(ia>0)
        alpha(x,y)=meanAph/ia;
    end
end
end
end

