function N=fa_vector(p1,p2,p3)

a=p2-p1;
b=p3-p1;

c=cross(b,a);

norm = sqrt(c(:,:,1).^2+c(:,:,2).^2+c(:,:,3).^2);
N(:,:,1)=c(:,:,1)./norm;
N(:,:,2)=c(:,:,2)./norm;
N(:,:,3)=c(:,:,3)./norm;

end
