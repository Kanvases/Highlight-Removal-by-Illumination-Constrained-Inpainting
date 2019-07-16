function u=Inpainting_TV(f,M,tao,lambda)

    epsilon=1;
    u=zeros(size(f));
    u_old=u;
    while(1)
        gradxu=gradx(u);
        gradyu=grady(u);
        norm_epsilon=sqrt(gradxu.^2+gradyu.^2+epsilon);
        
        increment=lambda*(f-u).*M+div(gradxu./norm_epsilon,gradyu./norm_epsilon);
        u=u+tao*increment;
        
        threshold=norm(u_old(:)-u(:))/(norm(u_old(:))+1e-14);
        if threshold<0.001,break;end
        
        u_old=u;
    
    end
    
end