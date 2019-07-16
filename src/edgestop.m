function s = edgestop(N,t)
s=ones(size(N));
s(N>t)=0;
end

