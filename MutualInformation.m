function [NMI]=MutualInformation(img1,img2)
        [Ma,Na] = size(img1);
        [Mb,Nb] = size(img2);
        m=min(Ma,Mb);
        n=min(Na,Nb); 
        ET=entropy(img1);
        ES=entropy(img2);%//???
        histq=zeros(256,256);%//????????
        %//?????
        for s=1:m
            for t=1:n
                x=img1(s,t)+1;y=img2(s,t)+1;%//??<—>??
                histq(x,y)=histq(x,y)+1;
            end
        end
        p=histq./sum(sum(histq));%//??????
        EST=-sum(sum(p.*log(p+eps)));
        NMI=(ES+ET)-EST;