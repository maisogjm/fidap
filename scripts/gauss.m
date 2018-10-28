function G=gauss(V)
% 08.01.95 : JM Maisog
%
% Function returns a vector containing a Gaussian with the specified
% mean (mu), standard deviation (sd), and length (N)
% 
mu=V(1);
sd=V(2);
N=V(3);
G=zeros(1,N);
if (sd~=0)
    factor1=(1/(sd*sqrt((2.*pi))));
else
    factor1=1;
end
for i=0:N-1
    kernel_index=abs(i-mu);
    if (i>(mu+(N-1)/2))
        kernel_index=N-i+mu;
    end
    if (sd~=0)
        factor2=exp(-(kernel_index*kernel_index)/(2.*sd*sd));
    else
        if (((mu>=0.) & (abs(mu-i)<0.5)) | ((mu<0.) & abs(N-i+mu)<=0.5) ...
	| ((mu<0.) & abs(mu-i)<0.5))
            % G(i+1)=1.;
            factor2=1;
        else
            % G(i+1)=0.;
            factor2=0;
        end
    end

    G(i+1)=factor1*factor2;
end
G=G/sum(G);
