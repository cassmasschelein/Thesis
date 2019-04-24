%This code provides a way to see opinions changing in time. It displays the
%continuous and discrete opinions of individuals for a randomly seeded initial condition
%opinion matrix over a specified time frame. NOTE: To be able to run this
%code, you must add the code 'redblue.m' to the path. This code is
%available from the same directory.

clear
n=10; %rows
m=10; %columns

for i=1:n
    for j=1:m
        A(i,j,1)=rand;
    end
end

alpha=0.5+(1-0.5).*rand(1,1);
beta=0.5+(1-0.5).*rand(1,1);
a=log(alpha/(1-beta));
b=log(beta/(1-alpha));
tend=10000;

B(2:(n+1),2:(m+1),1)=A(1:n,1:m,1);
B(1,2:(m+1),1)=A(2,1:m,1);
B(n+2,2:(m+1),1)=A(n-1,1:m,1);
B(2:(n+1),1,1)=A(1:n,2,1);
B(2:(n+1),m+2,1)=A(1:n,m-1,1);
B(1,1,1)=A(2,2,1);
B(1,m+2,1)=A(2,m-1,1);
B(n+2,1,1)=A(n-1,2,1);
B(n+2,m+2,1)=A(n-1,m-1,1);

for i=1:(n+2)
    for j=1:(m+2)
        if B(i,j,1)>0.5 %Republican
            C(i,j,1)=a;
        else
            C(i,j,1)=(-b);
        end
    end
end

%update loop

for t=1:tend
    i=randi([2 (n+1)],1,1);
    j=randi([2 (m+1)],1,1);
    N=C(i-1,j,t);
    S=C(i+1,j,t);
    E=C(i,j+1,t);
    W=C(i,j-1,t);
    NE=C(i-1,j+1,t);
    NW=C(i-1,j-1,t);
    SE=C(i+1,j+1,t);
    SW=C(i+1,j-1,t);
    V=log(A(i-1,j-1,t)/(1-A(i-1,j-1,t))); %log-odds vi are defined on the interval -ve infinity to +ve infinity and can be regarded as a continuous field over site i, with the choice of agent i being defined by a spin variable si(vi)-sign(vi)
    VP=V+0.01*N+0.01*S+0.01*E+0.01*W+0.0025*NW+0.0025*NE+0.0025*SW+0.0025*SE; %a bayesian update of individual i's belief/probability upon observation of its social neighbour's choices
    L=exp(VP);
    A(:,:,t+1)=A(:,:,t);
    A(i-1,j-1,t+1)=(L/(L+1));
    
    B(2:(n+1),2:(m+1),t+1)=A(1:n,1:m,t+1);
    B(1,2:(m+1),t+1)=A(2,1:m,t+1);
    B(n+2,2:(m+1),t+1)=A(n-1,1:m,t+1);
    B(2:(n+1),1,t+1)=A(1:n,2,t+1);
    B(2:(n+1),m+2,t+1)=A(1:n,m-1,t+1);
    B(1,1,t+1)=A(2,2,t+1);
    B(1,m+2,t+1)=A(2,m-1,t+1);
    B(n+2,1,t+1)=A(n-1,2,t+1);
    B(n+2,m+2,t+1)=A(n-1,m-1,t+1);
    
    for i=1:(n+2)
        for j=1:(m+2)
            if B(i,j,t+1)>0.5
                C(i,j,t+1)=a;
            else
                C(i,j,t+1)=-b;
            end
        end
    end
end

t=1:1:tend;
R=zeros(1,tend); D=zeros(1,tend);
R(1)=0; D(1)=0;
for j=1:tend
    R(j)=sum(sum(A(:,:,j)>0.5));
    D(j)=sum(sum(A(:,:,j)<=0.5));
end

%% Movie time
%in the left movie you see the continuous opinions of individuals changing
%as update interactions occur. On the right movie you see the discrete
%opinions of the neighbourhood - the way that a randomly selected
%individual is defined to view their neighbours

for t=1:tend
    lps=num2str(t);
    figure(2)
    subplot(1,2,1)
    imagesc(A(:,:,t));
    colormap(redblue);
    title(['Continuous Opinion ',lps]);
    subplot(1,2,2)
    imagesc(C(2:n+1,2:m+1,t));
    colormap(redblue);
    title(['Discrete Analysis ',lps]);
    M(t)=getframe;
end
