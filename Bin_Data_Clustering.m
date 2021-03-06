%This code provides a way to categorize how many trials converge to
%different behaviours (fully republican, fully democrat, or clustered).
%This code performs the runs with a opinion matrix with clustering

clear
trials=1000;
counter=NaN(2,trials); %fix the check that finds out when it converges


for z=1:trials
    n=5;
    m=20;
    
    for i=1:5
        for j=1:20
            if i<4 && i>1 && j<6 && j>1
                A(i,j,1)=0.5+(rand/2); %define where we want our clusters to be
            elseif i<3 && j<14 && j>9
                A(i,j,1)=0.5+(rand/2);
            elseif i>3 && j<19 && j>13
                A(i,j,1)=0.5+(rand/2);
            elseif i<3 && j<20 &&j>16
                A(i,j,1)=0.5+(rand/2);
            elseif i<4 && j<18 && j>13
                A(i,j,1)=0.5-(rand/2);
            elseif i<3 && j>5 && j<10
                A(i,j,1)=0.5-(rand/2);
            elseif i>3 && j<13 && j>6
                A(i,j,1)=0.5-(rand/2);
            elseif i>3 && j<4
                A(i,j,1)=0.5-(rand/2);
            else
                A(i,j,1)=rand;
            end
        end
    end
    
    alpha=0.5+(1-0.5).*rand(1,1);
    beta=0.5+(1-0.5).*rand(1,1);
    a=log(alpha/(1-beta));
    b=log(beta/(1-alpha));
    tend=1000000;
    
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
    
    % Figure to show popular opinion dynamics
    
    t=1:1:tend;
    R=zeros(1,tend); D=zeros(1,tend);
    R(1)=0; D(1)=0;
    for j=1:tend
        R(j)=sum(sum(A(:,:,j)>0.5));
        D(j)=sum(sum(A(:,:,j)<=0.5));
    end
    
    
    % hold on;
    % figure(1)
    % plot(R,'r', 'LineWidth', 3)
    % plot(D,'b','LineWidth', 3)
    % legend('% Republicans','% Democrats')
    % title('Dynamics of Voter Opinions')
    % xlabel('Number of neighbour interaction cycles')
    % ylabel('Percentage of the population in favour')
    % hold off;
    
    for t=1:(tend-10000)
        if R(t)==100 && R(t+50)==100 && R(t+100)==100 && R(t+200)==100 && R(t+1000)==100 && R(t+10000)==100
            counter(1,z)=1;
            counter(2,z)=t;
        end
        if R(t)==0 && R(t+50)==0 && R(t+100)==0 && R(t+200)==0 && R(t+1000)==0 && R(t+10000)==0
            counter(1,z)=0;
            counter(2,z)=t;
        end
    end
    z
end
dem=sum(counter(1,:)==0);
repub=sum(counter(1,:)==1);
cluster=trials-(dem+repub);
figure(2)
bar([dem,repub,cluster])
name={'Democrat', 'Republican', 'Clustering'}
set(gca,'xticklabel',name)
