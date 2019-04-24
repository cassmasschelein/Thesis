%This code provides a way to categorize how many trials converge to
%different behaviours (fully republican, fully democrat, or clustered).
%This code performs the runs with a opinion matrix seeded according to
%North Carolina voter data.

clear
trials=1000;
counter=NaN(2,trials); %fix the check that finds out when it converges


for z=1:trials
    n=5;
    m=20;
    
    D=[.264,.249,.235,.209,.34,.433,.403,.476,.616,.659,.629,.626,.682,.622,.445,.417,.347,.498,.256,.232,.485,.214,.181,.536,.587,.423,.74,.789,.584,.428,.491,.519,.655,.523,.491,.571,.415,.421,.372,.369,.324,.199,.207,.290,.235,.208,.296,.303,.245,.305,.241,.206,.536,.421,.366,.333,.432,.446,.461,.378,.356,.557,.238,.249,.385,.238,.359,.338,.421,.538,.366,.567,.409,.398,.405,.358,.226,.310,.336,.462,.423,.230,.365,.279,.345,.420,.374,.346,.346,.250,.337,.327,.633,.329,.557,.441,.53,.446,.383,.344];
    R=[.711,.727,.742,.765,.64,.548,.576,.501,.370,.326,.360,.364,.305,.371,.535,.560,.628,.476,.714,.730,.470,.766,.796,.434,.387,.552,.23,.185,.379,.546,.493,.463,.333,.450,.496,.419,.569,.561,.594,.624,.649,.783,.772,.684,.741,.768,.676,.671,.726,.672,.734,.773,.435,.552,.607,.640,.549,.544,.524,.596,.614,.411,.742,.726,.585,.740,.618,.634,.552,.430,.607,.407,.576,.590,.584,.624,.710,.656,.640,.503,.552,.747,.595,.694,.625,.539,.599,.626,.628,.729,.642,.643,.334,.640,.431,.542,.452,.514,.604,.631];
    V=zeros(100,1);
    
    for i=1:100
        if R(i)>D(i)
            V(i)=R(i); %format the North Carolina data
        elseif D(i)>R(i)
            V(i)=1-D(i);
        end
    end
    
    A0=zeros(5,20); %format the North Carolina data
    A0(1,1:20)=V(1:20);
    A0(2,1:20)=V(21:40);
    A0(3,1:20)=V(41:60);
    A0(4,1:20)=V(61:80);
    A0(5,1:20)=V(81:100);
    
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
