%% This code implements opinion dynamics on complex networks. It distributes the nodes and builds the edges of the social network.
%It assigns opinions to the nodes and allows for them to co-evolve.
%Throughout the code we also provide visualization of the network in
%various ways, these sections are labelled. It also checks the clusetring
%coefficient and the power law distribution of nodes
%% Assigns to each node in the network a r, theta, epsilon, zeta, and opinion
clear
close all
t_end=3000; %number of updates
N=5000; %agents in the network
nodes=zeros(5,N,t_end); %store information regarding each individual
alpha=2.5;

sizeOut = [1, 1]; % sample size
r1 = 0;  % lower bound
r2 = log(N); % upper bound


for i=1:N
    r=exprndBounded(alpha, sizeOut, r1, r2); %give radial value an exp distribution
    theta=2*pi*rand(1,1);
    zeta=(log(N)/3-1)*rand(1,1)+1; %individual curvature - effects how distance is measured
    %zeta=1.5;
    epsilon=(1-0.75)*rand(1,1)+0.75; %how willing to make friends
    nodes(1,i,1)=r;
    nodes(2,i,1)=theta;
    nodes(3,i,1)=zeta;
    nodes(4,i,1)=epsilon;
    nodes(5,i,1)=rand(1,1); %continuous opinion of an individual
end

for j=1:t_end
    nodes(1,:,j)=nodes(1,:,1); %remains constant
    nodes(2,:,j)=nodes(2,:,1); %remains constant
    nodes(3,:,j)=nodes(3,:,1); %remains constant
    nodes(4,:,j)=nodes(4,:,1); %remains constant
end


%% this section determines the distance between nodes (directed)
distance=zeros(N,N);

for i=1:N
    for j=1:N
        distance(i,j)= acosh(cosh(nodes(3,i,1)*nodes(1,i,1))*cosh(nodes(3,i,1)*nodes(1,j,1)) - sinh(nodes(3,i,1)*nodes(1,i,1))*sinh(nodes(3,i,1)*nodes(1,j,1))*cos(pi-abs(pi-abs(nodes(2,i,1)-nodes(2,j,1)))))/nodes(3,i,1);
    end
end

for i=1:N
    distance(i,i)=0;
end

%% this section determines which nodes are friends -i.e. the edges of the graph
request=distance;

for i=1:N
    for j=1:N
        if distance(i,j)<=(log(N)*nodes(4,i,1)) %send a friend request if person meets distance requirement
            request(i,j)=1;
        else
            request(i,j)=0;
        end
    end
end

friends=request;

for i=1:N
    for j=1:N
        if friends(i,j)==1 && friends(j,i)==0 %accept friend request if both people agree on distance threshold
            friends(i,j)=0;
        end
    end
end

for i=1:N
    friends(i,i)=0;
end
%A=nodes(1,:);
%plot(sort(A))

%% this plots the network and visualizes the connections of one person
[err i] = min(abs(sum(friends')-mean(sum(friends')))) %plot the friends of the person closest to the mean number of friends
x=NaN(1,N); %assume connection does not exist
for j=1:N
    if friends(i,j)==1
        x(j)=1; %add in a connection if they are friends/overwrite previous definition
    end
    thing1=[nodes(2,i,1) nodes(2,j,1)*x(j)]; %the two friends theta values
    thing2=[nodes(1,i,1) nodes(1,j,1)*x(j)]; %the two friends radial values
    figure(1)
    polarplot(nodes(2,:,1),nodes(1,:,1),'.','MarkerSize',15)
    polarplot(thing1,thing2,'-','Color','k')
    hold on;
end

%% this is a visualization of the network where individuals are represented by a dot with size corresponding to their popularity and colour representing their opinion
sz=zeros(1,N);
for i=1:N
    sz(i)=2*sum(friends(i,:))+20; %give people a size based on popularity/number of friends
end
c=nodes(5,:,1);
[err i] = min(abs(sum(friends')-mean(sum(friends')))) %plot the friends of the person closest to the mean number of friends
x=NaN(1,N);
colormap(redblue)
figure(2)
colormap(redblue)
polarscatter(nodes(2,:,1),nodes(1,:,1),sz,c,'filled','MarkerEdgeColor','k')

%% Monte Carlo Simulation for opinion update
for t=1:t_end
    p=randi([1 N],1,1);
    opinion=zeros(1,N);
    for i=1:N
        opinion(i)=nodes(5,i,t);
    end
    weight=friends(p,:)./distance(p,:);%make friendships weighted by taking inverse of distance measure
    weight(p)=0;
    new_opinion_almost=0;
    for j=1:N
        new_opinion_almost=weight(j)*log(opinion(j)/(1-opinion(j)))+new_opinion_almost;
    end
    new_opinion=log(opinion(p)/(1-opinion(p)))+new_opinion_almost;
    new_opinion=exp(new_opinion);
    nodes(5,:,t+1)=nodes(5,:,t);
    new_opinion=(new_opinion)/(new_opinion+1);
    nodes(5,p,t+1)=new_opinion;
    for g=1:N
        if nodes(5,g,t+1)>=1
            nodes(5,g,t+1)=0.9999;
        end
        if nodes(5,g,t+1)<=0
            nodes(5,g,t+1)=0.0001;
        end
    end
    
    sz=zeros(1,N);
    for i=1:N
        sz(i)=2*sum(friends(i,:))+20;
    end
    c=nodes(5,:,t);
    
    %         lps=num2str(t);
    %         figure(3)
    %         colormap(redblue)
    %         polarscatter(nodes(2,:,t),nodes(1,:,t),sz,c,'filled','MarkerEdgeColor','k')
    %         title(['Continuous Opinion ',lps]);
    %         M(t)=getframe;
end

%% This plots the percentage of people with democrat or republican opinion over time
t=1:1:t_end;
R=zeros(1,t_end); D=zeros(1,t_end);
R(1)=0; D(1)=0;
for j=1:t_end
    R(j)=(sum(nodes(5,:,j)>0.5)/N)*100;
    D(j)=(sum(nodes(5,:,j)<=0.5)/N)*100;
end

hold on;
figure(4)
plot(t,R,'r','LineWidth',3)
hold on;
plot(t,D,'b','LineWidth',3)
legend('% Republicans','% Democrats')
title('Opinion Dynamics')
xlabel('Number of interaction cycles')
ylabel('Percentage of the population in favour')
hold off;

%% Probability distribution check
k=zeros(1,N);
for i=1:N
    k(i)=sum(friends(i,:));
end

figure(5);
histogram(k)
xlabel('Node degree k')
ylabel('Number of agents with degree k')
title('Histogram for degree distribution')

mean(k)
std(k)

number=1:1:max(k);
y=zeros(1,max(k));
for i=1:max(k)
    y(i)=sum(k(:) == i-1);
end
% figure(6)
% plot(number,y,'o')
%% histogram with no zeros
begin=find(y==max(y))
z=find(y==0)
z1=find(z>begin)
thing=z1(1)
first=z(thing)-1;
y1=y(begin:first);
y1=y1/sum(y1)
p = polyfit(log(number(begin:first)),log(y1),1); 
a = p(1);
% Accounting for the log transformation
X=number(begin:first)
k = exp(p(2));
figure(7)
ezplot(@(X) k*X.^a,[X(1) X(end)])
title([])
xlabel('k')
ylabel('P(k)')
hold on;
plot(X,y1,'o')
legend('~exp(-2.6321)','Histogram Data')
%% fitting a  power law to the histogram
y1=y(begin:end)
y1=y1+0.001
p = polyfit(log(number(begin:end)),log(y1),1); 
a = p(1);
% Accounting for the log transformation
X=number(begin:end)
k = exp(p(2));
figure(7)
ezplot(@(X) k*X.^a,[X(1) X(end)])
hold on;
plot(X,y1,'o')


%% Clustering coefficient.

%calculate the number of triangles around a node i
triangles=zeros(1,N);
for i=1:N
    for j=1:N
        for h=1:N
            triangles(i)=friends(i,j)*friends(i,h)*friends(j,h)+triangles(i);
        end
    end
    i
end
triangles=0.5*triangles;

C=0;
for i=1:N
    if k(i)>2
    C=2*triangles(i)/(k(i)*(k(i)-1))+C;
    else
    C=C;
    end
end
C=C/N;

%% this is the function that exponentially distributes nodes within a defined radius
function r = exprndBounded(alpha, sizeOut, r1, r2)

minE = exp(r1*alpha);
maxE = exp(r2*alpha);

randBounded = minE + (maxE-minE).*rand(sizeOut);
r = (1/alpha).* log(randBounded);
end
