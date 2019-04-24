%This code generates the adjacancy matrix that can then be used in the R code
%to find shortest path length. To use this code, plug in the parameters
%wanted and then save the object called 'friends' in the .mat format. This
%will be used in the R code titled 'ShortestPathLength.R' on this same
%respository.

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
    %zeta=1;
    epsilon=(1-0.75)*rand(1,1)+0.75; %how willing to make friends
    %epsilon=1;
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


%% this section determines the distance between nodes (note that this distance is a quasidistance)

distance=zeros(N,N);

for i=1:N
    for j=1:N
        distance(i,j)= acosh(cosh(nodes(3,i,1)*nodes(1,i,1))*cosh(nodes(3,i,1)*nodes(1,j,1)) - sinh(nodes(3,i,1)*nodes(1,i,1))*sinh(nodes(3,i,1)*nodes(1,j,1))*cos(pi-abs(pi-abs(nodes(2,i,1)-nodes(2,j,1)))))/nodes(3,i,1);
    end
end

for i=1:N
    distance(i,i)=0;
end

%% this section determines which nodes are indeed friends -i.e. the edges of the graph

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
