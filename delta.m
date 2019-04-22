%This code calculates the delta value for the four point condition for a
%delta hyperbolic space. It takes the largest delta over 10,000 trials. The
%output of this code is a delta value which prints in the command window.
clear
close all
N=2500; %agents in the network
alpha=2.5;

sizeOut = [1, 1]; % sample size
r1 = 0;  % lower bound
r2 = log(N); % upper bound


for i=1:N
    r=exprndBounded(alpha, sizeOut, r1, r2); %give radial value an exp distribution
    theta=2*pi*rand(1,1);
   zeta=(log(N)/3-1)*rand(1,1)+1; %individual curvature - effects how distance is measured
    %zeta=1;
    epsilon=(1.75-0.75)*rand(1,1)+0.75; %how willing to make friends
    %epsilon=1;
    nodes(1,i,1)=r;
    nodes(2,i,1)=theta;
    nodes(3,i,1)=zeta;
    nodes(4,i,1)=epsilon;
end


distance=zeros(N,N);

for i=1:N
    for j=1:N
        distance(i,j)= acosh(cosh(nodes(3,i,1)*nodes(1,i,1))*cosh(nodes(3,i,1)*nodes(1,j,1)) - sinh(nodes(3,i,1)*nodes(1,i,1))*sinh(nodes(3,i,1)*nodes(1,j,1))*cos(pi-abs(pi-abs(nodes(2,i,1)-nodes(2,j,1)))))/nodes(3,i,1);
    end
end

for i=1:N
    distance(i,i)=0;
end
trials=10;
w=zeros(1,trials);
x=zeros(1,trials);
y=zeros(1,trials);
z=zeros(1,trials);
delta=0;
for k=1:trials
    w(k)=randi([1,N],1,1);
    x(k)=randi([1,N],1,1);
    y(k)=randi([1,N],1,1);
    z(k)=randi([1,N],1,1);
    newdelta=0.5*[distance(x(k),z(k))+distance(y(k),w(k))-max(distance(x(k),w(k))+distance(y(k),z(k)),distance(x(k),y(k))+distance(w(k),z(k)))];
    if newdelta>delta
        delta=newdelta;
    end
end

delta
%%
function r = exprndBounded(alpha, sizeOut, r1, r2)

minE = exp(r1*alpha);
maxE = exp(r2*alpha);

randBounded = minE + (maxE-minE).*rand(sizeOut);
r = (1/alpha).* log(randBounded);
end
