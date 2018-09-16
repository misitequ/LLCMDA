clear
load mfs;
load dss;
load knownre;
load pmi;
lodd pdi;
nd=383;
nm=495;
pp=5430;
interaction=knownre;
save pmi pmi;
save pdi pmi;
I=ones(495,1);
P=knownre;
Q=knownre;
for i=1:nm
    interaction=P;
    xm=interaction(:,i);
    Cm=(interaction-xm*I')'*(interaction-xm*I');
    w1=pinv((Cm+(diag(pmi(i,:)).^2)))*I;
    w=w1/(I'*w1);
    W(:,i)=w;
end
W(W<0)=0;
W=(W+W')/2;
I=ones(383,1);
for i=1:nd
    interaction=Q';
    xd=interaction(:,i);
    Cd=(interaction-xd*I')'*(interaction-xd*I');
    w2=pinv(Cd+(diag(pdi(i,:)).^2))*I;
    w=w2/(I'*w2);
    W1(:,i)=w;
end
W1(W1<0)=0;
W1=(W1+W1')/2;
M2=sum(W);
for i=1:495
    for j=1:495
        W(i,j)=W(i,j)/(((M2(i)*M2(j))^0.5));
    end
end
D2=sum(W1);
for i=1:383
    for j=1:383
        W1(i,j)=W1(i,j)/(((D2(i)*D2(j))^0.5));
    end
end
M1=sum(mfs);
for i=1:495
    for j=1:495
        mfs(i,j)=mfs(i,j)/(((M1(i)*M1(j))^0.5));
    end
end
D1=sum(dss);
for i=1:383
    for j=1:383
        dss(i,j)=dss(i,j)/(((D1(i)*D1(j))^0.5));
    end
end
gamma=0.5;
PM=knownre';
PM3=knownre';
PD=knownre;
PD3=knownre;
P0=knownre';
delta = 1;
while  (delta > 1e-6)
    PM1 = (1-gamma)*mfs*PM+gamma*P0;
    delta =abs(sum(sum((abs(PM1)-abs(PM)))));
    PM = PM1;
end
delta = 1;
while  (delta > 1e-6)
    PM2 = (1-gamma)*W*PM3+gamma*P0;
    delta =abs(sum(sum((abs(PM2)-abs(PM3)))));
    PM3 = PM2;
end
delta = 1;
while  (delta > 1e-6)
    PD1 = (1-gamma)*dss*PD+gamma*P0';
    delta =abs(sum(sum((abs(PD1)-abs(PD)))));
    PD =PD1;
end
delta = 1;
while  (delta > 1e-6)
    PD2 = (1-gamma)*W1*PD3+gamma*P0';
    delta =abs(sum(sum((abs(PD2)-abs(PD3)))));
    PD3 =PD2;
end
F=(PM1+PM2+PD1'+PD2');
F=F';
finalprediction=[];
for i=1:nm
    finalprediction=[finalprediction;F(:,i)];
end
finalpredictiondisease=[];
for i=1:nm
    finalpredictiondisease=[finalpredictiondisease,1:nd];
end
finalprediction(:,2)=finalpredictiondisease';
finalpredictionmicrobe=[];
for i=1:nm
    finalpredictionmicrobe=[finalpredictionmicrobe;i.*ones(nd,1)];
end
finalprediction(:,3)=finalpredictionmicrobe;
discard=finalprediction(find(finalprediction(:,1)==-10000),:);
finalprediction=setdiff(finalprediction,discard,'rows');
result=sortrows(finalprediction,-1);
[a1,~,a]=xlsread('disease');

[b1,~,b]=xlsread('miRNA');

m=length(result);
c=zeros(m,1);
f=zeros(m,1);
d={};
e={};
for i=1:m
    if(~isempty(find(a1==result(i,2))))
        c(i,1)=find(a1==result(i,2));
        d{i,1}=a{c(i,1),2};
    else
        d{i,1}=result(i,2);
    end
end
for i=1:m
    if(~isempty(find(b1==result(i,3))))
        f(i,1)=find(b1==result(i,3));
        e{i,1}=b{f(i,1),2};
    else
        e{i,1}=result(i,3);
    end
end
save result result