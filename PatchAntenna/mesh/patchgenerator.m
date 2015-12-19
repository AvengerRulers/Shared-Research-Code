clear;
clc;

Color1=[1.00 1.00 1.00];   
Color2=[0.45 0.45 0.45];    
Color3=[0.70 0.70 0.70];   
Color4=[0.90 0.90 0.90];   

h=0.01;

L=0.10;    
W=0.05;    
Nx=17;     
Ny=9;      

Number=3;

epsilon=1e-6; 
M=1;
for i=1:Nx+1
    for j=1:Ny+1
        X(M)=-L/2+(i-1)/Nx*L;
        Y(M)=-W/2+(j-1)/Ny*W-epsilon*X(M);
        M=M+1;
    end
end

TRI = delaunay(X,Y);
dt = delaunayTriangulation(X', Y');
t=dt.ConnectivityList'; t(4,:)=0;
clear new i a b c
p=[X; Y; -h*ones(1,length(X))];
save plate p t 

PatchNumber=[];
viewer plate; view(0,90); hold on
m=0;
while ~isempty(t)
    m=m+1;
    [xi,yi]=ginput(1); if isempty(xi|yi) break; end;
    TriangleNumber = dt.pointLocation(xi, yi);
    n=t(1:3,TriangleNumber);
    PatchNumber= [PatchNumber TriangleNumber];
    x= p(1,n);
    y= p(2,n);
    fill(x,y,Color4)
    clear xi yi
    clear x y n %
end
t(4,:)=3;
t(4,PatchNumber)=2;

save plate p t

hold off
viewer plate; view(0,90); hold on
FeedingTriangle=[];
TRI=t(1:3,:)';
while ~isempty(t)
    [xi,yi]=ginput(1);
    if isempty(xi|yi) break; end
    TriangleNumber = dt.pointLocation(xi, yi);
    n=t(1:3,TriangleNumber);
    FeedingTriangle= [FeedingTriangle TriangleNumber];
    x= p(1,n);
    y= p(2,n);
    fill(x,y,Color4)
    clear xi yi;
    clear x y n %
end

tbase=t; pbase=p;
tbase(4,:)=0;
p(3,:)=p(3,:)+h;

T=[tbase t+length(pbase)];
T(4,:)=[tbase(4,:) t(4,:)];
P=[pbase p];
p=P; t=T;

FeedingTriangle=[FeedingTriangle FeedingTriangle+length(tbase)];
for n=1:length(FeedingTriangle)/4
    FT=[FeedingTriangle(2*n-1) FeedingTriangle(2*n)];
    N=t(1:3,FT(1));
    M=t(1:3,FT(2));
    a=1-all([N-M(1) N-M(2) N-M(3)]);
    Edge_B=M(find(a)); 
    Edge_T  =[Edge_B'+length(pbase)];
    Edge_MM=Edge_B;
    for k=1:Number-1
        p(:,length(p)+1)=k/Number*(p(:,Edge_T(1))-p(:,Edge_B(1)))+p(:,Edge_B(1));
        p(:,length(p)+1)=k/Number*(p(:,Edge_T(2))-p(:,Edge_B(2)))+p(:,Edge_B(2));
        Edge_M=[length(p)-1,length(p)];
        tFeed1(:,k)  =[Edge_MM(1);Edge_MM(2);Edge_M(2);1];
        tFeed2(:,k)  =[Edge_MM(1);Edge_M(1);Edge_M(2);1];
        Edge_MM=Edge_M;
    end
        
    tFeed3  =[Edge_M(1);Edge_M(2);Edge_T(2);1];
    tFeed4  =[Edge_M(1);Edge_T(1);Edge_T(2);1];
    t=[t tFeed1 tFeed2 tFeed3 tFeed4];
end

new = [];
for i = 1:length(t)
    a = p(:, t(1, i));
    b = p(:, t(2, i));
    c = p(:, t(3, i));
    if (abs(a(2) - b(2)) < 1e-6) && (abs(a(2) - c(2)) < 1e-6)
        continue
    end
    if (abs(a(1) - b(1)) < 1e-6) && (abs(a(1) - c(1)) < 1e-6)
        continue
    end
    new = [new t(:, i)];
end
t = new;

save plate p t h FeedingTriangle
hold off
clear figure
viewer plate
    
    
    




