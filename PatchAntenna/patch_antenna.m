clear;
clc;

load('mesh/plate');
[s1 s2]=size(p);
if(s1==2)
    p(3,:)=0; 
end

Remove=find(t(4,:)>2);
t(:,Remove)=[];           
TrianglesTotal=length(t);

for m=1:TrianglesTotal
   N=t(1:3,m);
   Vec1=p(:,N(1))-p(:,N(2));
   Vec2=p(:,N(3))-p(:,N(2));
   Area(m) =norm(cross(Vec1,Vec2))/2;
   Center(:,m)=1/3*sum(p(:,N),2);
end

Edge_=[];
n=0;
for m=1:TrianglesTotal
    N=t(1:3,m);
    for k=m+1:TrianglesTotal
        M=t(1:3,k);      
        a=1-all([N-M(1) N-M(2) N-M(3)]);
        if(sum(a)==2) 
            n=n+1;
            Edge_=[Edge_ M(find(a))]; 
            TrianglePlus(n)=m;
            TriangleMinus(n)=k; 
        end; 
    end
end
EdgesTotal=length(Edge_);

Edge__=[Edge_(2,:); Edge_(1,:)];
Remove=[];
for m=1:EdgesTotal
    Edge_m=repmat(Edge_(:,m),[1 EdgesTotal]);
    Ind1=any(Edge_  -Edge_m);
    Ind2=any(Edge__ -Edge_m);
    A=find(Ind1.*Ind2==0);
    if(length(A)==3)  
        Out=find(t(4,TrianglePlus(A))==t(4,TriangleMinus(A)));
        Remove=[Remove A(Out)];
    end
end
Edge_(:,Remove)         =[];
TrianglePlus(Remove)    =[];
TriangleMinus(Remove)   =[];
EdgesTotal=length(Edge_)

EdgeIndicator=t(4,TrianglePlus)+t(4,TriangleMinus);

for m=1:EdgesTotal
   EdgeLength(m)=norm(p(:,Edge_(1,m))-p(:,Edge_(2,m)));
end

save mesh1  p ...
            t ...
            Edge_ ...
            TrianglesTotal ...
            EdgesTotal ...
            TrianglePlus ...
            TriangleMinus ...
            EdgeLength ...
            EdgeIndicator ...
            Area ...
            Center ...
            h                
        
clear all
load('mesh1')

IMT=[];
for m=1:TrianglesTotal
    n1=t(1,m);
    n2=t(2,m);
    n3=t(3,m); 
    M=Center(:,m);
    r1=    p(:,n1);
    r2=    p(:,n2);
    r3=    p(:,n3);
    r12=r2-r1;
    r23=r3-r2;
    r13=r3-r1;
    C1=r1+(1/3)*r12;
    C2=r1+(2/3)*r12;
    C3=r2+(1/3)*r23;
    C4=r2+(2/3)*r23;
    C5=r1+(1/3)*r13;
    C6=r1+(2/3)*r13;
    a1=1/3*(C1+C5+r1);
    a2=1/3*(C1+C2+M);
    a3=1/3*(C2+C3+r2);
    a4=1/3*(C2+C3+M);
    a5=1/3*(C3+C4+M);
    a6=1/3*(C1+C5+M);
    a7=1/3*(C5+C6+M);
    a8=1/3*(C4+C6+M);
    a9=1/3*(C4+C6+r3);
    Center_(:,:,m)=...
        [a1 a2 a3 a4 a5 a6 a7 a8 a9];
end

for m=1:EdgesTotal
    NoPlus=TrianglePlus(m);
    n1=t(1,NoPlus);
    n2=t(2,NoPlus);
    n3=t(3,NoPlus); 
    if((n1~=Edge_(1,m))&(n1~=Edge_(2,m))) NODE=n1; end;
    if((n2~=Edge_(1,m))&(n2~=Edge_(2,m))) NODE=n2; end;
    if((n3~=Edge_(1,m))&(n3~=Edge_(2,m))) NODE=n3; end;
    FreeVertex=p(:,NODE);
    
    RHO_Plus(:,m)   =+Center(:,NoPlus)-FreeVertex;
    RHO__Plus(:,:,m)  =...
        +Center_(:,:,NoPlus)-repmat(FreeVertex,[1 9]);
end

for m=1:EdgesTotal
    NoMinus=TriangleMinus(m);
    n1=t(1,NoMinus);
    n2=t(2,NoMinus);
    n3=t(3,NoMinus); 
    if((n1~=Edge_(1,m))&(n1~=Edge_(2,m))) NODE=n1; end;
    if((n2~=Edge_(1,m))&(n2~=Edge_(2,m))) NODE=n2; end;
    if((n3~=Edge_(1,m))&(n3~=Edge_(2,m))) NODE=n3; end;
    FreeVertex=p(:,NODE);
    
    RHO_Minus(:,m)   =-Center(:,NoMinus) +FreeVertex;
    RHO__Minus(:,:,m)=...
        -Center_(:,:,NoMinus)+repmat(FreeVertex,[1 9]);
end

save mesh2  p ...
            t ...            
            TrianglesTotal ...
            EdgesTotal ...
            Edge_ ...
            TrianglePlus ...
            TriangleMinus ...
            EdgeLength ...
            EdgeIndicator ...
            Area ...
            RHO_Plus ...
            RHO_Minus ...
            RHO__Plus ...
            RHO__Minus ...
            Center ...
            Center_ ...
            h           
        
clear all
load('mesh2');

NumberOfSteps=21;
FreqStart   =1.0e9;     
FreqStop    =6.0e9;    
step=(FreqStop-FreqStart)/(NumberOfSteps-1);

epsilon_    =8.854e-012;
epsilon_R   =1.0;
mu_         =1.257e-006;
c_=1/sqrt(epsilon_*mu_);
eta_=sqrt(mu_/epsilon_);

for m=1:EdgesTotal
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);  
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  
end

DP      =find(t(4,:)==0);
M       =length(DP);       
delta   =[0;0;1e-12];     

Middle=[0; 0; h/2];
for m=1:M
    N=t(1:3,m);
    Point(:,m)=Center(:,m)      +Middle; 
    IMT(:,:,m)=Center_(:,:,m)   +repmat(Middle,[1,9]);
end

for FF=1:NumberOfSteps    
    FF
    f(FF)       =FreqStart+step*(FF-1);
    omega       =2*pi*f(FF);
    k           =omega/c_;
    K           =j*k;
  
    Constant1   =mu_/(4*pi);
    Constant2   =1/(j*4*pi*omega*epsilon_);
    Factor      =1/9;    
    FactorA     =Factor*(j*omega*EdgeLength/4)*Constant1;
    FactorFi    =Factor*EdgeLength*Constant2;
    FactorA     =FactorA.';
    FactorFi    =FactorFi.';
    
    ZSS=  impmet( EdgesTotal,TrianglesTotal,...
            EdgeLength,K,...
            Center,Center_,...
            TrianglePlus,TriangleMinus,...
            RHO_P,RHO_M,...
            RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);   
    ZSS=ZSS.'; 

    ZDD = zeros(M,M)+j*zeros(M,M);
    for m=1:M
        OP      =Point(:,m)+delta;
        IP      =Point;
        E       =point_(OP,K,k,Constant2,Area(1:M),h,IP);
        ZDD(m,:)=E(3,:);        
        IP      =IMT(:,:,m);
        E       =point_(OP,K,k,Constant2,ones(1,9)*Area(m)/9,h,IP);
        ZDD(m,m)=sum(E(3,:));   
    end    
    ZDD=j*omega*epsilon_*(epsilon_R-1)*ZDD; 

    ZDS = zeros(EdgesTotal,M)+j*zeros(EdgesTotal,M);
    for m=1:EdgesTotal
        OPPlus  =Center(:,TrianglePlus(m))+delta;
        OPMinus =Center(:,TriangleMinus(m))+delta;
        EP      =point_(OPPlus,K,k,Constant2,Area(1:M),h,Point);
        EM      =point_(OPMinus,K,k,Constant2,Area(1:M),h,Point);
        ScalarPlus  =sum(EP.* repmat(RHO_Plus(:,m), [1 M]));
        ScalarMinus =sum(EM.* repmat(RHO_Minus(:,m),[1 M]));
        ZDS(m,:)    =EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2);    
    end

    ZSD = zeros(M,EdgesTotal)+j*zeros(M,EdgesTotal);
    C1=1/(2*epsilon_)/(-j*omega);
    C2=j*omega*epsilon_*(epsilon_R-1);
    C=C1*C2;
    for m=1:M         
        Q=sum(abs([Center(1,m)-Center(1,:); Center(2,m)-Center(2,:)]));
        T=find(Q<1.e-9);     
        for q=1:length(T)
            Plus    =find(TrianglePlus-T(q)==0);
            Minus   =find(TriangleMinus-T(q)==0);
            Ind     =(-1)^(q+1);
            for k=1:length(Plus)
                n=Plus(k);
                Charge=EdgeLength(n)/Area(T(q));
                ZSD(m,n)=ZSD(m,n)+C*Ind*Charge;
            end
            for k=1:length(Minus)
                n=Minus(k);
                Charge=-EdgeLength(n)/Area(T(q));
                ZSD(m,n)=ZSD(m,n)+C*Ind*Charge; 
            end 
        end
    end
  
    for m=1:EdgesTotal    
        V(m)=0;
    end

    Index=find(EdgeIndicator==1);

    V(Index)=1*EdgeLength(Index);   
    
    Z=[ZSS -ZDS; ZSD ZDD-eye(M,M)];

    V(EdgesTotal+1:EdgesTotal+M)=0;
    
    I=Z\V.';
    
    CURRENT(:,FF)=I(:);
    GapCurrent(FF)  =sum(I(Index).*EdgeLength(Index)');
    GapVoltage(FF)  =mean(V(Index)./EdgeLength(Index));
    Impedance(FF)   =GapVoltage(FF)/GapCurrent(FF);
    FeedPower(FF)   =1/2*real(GapCurrent(FF)*conj(GapVoltage(FF)));    
    Imp             =Impedance(FF)
end

FileName='current.mat'; 
save(FileName, 'f','NumberOfSteps','FreqStart','FreqStop','step',...
                'omega','mu_','epsilon_','c_', 'eta_',...
                'CURRENT','GapCurrent','GapVoltage','Impedance','FeedPower','M','h','Point','Index');   

clear all
load('current.mat');

a=figure
plot(f, real(Impedance),'.',f,imag(Impedance),'.');
hold on
plot(f, real(Impedance),f,imag(Impedance),'--');
xlabel ('Frequency, Hz')
ylabel('Input  resistance/reactance, Ohm')
title('Resistance-solid; reactance-dashed')
grid on

b=figure
Gamma=(Impedance-50)./(Impedance+50);
Out=20*log10(abs(Gamma));
plot(f, Out);
xlabel ('Frequency, Hz')
ylabel ('Return loss, dB')
grid on
hold on


 