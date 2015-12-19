clear all
%Load the data
load('mesh2');
load('current');
load('sphere');

p=100*p;    %sphere radius is 100 m

k=omega/c_;
K=j*k;

for m=1:EdgesTotal
    Point1=Center(:,TrianglePlus(m));
    Point2=Center(:,TriangleMinus(m));
    DipoleCenter(:,m)=0.5*(Point1+Point2);
    DipoleMoment(:,m)=EdgeLength(m)*I(m)*(-Point1+Point2); 
end

TotalPower=0;
%Sphere series
M=length(t);
for m=1:M
    N=t(1:3,m);
    ObservationPoint=1/3*sum(p(:,N),2);
    [E,H]=point(ObservationPoint,eta_,K,DipoleMoment,DipoleCenter);
    ET=sum(E,2); HT=sum(H,2);
    Poynting(:,m)=0.5*real(cross(ET,conj(HT)));
    U(m)=(norm(ObservationPoint))^2*norm(Poynting(:,m));    
    Vector1=p(:,N(1))-p(:,N(2));
    Vector2=p(:,N(3))-p(:,N(2));
    Area =0.5*norm(cross(Vector1,Vector2)); 
    TotalPower=TotalPower+norm(Poynting(:,m))*Area;
    %------------------------------
    X(1:3,m)=[p(1,N)]';
    Y(1:3,m)=[p(2,N)]';
    Z(1:3,m)=[p(3,N)]';      
end

TotalPower

GainLogarithmic     =10*log10(4*pi*max(U)/TotalPower)
GainLinear          =4*pi*max(U)/TotalPower
RadiationResistance =2*TotalPower/abs(GapCurrent)^2

FileName='gainpower.mat'; 
save(FileName, 'TotalPower','GainLogarithmic','GainLinear');

figure();
U=U/norm(U);
C=repmat(U,3,1);
h=fill3(X,Y,Z,C);
colormap gray;
axis('equal')
rotate3d on

