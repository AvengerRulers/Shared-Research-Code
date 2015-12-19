clear;
clc;

f = 6e7;
FeedPoint=[0; 0; 0];

load('mesh/dipole');
[s1 s2] = size(p);
if (s1 == 2)
    p(3,:) = 0;
end

Remove = find(t(4,:) > 1);   
t(:,Remove) = [];           
TrianglesTotal = length(t);

for m = 1:TrianglesTotal
   N = t(1:3,m);
   Vec1 = p(:,N(1)) - p(:,N(2));
   Vec2 = p(:,N(3)) - p(:,N(2));
   Area(m) = norm(cross(Vec1,Vec2)) / 2;
   Center(:,m) = 1/3 * sum(p(:,N),2);
end

Edge_ = [];
n = 0;
for m = 1:TrianglesTotal
    N = t(1:3,m);
    for k = m + 1:TrianglesTotal
        M = t(1:3,k);      
        a = 1 - all([N-M(1) N-M(2) N-M(3)]);
        if (sum(a) == 2)
            n = n+1;
            Edge_ = [Edge_ M(find(a))]; 
            TrianglePlus(n) = m;
            TriangleMinus(n) = k; 
        end; 
    end
end
EdgesTotal = length(Edge_);

Edge__ = [Edge_(2,:); Edge_(1,:)];
Remove = [];
for m = 1:EdgesTotal
    Edge_m = repmat(Edge_(:,m),[1 EdgesTotal]);
    Ind1 = any(Edge_  -Edge_m);
    Ind2 = any(Edge__ -Edge_m);
    A = find(Ind1.*Ind2 == 0);
    if (length(A) == 3) 
        Out = find(t(4,TrianglePlus(A)) == t(4,TriangleMinus(A)));
        Remove = [Remove A(Out)];
    end
end
Edge_(:,Remove) = [];
TrianglePlus(Remove) = [];
TriangleMinus(Remove) = [];
EdgesTotal = length(Edge_);

EdgeIndicator = t(4,TrianglePlus)+t(4,TriangleMinus);

for m = 1:EdgesTotal
   EdgeLength(m) = norm(p(:,Edge_(1,m)) - p(:,Edge_(2,m)));
end

IMT = [];
for m = 1:TrianglesTotal
    n1 = t(1,m);
    n2 = t(2,m);
    n3 = t(3,m); 
    M = Center(:,m);
    r1 = p(:,n1);
    r2 = p(:,n2);
    r3 = p(:,n3);
    r12 = r2 - r1;
    r23 = r3 - r2;
    r13 = r3 - r1;
    C1 = r1 + (1/3) * r12;
    C2 = r1 + (2/3) * r12;
    C3 = r2 + (1/3) * r23;
    C4 = r2 + (2/3) * r23;
    C5 = r1 + (1/3) * r13;
    C6 = r1 + (2/3) * r13;
    a1 = 1/3 * (C1 + C5 + r1);
    a2 = 1/3 * (C1 + C2 + M);
    a3 = 1/3 * (C2 + C3 + r2);
    a4 = 1/3 * (C2 + C3 + M);
    a5 = 1/3 * (C3 + C4 + M);
    a6 = 1/3 * (C1 + C5 + M);
    a7 = 1/3 * (C5 + C6 + M);
    a8 = 1/3 * (C4 + C6 + M);
    a9 = 1/3 * (C4 + C6 + r3);
    Center_(:,:,m) = [a1 a2 a3 a4 a5 a6 a7 a8 a9];
end

for m = 1:EdgesTotal
    NoPlus = TrianglePlus(m);
    n1 = t(1,NoPlus);
    n2 = t(2,NoPlus);
    n3 = t(3,NoPlus); 
    if((n1 ~= Edge_(1,m)) & (n1 ~= Edge_(2,m))) NODE = n1; end;
    if((n2 ~= Edge_(1,m)) & (n2 ~= Edge_(2,m))) NODE = n2; end;
    if((n3 ~= Edge_(1,m)) & (n3 ~= Edge_(2,m))) NODE = n3; end;
    FreeVertex = p(:,NODE);
    
    RHO_Plus(:,m) =+Center(:,NoPlus) - FreeVertex;
    RHO__Plus(:,:,m) =+Center_(:,:,NoPlus) - repmat(FreeVertex,[1 9]);
end

for m = 1:EdgesTotal
    NoMinus = TriangleMinus(m);
    n1 = t(1,NoMinus);
    n2 = t(2,NoMinus);
    n3 = t(3,NoMinus); 
    if((n1 ~= Edge_(1,m)) & (n1 ~= Edge_(2,m))) NODE = n1; end;
    if((n2 ~= Edge_(1,m)) & (n2 ~= Edge_(2,m))) NODE = n2; end;
    if((n3 ~= Edge_(1,m)) & (n3 ~= Edge_(2,m))) NODE = n3; end;
    FreeVertex = p(:,NODE);
    
    RHO_Minus(:,m) =-Center(:,NoMinus) + FreeVertex;
    RHO__Minus(:,:,m) =-Center_(:,:,NoMinus) + repmat(FreeVertex,[1 9]);
end

save mesh2 p t TrianglesTotal EdgesTotal Edge_ TrianglePlus TriangleMinus ...
     EdgeLength EdgeIndicator Area RHO_Plus RHO_Minus RHO__Plus RHO__Minus ...
     Center Center_

epsilon_ = 8.854e-012;
mu_ = 1.257e-006;
c_ = 1 / sqrt(epsilon_ * mu_);
eta_ = sqrt(mu_ / epsilon_);
omega = 2 * pi * f;                                            
k = omega / c_;
K = j * k;
Constant1 = mu_ / (4 * pi);
Constant2 = 1 / (j * 4 * pi * omega * epsilon_);
Factor = 1 / 9;    
FactorA = Factor * (j * omega * EdgeLength / 4) * Constant1;
FactorFi = Factor * EdgeLength * Constant2;

for m = 1:EdgesTotal
    RHO_P(:,:,m) = repmat(RHO_Plus(:,m),[1 9]);   
    RHO_M(:,:,m) = repmat(RHO_Minus(:,m),[1 9]); 
end
FactorA = FactorA.';
FactorFi = FactorFi.';

Z = impmet( EdgesTotal,TrianglesTotal,...
            EdgeLength,K,...
            Center,Center_,...
            TrianglePlus,TriangleMinus,...
            RHO_P,RHO_M,...
            RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);   

save impedance f omega mu_ epsilon_ c_ eta_ Z  

for m = 1:EdgesTotal
    V(m) = 0;
    Distance(:,m) = 0.5 * sum(p(:,Edge_(:,m)), 2) - FeedPoint;
end

[Y,INDEX] = sort(sum(Distance.*Distance));
Index = INDEX(1);                

V(Index) = 1 * EdgeLength(Index);    

I = Z \ V.';

GapCurrent = sum(I(Index).*EdgeLength(Index)');
GapVoltage = mean(V(Index)./EdgeLength(Index));
Impedance = GapVoltage / GapCurrent
FeedPower = 1/2 * real(GapCurrent * conj(GapVoltage))

save current f omega mu_ epsilon_ c_ eta_ I V GapCurrent GapVoltage Impedance FeedPower

clear;
load('mesh2');
load('current');

Index = find(t(4,:) <= 1);
Triangles = length(Index);

for k = 1:Triangles
    i = [0 0 0]';
    for m = 1:EdgesTotal
        IE = I(m) * EdgeLength(m);
        if(TrianglePlus(m) == k)
            i = i + IE * RHO_Plus(:,m) / (2 * Area(TrianglePlus(m)));
        end
        if(TriangleMinus(m) == k)
            i = i + IE * RHO_Minus(:,m) / (2*Area(TriangleMinus(m)));
        end
    end
    CurrentNorm(k) = abs(norm(i));
end

Jmax = max(CurrentNorm);
MaxCurrent = strcat(num2str(Jmax),'[A/m]')
CurrentNorm1 = CurrentNorm/max(CurrentNorm);
for m = 1:Triangles
    N = t(1:3,m);
    X(1:3,m) = [p(1,N)]';
    Y(1:3,m) = [p(2,N)]';
    Z(1:3,m) = [p(3,N)]';      
end
C = repmat(CurrentNorm1,3,1);
h = fill3(X, Y, Z, C);
colormap gray;
axis('equal');
rotate3d