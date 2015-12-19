function []=viewer(filename)
%VIEWER Visualizes the structure - all Chapters
%   Modification of Chapter 10: different antenna components
%   are shown using different colors
%   
%   Usage:  viewer('patch.mat') or
%           viewer('patch') or
%           viewer patch
%
%   Copyright 2002 AEMM. Revision 2002/03/22 Chapter 10

FeedingTriangle=[];
load(filename);

%Format:
%t(4,:)=0 triangles of metal ground plane
%t(4,:)=1 triangles of probe feed
%t(4,:)=2 triangles of the patch
%t(4,:)=3 triangles of the upper boundary of dielectric (not seen)

Color1=[1.00 1.00 1.00];    
Color2=[0.65 0.65 0.65];
Color3=[0.25 0.25 0.25];
Color4=[1 1 1];

TrianglesTotal=length(t)
MetalPlate      =find(t(4,:)==0);
MetalPatch      =find(t(4,:)==2);
Feed            =find (t(4,:)==1);
Feed            =[Feed FeedingTriangle];
Dielectric      =find (t(4,:)==3);

for k=1:length(MetalPlate)
    m=MetalPlate(k);
    n=t(1:3,m);
    X(1:3,m)=p(1,n)';
    Y(1:3,m)=p(2,n)';
    Z(1:3,m)=p(3,n)';
end
for k=1:length(MetalPatch)
    m=MetalPatch(k);
    n=t(1:3,m);
    X_(1:3,m)=p(1,n)';
    Y_(1:3,m)=p(2,n)';
    Z_(1:3,m)=p(3,n)';
end
for k=1:length(Dielectric)
    m=Dielectric(k);
    n=t(1:3,m);
    X1(1:3,m)=p(1,n)';
    Y1(1:3,m)=p(2,n)';
    Z1(1:3,m)=p(3,n)';
end
for k=1:length(Feed)
    m=Feed(k);
    n=t(1:3,m);
    X2(1:3,m)=p(1,n)';
    Y2(1:3,m)=p(2,n)';
    Z2(1:3,m)=p(3,n)';
end

if(~isempty(Feed))
    g=fill3(X,Y,Z,Color2,X_,Y_,Z_,Color3,X2,Y2,Z2,Color4);
    %g=fill3(X,Y,Z,Color2,X1,Y1,Z1,Color1,X2,Y2,Z2,Color4);
elseif (~isempty(MetalPatch))
    g=fill3(X1,Y1,Z1,Color2,X_,Y_,Z_,Color3);
else
    g=fill3(X,Y,Z,Color2);
end

axis('equal')
view(148,34)
rotate3d on
