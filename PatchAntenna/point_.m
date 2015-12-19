function[E]=point_(ObservationPoint,K,k,Constant2,Area,h,DipoleCenter)
%POINT_ Radiated electric field of volume current element(s) 
%   Replaces the volume current element by the equivalent 
%   finite-length dipole
%
%   Observation point                       ObservationPoint(1:3)         
%   Cross-section of the volume element     Area(1:M) 
%   Center of the volume element            DipoleCenter(1:3,1:M)
%   Height of the volume element            h
%   E-field at the observation point        E
%
%   Copyright 2002 AEMM. Revision 2002/03/17 
%   Chapter 10

c       =DipoleCenter;
r       =repmat(ObservationPoint,[1 length(c)])-c(1:3,:);
R       =sqrt(sum(r.*r));
G       =exp(-K*R)./R;

EZ=Constant2*k^2*h*Area.*G;

r       =repmat(ObservationPoint+[0; 0; h/2],[1 length(c)])-c(1:3,:);
R       =sqrt(sum(r.*r));
OverR   =1./R;
G       =exp(-K*R).*OverR; 

Common=Constant2*Area.*G.*(OverR.*OverR+K.*OverR);
EXPlus=Common.*(-r(1,:));
EYPlus=Common.*(-r(2,:));
EZPlus=Common.*(-r(3,:));

r       =repmat(ObservationPoint+[0; 0; -h/2],[1 length(c)])-c(1:3,:);
R       =sqrt(sum(r.*r));
OverR   =1./R;
G       =exp(-K*R).*OverR; 

Common=Constant2*Area.*G.*(OverR.*OverR+K.*OverR);
EXMinus=Common.*(-r(1,:));
EYMinus=Common.*(-r(2,:));
EZMinus=Common.*(-r(3,:));

E(1,:)=EXPlus-EXMinus;
E(2,:)=EYPlus-EYMinus;
E(3,:)=EZPlus-EZMinus+EZ;

