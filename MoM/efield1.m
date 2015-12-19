clear all
%Load the data
load('mesh2');
load('current');

k=omega/c_;
K=j*k;

for m=1:EdgesTotal
    Point1=Center(:,TrianglePlus(m));
    Point2=Center(:,TriangleMinus(m));
    DipoleCenter(:,m)=0.5*(Point1+Point2);
    DipoleMoment(:,m)=EdgeLength(m)*I(m)*(-Point1+Point2); 
end

ObservationPoint=[0; 0; 100];
[E,H]=point(ObservationPoint,eta_,K,DipoleMoment,DipoleCenter);

%find the sum of all dipole contributions
EField=sum(E,2); HField=sum(H,2);

%Common
EField                  %Radiated/scattered electric field 
                        %(complex vector at a point, V/m)

HField                  %Radiated/scattered magnetic field 
                        %(complex vector at a point, A/m)            

Poynting=0.5*real(cross(EField,conj(HField)))           
                        %Poynting vector (W/m^2) for radiated/scattered field

W=norm(Poynting)        %Radiation density (W/m^2) for radiated/scattered field
   
U=norm(ObservationPoint)^2*W                            
                        %Radiation intensity (W/unit solid angle)                     

%Only scattering
RCS=4*pi*(norm(ObservationPoint))^2*sum(EField.*conj(EField));     
                        %Backscattering radar cross-section (scattering)


         
