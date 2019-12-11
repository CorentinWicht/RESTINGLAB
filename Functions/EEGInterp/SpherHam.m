
function[Ylm] = SpherHam (L,M, Pos)

[Theta,Phi]=cart2sph(Pos(:,1),Pos(:,2),Pos(:,3));

% Part of this routine taken from Real-Valued spherical harmonics by Bing
% Jian

Plm = legendre(L,cos(Theta));

if L~=0
  Plm = squeeze(Plm(norm(M)+1,:,:));
end

a1 = ((2*L+1)/(4*pi));
a2 = factorial(L-norm(M))/factorial(L+norm(M));
C = sqrt(a1*a2);

if M <=-1
    
Ylm=(-1)^M*sqrt(2)*C*Plm'.*sin(norm(M)*Phi);

else 
if M>0
 Ylm=(-1)^M*sqrt(2)*C*Plm'.*cos(M*Phi);

else 
    Ylm=sqrt(a1)*Plm';
end
 
end
    





