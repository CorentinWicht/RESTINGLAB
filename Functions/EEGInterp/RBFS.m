function [Eval,l] = RBFS(M,RBF,c)
%RBFS.m
%
%Computes the values of a specific radial basis functions at given values or a distance matrix
%
%INPUT::
% M := Points to evaluate the RBF usually a distance matrix
% RBF := Which radial basis function to be used, possible choices:
%   'TPS':= Thin-plate-splines / 'r^3' := Konst^3 
%   'MQ':= Multiquadric/ 'IMQ':= Inverse MQ / 'Gau':= Gaussian
% c := Parameter for the radial basis functions (only needed for TPS MQ IMQ Gau)
% For Gaussian c>0 
%
%OUTPUT:
% Eval := Evaluation of the basis function in M
% l    := Order of spherical harmonics which need to be added to the interpolation
%

switch RBF




%Shifted Thin-Plate-Splines
    case {'TPS'}
if c==0 
	for j=1:size(M,1)
	 for n=1:size(M,2)
	 	if M(j,n)==1
		M(j,n)=0.5; 
		end
	 end
  end; 
  end;
Eval=(2-2.*M+c^2).*log(2-2.*M+c^2);
l=2;


%Linear
case{'r^3'}
Eval=M.^3;
l=0;


%Multiquadric
    case{'MQ'}
Eval=sqrt(2-2.*M+c^2);
l=0;


%Inverse Multiquadric
    case{'IMQ'}
Eval=(2-2.*M+c^2).^(-0.5);
l=0;



%Gaussian
    case{'Gau'}
Eval=exp(-c.*(2-2.*M));
l=0;

end
end



