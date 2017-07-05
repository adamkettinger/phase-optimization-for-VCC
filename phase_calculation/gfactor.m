function [out]=gfactor(cmap,af)

% g-factor calculation in SENSE reconstruction, as used in the manuscript.
%
% Input:    coil sensitivity matrix 'CMAP'
%           acceleration factor 'AF'
% Output:   g-Faktor values of the overlapping voxel group 'G'
%

[nc,ylines,nx]=size(cmap);

ny=ylines/af;

for y=1:ny
    for x=1:nx
        C=cmap(:,y:ny:af*ny,x);    
        
        bild=sqrt(diag(C'*C).*diag(pinv(C'*C)));
        
        g(y:ny:af*ny,x)=bild;
        
    end
end
g(g<1) = 1;       
out=abs(g);