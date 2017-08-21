function s = simpson3d( z, dx, dy, dz )
[nx ny nz] = size(z);
nx = nx-1;
ny = ny-1;
nz = nz-1;
s1 = ( z(1,1,1) + z(1,1,nz+1) + z(nx+1,1,1) + z(1,ny+1,1)+ z(nx+1,ny+1, nz+1)+ z(1,ny+1,nz+1) + z(nx+1,ny+1,1) + z(nx+1,1,nz+1));
ixo = 2:2:nx;
ixe = 3:2:nx-1;
iyo = 2:2:ny;
iye = 3:2:ny-1;
izo = 2:2:nz;
ize = 3:2:nz-1;
s2 = 2*( sum(z(1,iye,1)) + sum(z(ixe,1,1)) +sum(z(1,1,ize)) + sum(z(nx+1,iye,1)) +...
    sum(z(1,iye,nz+1)) + sum(z(nx+1,iye, nz+1))+ sum(z(1,ny+1,ize))+ sum(z(nx+1,1,ize)) +...
    sum(z(nx+1,ny+1, ize)) +sum(z(ixe,ny+1,1)) + sum(z(ixe,1,nz+1)) +sum(z(ixe,ny+1, nz+1)));
s3 = 4*( sum(z(1,iyo,1)) + sum(z(ixo,1,1)) +sum(z(1,1,izo)) + sum(z(nx+1,iyo,1)) +...
    sum(z(1,iyo,nz+1)) + sum(z(nx+1,iyo, nz+1))+ sum(z(1,ny+1,izo))+ sum(z(nx+1,1,izo)) +...
    sum(z(nx+1,ny+1, izo)) +sum(z(ixo,ny+1,1)) + sum(z(ixo,1,nz+1)) +sum(z(ixo,ny+1, nz+1)));

s4 = 64*sum(sum( sum( z(ixo,iyo,iyo) ) )) + 8*sum(sum( sum( z(ixe,iye, ize) )) );
s5 =  32*sum(sum( sum( z(ixe,iyo, izo) ) )) + 32*sum( sum(sum( z(ixo,iye, izo) ) )) + 16*sum(sum( sum( z(ixe,iye, izo)) ) ) + 16*sum( sum(sum( z(ixo,iye, ize) ) )) + 16*sum(sum( sum( z(ixe,iyo, ize) )) ) + 32*sum( sum(sum( z(ixo,iyo, ize) )) );
s = s1 + s2 + s3 + s4 + s5;
s = s*dx*dy*dz/27.0;
end
