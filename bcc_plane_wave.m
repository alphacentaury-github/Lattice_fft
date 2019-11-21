L=5;
L3=L**3;
nr = 0:(L3-1);
lx = mod(nr,L);
ly = mod((nr - lx)/L,L);
lz = (nr -lx -ly*L)/L**2; 
lx_bcc(1:L3) = lx ;
lx_bcc(L3+1:2*L3) =lx;
ly_bcc(1:L3) = ly ;
ly_bcc(L3+1:2*L3) =ly;
lz_bcc(1:L3) = lz ;
lz_bcc(L3+1:2*L3) =lz;
parity_bcc(1:L3) = 0;
parity_bcc(L3+1:2*L3) = 1;
nx = 2*lx_bcc +parity_bcc ;
ny = 2*ly_bcc +parity_bcc ;
nz = 2*lz_bcc +parity_bcc ;
x_bcc(1:L3) = mod(lx+L/2,L)-L/2;
x_bcc(L3+1:2*L3) =x_bcc(1:L3)+0.5;
y_bcc(1:L3) = mod(ly+L/2,L)-L/2;
y_bcc(L3+1:2*L3) =y_bcc(1:L3)+0.5;
z_bcc(1:L3) = mod(lz+L/2,L)-L/2;
z_bcc(L3+1:2*L3) =z_bcc(1:L3)+0.5;

dx = x_bcc; dy=y_bcc; dz = z_bcc;
dr2 = dx.^2+dy.^2+dz.^2;
dr =sqrt(dr2); 
[dr_sorted  idx] = sort(dr) ;
for i=1:66
  k=idx(i);
  wave1(i,nx,ny,nz) = exp(1j*2*pi/(2*L)*(nx(k)*nx+ny(k)*ny+nz(k)*nz));
endfor   