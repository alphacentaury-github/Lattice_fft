# lattice dispersion relation.
# Cubic lattice and BCC lattice .  
# momentum p=2*pi/La *(nx,ny,nz) 
# Note here that La is a length. 
clear

function [ plist, psqr] = gen_momentum(L)
  # Generate momentum in a cubic box L 
  # px= 2*pi/L*nx  with nx in [-L,L-1] range 
  # sort them in increasing momentum square 
  # Note: no periodic condition is checked. 
  # construct points in cubic [0,L] 
  nn = 0: L^3-1 ;
  nx = mod(nn,L);           # range [0,L-1]
  ny = mod((nn - nx)/L, L);
  nz = (nn - nx -ny*L)/L^2 ; 
  rx = mod(nx +L/2,L) - L/2; # range [-L/2,L/2-1]   
  ry = mod(ny +L/2,L) - L/2;   
  rz = mod(nz +L/2,L) - L/2;   
  plist = [rx' ry' rz'];  
  psqr = rx.^2+ry.^2+rz.^2;
  # sort 
  [sorted,ind] = sort(psqr) ;
  # sorted results 
  plist = plist(ind,:) ; 
  psqr = sorted ; 
endfunction  

function psqr = kin_energy(l_type,korder,nx,ny,nz,La)
# compute free particle energy of one momentum 
# using lattice dispersion relation 
# Input:
#   l_type = type of lattice 'free', 'cubic' or 'bcc' 
#   korder = order of improved action. 
#          For the moment, only 3 or 4 is defined. 
#   (nx,ny,nz) = integer value corresponding to momentum 
#   La = lattice size in lattice unit. La= L*a 
#      momentum is defined as
#      (px,py,pz)=(2*pi)/(L*a)*(nx,ny,nz) 
# Output: 
#   psqr = p**2 value in lattice unit 
  if strcmp(l_type,'free')
     p = 2*pi/La *[nx ny nz] ;
     psqr = p(1)^2+p(2)^2+p(3)^2;
     return 
  endif 
  if strcmp(l_type,'cubic')
    if korder==3
      c0=49/6;
      c = [-3/2 3/20 -1/90];
      p = 2*pi/La *[nx ny nz];
      psqr = c0 + c(1)*2*(cos(p(1))+cos(p(2))+cos(p(3)))+...
            +c(2)*2*(cos(2*p(1))+cos(2*p(2))+cos(2*p(3)))+...
            +c(3)*2*(cos(3*p(1))+cos(3*p(2))+cos(3*p(3)));
      return       
    endif 
    if korder==4
      c0=205/24;
      c = [-8/5 1/5 -8/315 1/560];
      p = 2*pi/La *[nx ny nz];
      psqr = c0 + c(1)*2*(cos(p(1))+cos(p(2))+cos(p(3)))+...
            +c(2)*2*(cos(2*p(1))+cos(2*p(2))+cos(2*p(3)))+...
            +c(3)*2*(cos(3*p(1))+cos(3*p(2))+cos(3*p(3)))+...
            +c(4)*2*(cos(4*p(1))+cos(4*p(2))+cos(4*p(3)));
      return       
    endif   
  endif   
  if strcmp(l_type,'bcc')
    if korder==3
      c0= 98/9;
      c = [-3/2 3/20 -1/90];
      p = 2*pi/La *[nx ny nz];
      psqr = c0 + c(1)*8*cos(p(1)/2)*cos(p(2)/2)*cos(p(3)/2)+...
            + c(2)*8*cos(p(1))*cos(p(2))*cos(p(3))+...
            + c(3)*8*cos(3/2*p(1))*cos(3/2*p(2))*cos(3/2*p(3));
      return       
    endif 
    if korder==4 
      c0= 205/18;
      c = [-8/5 1/5 -8/315 1/560];
      p = 2*pi/La *[nx ny nz];
      psqr = c0 + c(1)*8*cos(p(1)/2)*cos(p(2)/2)*cos(p(3)/2)+...
            + c(2)*8*cos(p(1))*cos(p(2))*cos(p(3))+...
            + c(3)*8*cos(3/2*p(1))*cos(3/2*p(2))*cos(3/2*p(3))+...
            + c(4)*8*cos(2*p(1))*cos(2*p(2))*cos(2*p(3));
      return       
    endif         
  endif 
  printf('Error: Unknown option!')   
endfunction 

# Energy of free N Fermion Gas in a box L 
# N fermion is N/2 spin up and N/2 spin down 
cutoff = 100 ;% MeV 
a=1/cutoff ; % MeV^-1
tcutoff = 1000; % MeV 
at = 1/tcutoff;
atovera = at/a ; % unit of a 
for L=4:11
%L=5; 
N=66;
NMass= 938.92/cutoff ; %a unit  
h0=1./(2*NMass) ; % a unit 
[ plist, psqr] = gen_momentum(L) ; 

%-----because no momentum equivalence is considered
%     following energy calculation is only valid when
%     number of particle is sufficiently smaller than L^3. 
sum_free= 0;         
sum_cubic= 0;
sum_cubic2=0;
sum_bcc= 0;
sum_bcc2=0;
La = L*a ;
sum_bcc_tr = 0;
sum_cubic2_tr = 0;
for i = 1:N/2
  nn = plist(i,:);
  sum_free = sum_free + 2*kin_energy('free',3,nn(1),nn(2),nn(3),L);
  sum_cubic = sum_cubic + 2*kin_energy('cubic',3,nn(1),nn(2),nn(3),L);
  sum_cubic2 = sum_cubic2 + 2*kin_energy('cubic',4,nn(1),nn(2),nn(3),L);
  sum_cubic2_tr = sum_cubic2_tr - 1/atovera*2*log(1-h0*kin_energy('cubic',4,nn(1),nn(2),nn(3),L)*atovera);
  sum_bcc = sum_bcc + 2*kin_energy('bcc',3,nn(1),nn(2),nn(3),L);
  sum_bcc2 = sum_bcc2 + 2*kin_energy('bcc',4,nn(1),nn(2),nn(3),L);
  sum_bcc_tr = sum_bcc_tr - 1/atovera*2*log(1-h0*kin_energy('bcc',3,nn(1),nn(2),nn(3),L)*atovera);
  E_free(i) = h0*sum_free*cutoff ; %convert to MeV unit 
  E_cubic(i) = h0*sum_cubic*cutoff;
  E_cubic2(i) = h0*sum_cubic2*cutoff;
  E_bcc(i) = h0*sum_bcc*cutoff ;
  E_bcc2(i) = h0*sum_bcc2*cutoff ;
  E_bcc_tr(i) = sum_bcc_tr*cutoff;  
  E_cubic2_tr(i) = sum_cubic2_tr*cutoff;  
endfor   

#-----Thermodynamic limit expression------------------
Nlist = [ 14 38 54 66 114];
EfgN14 = h0*(2*pi/(L))**2*6*2;
EfgN38 = h0*(2*pi/(L))**2*(6*2+2*12*2);
EfgN54 = h0*(2*pi/(L))**2*(6*2+2*12*2+3*8*2);
EfgN66 = h0*(2*pi/(L))**2*(6*2+2*12*2+3*8*2+4*6*2);
EfgN114 =h0*(2*pi/(L))**2*(6*2+2*12*2+3*8*2+4*6*2+5*24*2);
Efglist = [EfgN14 EfgN38 EfgN54 EfgN66 EfgN114];
for i = 1:length(Nlist)
  NN=Nlist(i); 
  Efree(i) = Efglist(i)*cutoff; 
  kf = (3*pi**2*NN)**(1./3.)/(L); %MeV  
  EfgTherm(i) = 3*NN*kf**2/(10*NMass)*cutoff ; % MeV 
  gamma(i) = abs(Efree(i)-EfgTherm(i))/EfgTherm(i) ;
endfor 

%printf('****** N=66, L= %i \n',L) 
%printf('     exact   cubic    cubic2    bcc      bcc2     FG \n')
%disp([E_free(33) E_cubic(33) E_cubic2(33) E_bcc(33) E_bcc2(33) EfgTherm(4)])

printf('N    L     exact   bcc      bcc_tr     FG  cubic_3 \n')
printf('%i %i %f  %f  %f   %f\n',N,L,E_free(33), E_bcc(33),E_bcc_tr(33), EfgTherm(4)...
                             ,E_cubic(33)  )
#disp([N  L  E_free(33) E_bcc(33)   E_bcc_tr(33) EfgTherm(4)])
          
endfor 

# Unitary Fermion Gas in BCC lattice 
