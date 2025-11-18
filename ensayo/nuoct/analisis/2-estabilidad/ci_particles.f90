!!!   ci_particles.f90  :  
!!!      
!!!   Condiciones iniciales para particulas (ncorp.f90). 
!!!   Chequear que tipo de coord/elementos pide el ncorp.in . 
!!!    
  implicit real*8 (a-h,k-z)
  integer npl,nx,ny
  real*8 xran(6)
  !
  twopi  = 8.0d0*datan(1.0d0)
  mearth = 3.04043d-6       ! earth mass in solar masses
  !
  open (3,file='particles.in',status='replace')

!!! array de numeros aleatorios xran(2).
  call random_seed ()
  call random_number (xran)

!!! minimum and maximum semimajor axis.
  amin = 1.276
  amax = 2.679
  
!!! minimum and maximum eccentricity.
  emin = 0
  emax = 0
  
!!! number of initial conditions in each orbital element.
  nx = 100
  ny = 1
  
!!! initial conditions for all particles.
  do ix = 1,nx
    do iy = 1,ny
      a = amin
      e = emin
      if (nx > 1) a = amin + (amax-amin)*dfloat(ix-1)/dfloat(nx-1)
      if (ny > 1) e = emin + (emax-emin)*dfloat(iy-1)/dfloat(ny-1)
      write (3,100) a,e,0.0d0,112.5,350.0,0.0d0
    end do
  end do
!
 100  format (1p8d15.4)
!
end program


function ran0 (iseed)
  implicit real*8 (a-h,k-z)
  implicit integer*4 (i-j)
  parameter (ia=48271,im=2147483647,iq=44488,ir=3399)
  parameter (scale=1./im,eps=1.2e-7,rnmx=1.-eps)
  !
  j = iseed/iq
  iseed = ia*(iseed-j*iq) - ir*j
  if (iseed < 0) iseed = iseed + im
  ran0 = min(iseed*scale,rnmx)
  !
end function ran0
