c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'          ! UDIFF, UTRANS

      integer ix,iy,iz,ieg

      UDIFF =0.
      UTRANS=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'SOLN'            ! jp
      include 'NEKUSE'          ! FFX, FFY, FFZ
      include 'PARALLEL'        ! GLLEL

      integer ix,iy,iz,ieg,e

      e = gllel(ieg)

!----------------------------------------------------------------------
!     Build the forcing for the perturbation evolution equation
!
!     du_i/dt = L_{NS} (u_i) - < u_j , L_{NS}(u_i) > u_j
!
! =>  du_i/dt = L_{NS} (u_i) - Lr_{ij}*u_j
!     \____________________/ \___________/
!                  |               |
! reg. pert. eq. __|               |______ user defined forcing FFi
!
      FFX = 0.0 
      FFY = 0.0 
      FFZ = 0.0 
     
      call set_OTD_forcing(FFX,FFY,FFZ,ix,iy,iz,ieg)

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'          ! QVOL

      integer ix,iy,iz,ieg

      QVOL   = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,eg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y
      include 'SOLN'            ! jp

      integer ix,iy,iz,iside,eg

      if (jp.eq.0) then
c     baseflow
        if (cbu.eq.'v  ') then
          UX = 0.0
          UY = 0.0
          UZ = 0.0
        endif
      else
c     perturbations
        UX = 0.0
        UY = 0.0
        UZ = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, Z
      include 'SOLN'            ! JP
      include 'PARALLEL'        ! gllel
      include 'OTD'

      integer ix,iy,iz,ieg

! IC check
      real ubic(lx1,ly1,lz1,lelv)
     $  ,  vbic(lx1,ly1,lz1,lelv)
     $  ,  wbic(lx1,ly1,lz1,lelv)
     $  ,  upic(lx1*ly1*lz1*lelv,lpert)
     $  ,  vpic(lx1*ly1*lz1*lelv,lpert)
     $  ,  wpic(lx1*ly1*lz1*lelv,lpert)
      common /bpic/ ubic,vbic,wbic
     $  ,           upic,vpic,wpic

      integer ie, ijke

      real amp, pi, kx(lpert), ky(lpert) 

      real mth_rand
      real xl(LDIM)
      real fcoeff(3)

      ie = gllel(ieg)

      pi = 4.0*atan(1.0)

c     velocity
c     base flow
      if (jp.eq.0) then
        UX = (1.0-Y**2)
        UY = 0.0
        UZ = 0.0
        
        ubic(ix,iy,iz,ie) = UX
        vbic(ix,iy,iz,ie) = UY
        wbic(ix,iy,iz,ie) = UZ
      else
c     perturbation
        ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(ie-1)))
        UX = vxpic(ijke,jp)
        UY = vypic(ijke,jp)
        UZ = vzpic(ijke,jp)
        temp = 0.0

        upic(ijke,jp) = UX
        vpic(ijke,jp) = UY
        wpic(ijke,jp) = UZ
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'            ! NX1, NY1, NZ1, NELV, NIO
      include 'FRAMELP'
      include 'TSTEP'           ! ISTEP, IOSTEP, TIME, LASTEP
      include 'SOLN'            ! V[XYZ], V[XYZ]P, PRP, JP, VMULT
      include 'MASS'            ! BM1
      include 'OTD'

! IC check
      real ubic(lx1,ly1,lz1,lelv)
     $ ,   vbic(lx1,ly1,lz1,lelv)
     $ ,   wbic(lx1,ly1,lz1,lelv)
     $ ,   upic(lx1*ly1*lz1*lelv,lpert)
     $ ,   vpic(lx1*ly1*lz1*lelv,lpert)
     $ ,   wpic(lx1*ly1*lz1*lelv,lpert)
      common /bpic/ ubic,vbic,wbic
     $ ,            upic,vpic,wpic

      integer     i
      character*1 str
      character*3 oname

!---------------------------------------------------------
!    Setup framework
!
      if (ISTEP.eq.0) then
        call frame_start
      endif

!---------------------------------------------------------
!    Setup monitoring
!
      call frame_monitor

!     Outpost initial conditions
      if (istep.eq.0) then
        call outpost(ubic,vbic,wbic,pr,t,'ibf')
        do i=1,npert
          write(str,'(I1)') i
          oname = 'ip'//trim(str)
          call outpost(upic(1,i),vpic(1,i),wpic(1,i),pr,t,oname)
        enddo
      endif

!-----------------------------------------------------------------------
!     Run OTD module 
!
      call run_OTD

!-----------------------------------------------------------------------
!     finalise framework
!
      if (ISTEP.eq.NSTEPS.or.LASTEP.eq.1) then
         call frame_end
      endif

      return
      end

c-----------------------------------------------------------------------
c     This routine to modify element vertices
      subroutine usrdat

      implicit none

      include 'SIZE'
      include 'TOTAL'

      real glmin,glmax
           
      real fact,x_min,x_max,y_max,y_min,z_min,z_max
      integer n

      n = 8*nelv
      
      fact = 4.*atan(1.)
      call cmult(xc,fact,n)
      if (if3d) then    
        call cmult(zc,fact,n)
      endif     

      x_min=glmin(xc,n)
      y_min=glmin(yc,n)
      x_max=glmax(xc,n)
      y_max=glmax(yc,n)
      if (if3d) then
        z_min=glmin(zc,n)
        z_max=glmax(zc,n)
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      implicit none
      
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      implicit none

      return
      end
!=======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     register modules
      call io_register
      call otd_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     initialise modules
      call otd_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------

      
      return
      end subroutine

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
