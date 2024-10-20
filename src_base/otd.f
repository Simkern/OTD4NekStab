!
!     OTD implementation by Simon Kern
!     Based on the method proposed in
!     Babaee,H. and Sapsis,T.P., "A minimization principle for the
!      description of modes associated with finite-time instabilities",
!      Proc. R. Soc. A 474: 20150779,
!      https://dx.doi.org/10.1098/rspa.2015.0779  
!      
!     Email:     skern@mech.kth.se      
!
!     List of subroutines:
!       FRAMEWORK:
!         otd_register
!         otd_init
!       MAIN INTERFACE:
!         run_otd
!       SETUP:
!         white_noise_IC 
!         read_OTDIC
!         OTD_ON
!       LINEARIZED OPERATOR:
!         gen_LU
!         proj_OTDmodes
!         gen_OTD_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)       -> userbc
!       OUTPOSTING
!         outpost_projmodes
!         outpost_OTDbasis
!       FTLEs
!         get_FTLE
!         reset_FTLE
!         compute_FTLE_Phi
!         compute_FTLE_Blanchard
!
!======================================================================

!======================================================================
!> @brief Register OTD module
!! @note This routine should be called in frame_usr_register

      subroutine otd_register()     
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'OTD'

      ! local variables
      integer lpmid
      real ltim

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,otd_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(otd_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'Parent module ['//'FRAME'//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(otd_id,lpmid,otd_name,
     $          'Optimally time-dependent (OTD) basis.')

      ! register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      ! total time
      call mntr_tmr_reg(otd_ttot_id,lpmid,otd_id,
     $     'OTD_TOT','OTD module total time',.false.)
      ! Initialisation time
      call mntr_tmr_reg(otd_tini_id,otd_ttot_id,otd_id,
     $     'OTD_INI','OTD initialisation time',.true.)
      ! Orthonormalisation time
      call mntr_tmr_reg(otd_tGS_id,otd_ttot_id,otd_id,
     $     'OTD_GS','Orthonormalisation time',.true.)
      ! LU & forcing building time
      call mntr_tmr_reg(otd_tgen_id,otd_ttot_id,otd_id,
     $     'OTD_gen','LU and foring computation time',.true.)
      ! OTD mode computation time
      call mntr_tmr_reg(otd_tget_id,otd_ttot_id,otd_id,
     $     'OTD_get','OTD mode computation time',.true.)
      ! FTLE computation time
      call mntr_tmr_reg(otd_tFTLE_id,otd_ttot_id,otd_id,
     $     'OTD_FLTE','FTLE computation time',.true.)
      ! IO
      call mntr_tmr_reg(otd_tIO_id,otd_ttot_id,otd_id,
     $     'OTD_IO','OTD module IO time',.true.)
       
      ! register and set active section
      call rprm_sec_reg(otd_sec_id,otd_id,
     $     '_'//adjustl(otd_name),
     $     'Runtime parameter section for FST module')
      call rprm_sec_set_act(.true.,otd_sec_id)

      ! register parameters
      ! otd_nusrIC
      call rprm_rp_reg(otd_nusric_id,otd_sec_id,'OTD_NUSRIC',
     $     'Number of IC fields to be read ',
     $     rpar_int,0,0.0,.false.,'  ')
      ! otd_IOstep
      call rprm_rp_reg(otd_iostep_id,otd_sec_id,'OTD_IOSTEP',
     $     'IO frequency for the OTD module ',
     $     rpar_int,5,0.0,.false.,'  ')
      ! otd_rststp
      call rprm_rp_reg(otd_rststp_id,otd_sec_id,'OTD_RSTSTP',
     $     'IO frequency for OTD basis (restart) ',
     $     rpar_int,1000,0.0,.false.,'  ')
      ! otd_prntsp
      call rprm_rp_reg(otd_prntsp_id,otd_sec_id,'OTD_PRNTSP',
     $     'Print frequency for the OTD module ',
     $     rpar_int,100,0.0,.false.,'  ')
      ! otd_GSstep
      call rprm_rp_reg(otd_gsstep_id,otd_sec_id,'OTD_GSSTEP',
     $     'Frequency of GS orthonormalization of OTD basis ',
     $     rpar_int,5,0.0,.false.,'  ')
      ! otd_FTLEpd
      call rprm_rp_reg(otd_FTLEpd_id,otd_sec_id,'OTD_FTLEpd',
     $     'Time horizon for FTLE computation ',
     $     rpar_real,0,0.0,.false.,'  ')
      ! otd_ifFTLE
      call rprm_rp_reg(otd_ifFTLE_id,otd_sec_id,'OTD_IFFTLE',
     $     'Compute FTLEs? ',
     $     rpar_log,0,0.0,.false.,'  ')
      ! otd_ifvlog
      call rprm_rp_reg(otd_ifvlog_id,otd_sec_id,'OTD_IFVLOG',
     $     'Verbose logging for OTD module ',
     $     rpar_log,0,0.0,.false.,'  ')
      ! otd_debug
      call rprm_rp_reg(otd_debug_id,otd_sec_id,'OTD_DEBUG',
     $     'Debugging mode for OTD module ',
     $     rpar_log,0,0.0,.false.,'  ')

      ! set initialisation flag
      otd_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_ttot_id,1,ltim)

      return
      end subroutine otd_register

!======================================================================
!> @brief Initilise OTD module
!! @note This routine should be called in frame_usr_init
!! @remark This routine uses global scratch space \a SCRUZ      

      subroutine otd_init()
      implicit none

      include 'SIZE'
      include 'INPUT'               ! param(59), initc
      include 'TSTEP'               ! time
      include 'GEOM'                ! [xyz]m1
      include 'FRAMELP'
      include 'OTD'

      ! local variables
      integer       itmp
      real          rtmp, ltim
      logical       ltmp
      character*20  ctmp
      character*2   str1, str2
      character*200 lstring
      real          xtmp(lx1,ly1,lz1,lelv)
      real          ytmp(lx1,ly1,lz1,lelv)
      real          ztmp(lx1,ly1,lz1,lelv)

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the module was already initialised
      if (otd_ifinit) then
         call mntr_warn(otd_id,
     $        'module ['//trim(otd_name)//'] already initiaised.')
         return
      endif

      ! get runtime parameters
! otd_nusric
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_nusric_id,rpar_int)
      otd_nusric = itmp
! otd_IOstep
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_iostep_id,rpar_int)
      otd_iostep = itmp
! otd_rststp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_rststp_id,rpar_int)
      otd_rststp = itmp
! otd_prntsp
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_prntsp_id,rpar_int)
      otd_prntsp = itmp
! otd_GSstep
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_gsstep_id,rpar_int)
      otd_gsstep = itmp
! otd_FTLEpd
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_FTLEpd_id,rpar_real)
      otd_FTLEpd = rtmp
! otd_ifFTLE
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_ifFTLE_id,rpar_log)
      otd_ifFTLE = ltmp
! otd_ifvlog
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_ifvlog_id,rpar_log)
      otd_ifvlog = ltmp
! otd_debug
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,otd_debug_id,rpar_log)
      otd_debug = ltmp

!---------------------------------------------------------
!    Initialise all IC fields to white noise
!
      call white_noise_IC()

!---------------------------------------------------------
!    Read IC fields 
!
      if (otd_nusric.lt.0) then
        call mntr_abort(otd_id, 'Choose a valid number of IC files.')
      else if (otd_nusric.eq.0) then
        call mntr_log(otd_id,lp_inf,
     $ 'No IC files read. ICs will be white noise.')
      else
        if (otd_nusric.gt.npert) then
          call mntr_abort(otd_id,
     $ 'otd_nusric > npert! Increase the number of perturbations.')
        endif
        ! save mesh
        call opcopy(xtmp,ytmp,ztmp,xm1,ym1,zm1)
        ! read ICs
        call read_OTDIC()
        ! restore mesh
        call opcopy(xm1,ym1,zm1,xtmp,ytmp,ztmp)
      endif

!---------------------------------------------------------
!    Sanity check
!
      if (otd_ifFTLE) then
        if (otd_FTLEpd.eq.0.0) then
          call mntr_warn(otd_id,
     $ 'The FLTEs will be computed continuously. Set otd_FTLEpd in the
     $ par-file to define a finite time horizon.')
        endif
        pcount = 0
      endif

!---------------------------------------------------------
!    Ensure right variable size
!
      if (npert.gt.lpert) then
        call mntr_abort(otd_id,
     $ 'npert > lpert! Set lpert=npert for the tool to work.')
      endif

!---------------------------------------------------------
!    Set ICs 
!
      call blank(initc,132)
      call setics
      call mntr_log(otd_id,lp_inf,'INIT - set ICs :: done')
      
!---------------------------------------------------------
!    Set default timestep at which to start OTD computation
!
!    we use the first timestep to make sure that the perturbations
!     1. are divergence free
!     2. satisfy the orthonormality constraint
!     3. satisfy the boundary conditions of the problem
!
      startOTD = 1

      ! everything is initialised
      otd_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tini_id,1,ltim)

      return
      end subroutine otd_init

!=======================================================================
!> @brief Check if module was initialised
!! @return otd_is_initialised

      logical function otd_is_initialised()
      implicit none

      include 'SIZE'
      include 'OTD'
!-----------------------------------------------------------------------
      otd_is_initialised = otd_ifinit

      return
      end function otd_is_initialised

!=======================================================================
!> @brief Main interface for the OTD module
!  
      subroutine run_OTD()
      implicit none

      include 'SIZE'
      include 'TSTEP'           ! istep,iostep
      include 'FRAMELP'
      include 'OTD'

      include 'MASS'

      real op_glsc2_wt
      integer n9,ntot,i
      real iz(lpert)

!----------------------------------------------------------------------

      if (istep.eq.0) then              ! OTD initialisation
!-----------------------------------------------------------------------
!       Perform orthonormalisation of the basis vectors
!       Initialise fundamental solution matrix for FTLE computation
!       Compute one time step 
!
        gsstep_override = .true.        ! initialisation
        call OTD_ON()
!        call outpost_OTDbasis()
        if (otd_ifFTLE) then
          call reset_FTLE
        endif
!
      elseif (istep.ge.startOTD) then   ! standard OTD

!-----------------------------------------------------------------------
!       Orthonormalize perturbations (if necessary)     
!
        gsstep_override = .false.
        call OTD_ON()

!-----------------------------------------------------------------------
!       Build the action of the linearized NS-operator and compute the
!       reduced operator Lr_{ij} 
!
        call gen_Lu()

!-----------------------------------------------------------------------
!       Compute additional forcing for perturbation equation 
!
        call gen_OTD_forcing()
        
!-----------------------------------------------------------------------
!       Compute the OTD modes
!
        if (mod(istep,otd_prntsp).eq.0) then
          call proj_OTDmodes()
        endif

!-----------------------------------------------------------------------
!       Output OTD modes
!
        if (mod(istep,otd_iostep).eq.0) then
          call outpost_projmodes()
        endif

!-----------------------------------------------------------------------
!       Output perturbation fields (OTD basis) for restart
!
        if (mod(istep,otd_rststp).eq.0
     $      .or.istep.eq.nsteps.or.lastep.eq.1) then
          call outpost_OTDbasis()
        endif

!-----------------------------------------------------------------------
!       FTLEs
!
        if (otd_ifFTLE) then
          call get_FTLE()
        endif
!
      endif     ! istep 0 --> init


      return
      end subroutine run_OTD

!=======================================================================
!> @brief Construct action of linearized NS operator on perturbation
!  field
!
      subroutine gen_Lu()
      implicit none
              
      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'FRAMELP'
      include 'MASS'            ! BM1
      include 'SOLN'            ! V[XYZ]P
      include 'TSTEP'           ! istep
      include 'OTD'             ! conv[xyz],diff[xyz],gradp[xyz]
                                ! Lu[xyz], Lr

      ! local variables
      integer ipert, jpert, ntot, i

      ! timing
      real ltim

      ! function
      real op_ip
      real dnekclock

!#define DEBUG
#ifdef DEBUG
      character*1 str
      character*3 oname
#endif

      ! timing
      ltim = dnekclock()

!----------------------------------------------------------------------
!     Build the elements of the linearized NS-operator L_{NS} (u_j)
!
!     L_{NS} (u_j) = 1/Re (grad^2 u)_j - (grad p)_j - (Ub.grad) u_j - (u_j.grad) Ub
!
      do ipert=1,npert

        ! Convective terms 
        call Lu_op_conv(vxp(1,ipert),vyp(1,ipert),vzp(1,ipert),ipert)

        ! Perturbation pressure gradient
        call Lu_op_gradp(prp(1,ipert),ipert)

        ! Diffusive term
        call Lu_op_diff(vxp(1,ipert),vyp(1,ipert),vzp(1,ipert),ipert)
      
      enddo

#ifdef DEBUG
! DIAGNOSTICS of the elements of LU_{NS}
      do i=1,npert
        if (mod(istep,otd_iostep).eq.0) then
          write(str,'(I1)') i
          oname = 'p'//trim(str)//'d'
          call outpost(diffx(1,i),diffy(1,i),diffz(1,i),pr,t,oname)
          oname = 'p'//trim(str)//'p'
          call outpost(gradpx(1,i),gradpy(1,i),gradpz(1,i),pr,t,oname)
          oname = 'p'//trim(str)//'c'
          call outpost(convx(1,i),convy(1,i),convz(1,i),pr,t,oname)
! To compute the individual convective terms, activate in Lu_op_conv
!          oname = 'c'//trim(str)//'1'
!          call outpost(convx1(1,i),convy1(1,i),convz1(1,i),pr,t,oname)
!          oname = 'c'//trim(str)//'2'
!          call outpost(convx2(1,i),convy2(1,i),convz2(1,i),pr,t,oname)
        endif
      enddo
#endif
#undef DEBUG

!----------------------------------------------------------------------
!     Assemble the action of the operator
!
!     L_{NS} (u_j) = 1/Re grad^2 u_j - grad p - (Ub.grad) u_j - (u_j.grad) Ub
!
      ntot = lx1*ly1*lz1*nelv
!      
      do jpert=1,npert
        do i=1,ntot
          Lux(i,jpert) = diffx(i,jpert)- gradpx(i,jpert)- convx(i,jpert)
          Luy(i,jpert) = diffy(i,jpert)- gradpy(i,jpert)- convy(i,jpert)
          if (if3d) then
            Luz(i,jpert)=diffz(i,jpert)- gradpz(i,jpert)- convz(i,jpert)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     Compute the innner product < L_{NS}(u_i),u_j > with i,j = 1,...,r
!
      call rzero(Lr,lpert*lpert)
      do ipert=1,npert
        do jpert=1,npert
          Lr(ipert,jpert) = op_ip(ipert,jpert,2)
        enddo
      enddo
!
! Debug output: Print reduced operator Lr_{ij}
!
      if (nid.eq.0.and.otd_debug) then
        if (mod(istep,otd_prntsp).eq.0) then
          ! write out Lr
          call outmat(Lr,npert,npert,'Lrmat  ',istep)
        endif    ! otd_prntsp
      endif      ! otd_debug.and.nid.eq.0

!----------------------------------------------------------------------
!     Here you could add internal rotations phi_rot into the method.
!     This does NOT change the subspace.
!
!       dU / dt = L_{NS}U - U (Lr - phi_rot)
!
      call rzero(phi_rot,lpert*lpert)
!
!     The rotation matrix Phi_ij must be skew-symmetric (but is otherwise
!     arbitrary
!
!       phi_rot(i,j) = -phi_rot(j,i)
!      
!     e.g. to obtain an evolution that corresponds to continuously
!     performing Gram-Schmidt on the basis (i.e. turning Lr into a lower
!     triangular matrix), set the rotation matrix to
! 
!                  / -<Lu_j,u_i>     j < i
!       phi_rot = {   0              j = i
!                  \  <Lu_j,u_i>     j > i
!
      if (npert.gt.1) then
        do jpert=1,npert
          do ipert=jpert+1,npert
            phi_rot(ipert,jpert) = op_ip(ipert,jpert,2)
            phi_rot(jpert,ipert) = -phi_rot(ipert,jpert)
          enddo
        enddo
      endif
      
      if (nid.eq.0.and.otd_debug) then
        if (mod(istep,otd_prntsp).eq.0) then
          call outmat(Lr,npert,npert,'Lr-mat',istep)
          call outmat(phi_rot,npert,npert,'phimat',istep)
        endif
      endif
      ! add internal rotation if defined
      call sub2(Lr,phi_rot,lpert*lpert)
      if (nid.eq.0.and.otd_debug) then
        if (mod(istep,otd_prntsp).eq.0) then
          call outmat(Lr,npert,npert,'Lrpmat',istep)
        endif
      endif

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tgen_id,1,ltim)

      return
      end subroutine gen_Lu

!=======================================================================
!> @brief Compute eigenspectrum of the reduced operator Lr_{ij} and
!         project the velocity perturbations onto the eigendirections to
!         obtain the most unstable modes
! 
      subroutine proj_OTDmodes()
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'FRAMELP'
      include 'TSTEP'           ! istep
      include 'SOLN'            ! v[xyz]p
      include 'WLAPACK'
      include 'OTD'             ! EIG[RI],EV[RL],Lr,OTDmr[xyz]

      ! local variables
      real    tmp(lpert,lpert)
      integer ntot,i,j
      character*1 str
      character*3 oname
      character*20 fmtr, fmti
      ! timing
      real ltim

      ! function
      real dnekclock
!----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

!-----------------------------------------------------------------------
!     Compute the eigenspectrum of the reduced operator, sort the
!     eigenvalues in decreasing order
!
      ! LAPACK
      ! compute eigenvalues of Lsym = (Lr+Lr^T)/2
      call copy(tmp,Lr,lpert*lpert)                   ! Save Lr
      do i=1,npert
        do j=1,npert
          Lr(i,j) = 0.5*(tmp(i,j)+tmp(j,i))           ! Compute Lsym
        enddo
      enddo
      call eig_wrapper(npert,'r')                     ! Compute lambdas
      call sorteigs('Ls ')
! Output
      if (mod(istep,otd_prntsp).eq.0) then
        if (nid.eq.0) then
          write(fmtr,'("(", I0, "(E15.7,1X))")') npert
          write(6,100,ADVANCE='NO') istep, time, 'Ls | Re'
          write(6,fmtr) (EIGR(i), i=1,npert)
          write(6,100,ADVANCE='NO') istep, time, 'Ls | Im'
          write(6,fmtr) (EIGI(i), i=1,npert)
        endif
      endif
      
      ! compute eigenvalues of Lr
      call copy(Lr,tmp,lpert*lpert)                   ! Restore Lr
      call eig_wrapper(npert,'r')                     ! Compute lambdas
      call sorteigs('Lr ')
      call copy(Lr,tmp,lpert*lpert)                   ! Restore Lr

! Output
      if (mod(istep,otd_prntsp).eq.0) then
        if (nid.eq.0) then
          write(6,100,ADVANCE='NO') istep, time, 'Lr | Re '
          write(6,fmtr) (EIGR(i), i=1,npert)
          write(6,100,ADVANCE='NO') istep, time, 'Lr | Im '
          write(6,fmtr) (EIGI(i), i=1,npert)
          ! print order for reference
          write(fmti,'("(", I0, "(I4,1X))")') npert
          write(6,100,ADVANCE='NO') istep, time, 's-idx   '
          write(6,fmti) (idx(i), i=1,npert)
          ! print out non-zero elements of rotated Lr
          write(fmtr,'("(", I0, "(E15.7,1X))")') npert*(npert+1)/2
          write(6,100,ADVANCE='NO') istep, time, 'Lrmat   '
          write(6,fmtr) ( ( Lr(i,j) , j=i,npert ), i=1,npert )
        endif
      endif
  100 format('  [OTD] ',I7,1X,E14.7,1X,A8)

!-----------------------------------------------------------------------
!     Project the perturbation velocity field (OTD basis) onto the 
!     eigendirections of the reduced operator to obtain the most
!     unstable directions
!
      ntot = lx1*ly1*lz1*lelv
      call mxm(vxp,ntot,EVRr,lpert,OTDmrx,npert)
      call mxm(vxp,ntot,EVRi,lpert,OTDmix,npert)
      call mxm(vyp,ntot,EVRr,lpert,OTDmry,npert)
      call mxm(vyp,ntot,EVRi,lpert,OTDmiy,npert)
      if (if3d) then
        call mxm(vzp,ntot,EVRr,lpert,OTDmrz,npert)
        call mxm(vzp,ntot,EVRi,lpert,OTDmiz,npert)
      endif

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tget_id,1,ltim)

      return
      end subroutine proj_OTDmodes

!=======================================================================
!> @brief Output the projection of the OTD basis on the eigendirection
!     as fields
!
      subroutine outpost_projmodes()
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'OTD'

      ! local variables
      integer ipert
      character*2 str
      character*3 oname
      ! timing
      real ltim

      ! function
      real dnekclock
!-----------------------------------------------------------------------
      if (mod(istep,otd_prntsp).ne.0) then
        call proj_OTDmodes()
      endif

      ! timing
      ltim = dnekclock()

      do ipert=1,npert
        if (lpert .ge. 10) then
          write(str,'(I2.2)') ipert
          oname = 'o'//trim(str)
        else
          write(str,'(I1)') ipert
          oname = 'ot'//trim(str)
        endif
        call outpost(OTDmrx(1,ipert),OTDmry(1,ipert),OTDmrz(1,ipert)
     $ ,             prp(1,ipert),t,oname)
      enddo

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tIO_id,1,ltim)

      return
      end subroutine outpost_projmodes

!=======================================================================
!> @brief Output the OTD basis directly to restart.
!     We could alternatively reconstruct the OTD basis from the modes
!     but for this we would need both real and imaginary part. Since we
!     currently only outpost the real part, it's cheaper to just outpost
!     the OTD basis directly when we also outpost the baseflow.     
!
      subroutine outpost_OTDbasis()
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'OTD'

      ! local variables
      integer ipert
      character*2 str
      character*3 oname
      ! timing
      real ltim

      ! function
      real dnekclock
!-----------------------------------------------------------------------
      ! orthonormalize
      gsstep_override = .true.
      call OTD_ON()
      ! timing
      ltim = dnekclock()

      do ipert=1,npert
        write(str,'(I2.2)') ipert
        oname = 'r'//trim(str)
        call outpost(vxp(1,ipert),vyp(1,ipert),vzp(1,ipert)
     $ ,             prp(1,ipert),t,oname)
      enddo

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tIO_id,1,ltim)

      return
      end subroutine outpost_OTDbasis

!=======================================================================
!> @brief Initialise all perturbation fields to random noise
!
      subroutine white_noise_IC()
      implicit none

      include 'SIZE'
      include 'INPUT'           ! if3d
      include 'FRAMELP'
      include 'NEKUSE'
      include 'GEOM'            ! [xyz]m1
      include 'PARALLEL'        ! lglel
      include 'OTD'

      ! local variables
      integer ix,iy,iz,ie,ieg,i,ijke
      real    mth_rand,xl(LDIM),fcoeff(3)
!-----------------------------------------------------------------------
      call mntr_log(otd_id,lp_inf,'INIT - creating white noise IC')
      do i=1,npert
        do ix=1,lx1
        do iy=1,ly1
        do iz=1,lz1
          do ie=1,lelv
            xl(1) = XM1(ix,iy,iz,ie)
            xl(2) = YM1(ix,iy,iz,ie)
            if (if3d) then
              xl(NDIM) = ZM1(ix,iy,iz,ie)
            endif
            ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(ie-1)))
            ieg=lglel(ie)
            fcoeff(1)= sin(real(i))**2* 3.0e4
            fcoeff(2)= sin(real(i))**2*(-1.5e3)
            fcoeff(3)= sin(real(i))**2* 0.5e5
            vxpic(ijke,i)=mth_rand(ix,iy,iz,ieg,xl,fcoeff)
            fcoeff(1)= sin(real(i))**2* 2.3e4
            fcoeff(2)= sin(real(i))**2* 2.3e3
            fcoeff(3)= sin(real(i))**2*(-2.0e5)
            vypic(ijke,i)=mth_rand(ix,iy,iz,ieg,xl,fcoeff)
            if (if3d) then
              fcoeff(1)= sin(real(i))**2*2.e4
              fcoeff(2)= sin(real(i))**2*1.e3
              fcoeff(3)= sin(real(i))**2*1.e5
              vzpic(ijke,i)=mth_rand(ix,iy,iz,ieg,xl,fcoeff)
            endif
          enddo
        enddo
        enddo
        enddo
      enddo
      call mntr_log(otd_id,lp_inf,'INIT - white noise :: done')

      return
      end subroutine white_noise_IC

!=======================================================================
!> @brief Read OTD IC fields and run setics 
!
      subroutine read_OTDIC()
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'INPUT'
      include 'SOLN'
      include 'OTD'
      include 'TSTEP'

      ! local variables
      logical       exist_IC
      integer       i
      character*2   istr
      character*132 ifile
!-----------------------------------------------------------------------

      call mntr_log(otd_id,lp_inf,'INIT - read ICs')
      do i=1,otd_nusric
        write(istr,'(I0.2)') i
        ifile='OTDIC_'//trim(istr)//'.fld'
        inquire (file=ifile,exist=exist_IC)
        if (exist_IC) then
          call mntr_log(OTD_id,lp_inf,
     $             'INIT - Reading IC file '//trim(ifile))
          call load_fld(ifile)
          ! Copy the initial conditions into the fields v[xyz]pic
          call opcopy(vxpic(1,i),vypic(1,i),vzpic(1,i),vx,vy,vz)
          OTDrsttime = time
          call mntr_logr(otd_id,lp_inf,
     $                   'IC file '//trim(ifile)//' time:',time)
        else
          call mntr_abort(otd_id,'Cannot open '//trim(ifile)//' !')
        endif
      enddo
      call mntr_log(otd_id,lp_inf,'INIT - read ICs :: done')

      return
      end subroutine read_OTDIC

!=======================================================================
!> @brief Create forcing for OTD evolution equation 
!
      subroutine gen_OTD_forcing()
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'OTD'
      include 'INPUT'   ! if3d
      include 'TSTEP'

      ! local variable list
      integer ntot

      ! timing
      real ltim

      ! function 
      real dnekclock

!#define DEBUG
#ifdef DEBUG
      character*1 str
      character*3 oname
      integer i
      real glmax 
      real gl(3)
      real dudtx (lx1*ly1*lz1*lelv,lpert)
     $ ,   dudty (lx1*ly1*lz1*lelv,lpert)
     $ ,   dudtz (lx1*ly1*lz1*lelv,lpert)
#endif
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ntot = lx1*ly1*lz1*lelv

      call mxm(VXP,ntot,Lr,lpert,OTDfx,npert)
      call mxm(VYP,ntot,Lr,lpert,OTDfy,npert)
      call mxm(VZP,ntot,Lr,lpert,OTDfz,npert)

#ifdef DEBUG
      if (mod(istep+1,otd_iostep).eq.0) then
        do i=1,npert
          write(str,'(I1)') i
          oname = 'p'//trim(str)//'f'
          call outpost(OTDfx(1,i),OTDfy(1,i),OTDfz(1,i),pr,t,oname)
          call opcopy(dudtx(1,i),dudty(1,i),dudtz(1,i)
     $ ,              Lux(1,i),Luy(1,i),Luz(1,i))
          call sub2(dudtx(1,i),OTDfx(1,i),ntot)
          call sub2(dudty(1,i),OTDfy(1,i),ntot)
          call sub2(dudtz(1,i),OTDfz(1,i),ntot)
          oname = 'p'//trim(str)//'u'
          call outpost(dudtx(1,i),dudty(1,i),dudtz(1,i),prp,t,oname)
          oname = 'p'//trim(str)//'l'
          call outpost(Lux(1,i),Luy(1,i),Luz(1,i),prp,t,oname)
        enddo
      endif
#endif
#undef DEBUG
      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tgen_id,0,ltim)

      return 
      end subroutine gen_OTD_forcing

!=======================================================================
!> @brief Set forcing for OTD evolution equation 
!     This function is called in userf for each GLL point. Therefore we
!     need a manual switch for when to start setting the forcing.
!
      subroutine set_OTD_forcing(FFX,FFY,FFZ,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'
      include 'SOLN'           ! jp
      include 'TSTEP'          ! istep
      include 'PARALLEL'       ! gllel
      include 'OTD'

      ! argument list
      integer ix,iy,iz,ieg
     
      ! output 
      real    ffx,ffy,ffz

      ! local variable list
      integer ijke,e
!-----------------------------------------------------------------------

      e = gllel(ieg)

      ijke = ix + lx1*((iy-1) + ly1*((iz-1) + lz1*(e-1)))
      if (jp.ne.0) then
!       only for the perturbations
        FFX = FFX - OTDfx(ijke,jp)
        FFY = FFY - OTDfy(ijke,jp)
        FFZ = FFZ - OTDfz(ijke,jp)
      endif

      return 
      end subroutine set_OTD_forcing
 
!=======================================================================
!> @brief Orthonormalize OTD basis 
!
      subroutine OTD_ON()
      implicit none

      include 'SIZE'
      include 'OTD'
     
      include 'SOLN'    ! v[xyz]p 
      include 'TSTEP'   ! istep

      ! local variables
      real N,O
      logical runON

      ! timing
      real ltim

      ! function
      real dnekclock
!----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()
!
      call compute_NO(N,O,'pre   ',.true.)
!
      runON = .false.
      if (otd_gsstep.ne.0) then                 ! gsstep=0 => no GS
        if (mod(istep,otd_gsstep).eq.0) then
          runON = .true.
        endif
      endif
!     override
      if (gsstep_override) runON = .true.       ! override when needed

!     runON
      if (runON) then
!        call CGS()        ! Classical Gram-Schmidt
        call MGS()        ! Modified Gram-Schmidt
      endif
!
      call compute_NO(N,O,'post  ',.false.)
!
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tGS_id,1,ltim)
 
      return     
      end subroutine OTD_ON

!=======================================================================
!> @brief Reset FTLE computation 
!
      subroutine reset_FTLE()
      implicit none

      include 'SIZE'
      include 'OTD'
      
      call rzero(FTLEv,lpert)
      ! Phi
      !call ident(Phi,lpert)
      ! Blanchard
      call rzero(trapz,lpert)

      end subroutine reset_FTLE

!=======================================================================
!> @brief Compute FTLEs (Based on Babaee et al., 2017) 
!
!     1. Advect fundamental solution matrix Phi
!
!       dPhi/dt = Lr    , with Phi(t=t0) = I(rxr)
!
!     2. Compute FTLEs
!
!       FTLE(i) = 1/T * log(svd(Phi(t)))
!
      subroutine compute_FTLE_Phi(Lrp,Lrc,deltat,ftledt)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'

      ! argument list
      real deltat                       ! time interval for FTLE computation 
      real ftledt                       ! dt for Phi advection
      real Lrp   (lpert,lpert)          ! Lr matrix from previous step
      real Lrc   (lpert,lpert)          ! Lr matrix for current step

      ! local variables
      real lhs   (lpert,lpert)          ! lhs of Phi advection equation
     $ ,   invlhs(lpert,lpert)          ! inverse of lhs
     $ ,   tmp   (lpert,lpert)           
     $ ,   rhs   (lpert,lpert)          ! rhs of Phi advection equation
      real fact
      integer i

      ! Advect the fundamental solution matrix (using the implicit CN scheme)
      fact = 0.5*ftledt
      ! build LHS
      call ident(lhs,lpert)
      call add2s2(lhs,Lrc,-fact,lpert*lpert)
      ! build RHS
      call ident(tmp,npert)
      call add2s2(tmp,Lrp, fact,lpert*lpert)
      call mxm(tmp,lpert,Phi,lpert,rhs,lpert)
      ! invert LHS and solve system
      call invmt(lhs,invlhs,tmp,npert)
      call mxm(invlhs,lpert,rhs,lpert,Phi,lpert)

      ! compute SVD
      call copy(VMATX,Phi,lpert*lpert)
      call svd_wrapper(lpert,lpert,'A')

      ! compute FTLEs
      do i=1,npert
        FTLEv(i) = log(OSIGMA(i))/deltat
      enddo

      end subroutine compute_FTLE_Phi

!=======================================================================
!> @brief Compute FTLEs (Based on Blanchard & Sapsis, 2019) 
!
!       FTLE(i) = 1/T * ( int_(t_0)^t <Lu_i,u_i> d tau )
!
      subroutine compute_FTLE_Blanchard(Lrp,Lrc,deltat,ftledt)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'

      ! argument list
      real deltat                       ! time interval for FTLE computation 
      real ftledt                       ! dt for Phi advection
      real Lrp  (lpert,lpert)           ! Lr matrix from previous step
      real Lrc  (lpert,lpert)           ! Lr matrix for current step

      ! local variables
      integer i

      ! compute FTLEs
      do i=1,npert
        ! Trapezoid rule, this could be done better using the AB scheme
        trapz(i) = trapz(i) + 0.5*ftledt*(Lrp(i,i) + Lrc(i,i))
        FTLEv(i) = trapz(i)/deltat
      enddo

      end subroutine compute_FTLE_Blanchard

!=======================================================================
!> @brief Get FTLEs
!      
      subroutine get_FTLE()

      include 'SIZE'
      include 'TSTEP'
      include 'OTD'             ! Lr, FTLEv

      ! local variables
      real pfrac                ! current fraction of the FTLE comp. period
      real ftledt               ! dt for FTLE computation
      real Lrp(lpert,lpert)     ! Lr from previous step
      real Lrc(lpert,lpert)     ! (linear approx.) of Lr at end of int. interval
      real fact
      real t0
      integer icalld
      ! timing
      real ltim
      ! function
      real dnekclock
      ! save variables
      save Lrp
      data t0 /0.0/
      save t0
      data icalld /0/
      save icalld

!----------------------------------------------------------------------

      ! timing
      ltim = dnekclock()

      ! determine FTLE horizon
      if (otd_FTLEpd.eq.0.0) then
        if (icalld.eq.0) then
          t0 = time
          icalld = 1
        endif
        period = time-t0
        pfrac = time-t0
      else
        period = otd_FTLEpd
        pfrac = mod(time,period)
      endif
      ! initialisation
      if (istep.eq.startOTD) then
        call copy(Lrp,Lr,lpert*lpert)
      endif

      if (pfrac.lt.dt) then
        ftledt = dt-pfrac
        ! compute approximation of Lrc at end of period (linear interp.)
        call copy(Lrc,Lrp,lpert*lpert)
        fact = ftledt/dt
        call add2s2(Lrc,Lrp,-fact,lpert*lpert)
        call add2s2(Lrc,Lr , fact,lpert*lpert)
        !call compute_FTLE_Phi(Lrp,Lrc,period,ftledt)
        call compute_FTLE_Blanchard(Lrp,Lrc,period,ftledt)
        pcount = pcount + 1
        if (nid.eq.0) then
          write(6,100,ADVANCE='NO') pcount,istep,time
          do i=1,npert
            write(6,102,ADVANCE='NO') FTLEv(i)
          enddo
          write(6,*)
        endif
        call reset_FTLE
        call copy(Lrp,Lrc,lpert*lpert)
        call copy(Lrc,Lr ,lpert*lpert)
        !call compute_FTLE_Phi(Lrp,Lrc,period,pfrac)
        call compute_FTLE_Blanchard(Lrp,Lrc,period,pfrac)
      else
        call copy(Lrc,Lr,lpert*lpert)
        !call compute_FTLE_Phi(Lrp,Lrc,pfrac,dt)
        call compute_FTLE_Blanchard(Lrp,Lrc,pfrac,dt)
      endif
      ! update Lrp
      call copy(Lrp,Lr,lpert*lpert)
 
      if (nid.eq.0.and.mod(istep,otd_prntsp).eq.0) then
        write(6,101,ADVANCE='NO') istep,time,pfrac
        do i=1,npert
          write(6,102,ADVANCE='NO') FTLEv(i)
        enddo
        write(6,*)
      endif
  100 format(' [OTD] FTLE PRD',1X,I5,1X,I7,' t=',1X,E15.7)
  101 format(' [OTD] FTLE (t)',1X,I7,' t=',2(1X,E15.7))
  102 format(1X,E15.7)

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(otd_tFTLE_id,1,ltim)

      end subroutine get_FTLE
