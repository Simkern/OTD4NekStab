!> @file spng_otd.f
!! @brief Sponge/fringe for localized OTD
!=======================================================================
!> @brief Register spng_otd module
!! @note This routine should be called in frame_usr_register
      subroutine spng_otd_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'SPNGOTD'

      ! local variables
      integer lpmid
      real ltim

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,spng_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(spng_name)//'] already registered')
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
      call mntr_mod_reg(spng_id,lpmid,spng_name,
     $          'Sponge/fringe for rectangular domain')

      ! register timer
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
      call mntr_tmr_reg(spng_tmr_id,lpmid,spng_id,
     $     'SPNG_INI','Sponge calculation initialisation time',.false.)

      ! register and set active section
      call rprm_sec_reg(spng_sec_id,spng_id,'_'//adjustl(spng_name),
     $     'Runtime paramere section for sponge_box module')
      call rprm_sec_set_act(.true.,spng_sec_id)

      ! register parameters
      call rprm_rp_reg(spng_str_id,spng_sec_id,'STRENGTH',
     $     'Sponge strength',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spng_ws_id,spng_sec_id,'WIDTHS',
     $     'Width of whole sponge region',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spng_wt_id,spng_sec_id,'WIDTHT',
     $     'Width of transition region',rpar_real,0,0.0,.false.,' ')
      
      call rprm_rp_reg(spng_pt_id(1),spng_sec_id,'PX',
     $     'Corner point: coordinate X',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spng_pt_id(2),spng_sec_id,'PY',
     $     'Corner point: coordinate Y',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spng_a_id(1),spng_sec_id,'ATANG',
     $     'Slope of tangential wall',rpar_real,0,0.0,.false.,' ')

      call rprm_rp_reg(spng_a_id(2),spng_sec_id,'ANORM',
     $     'Slope of normal wall',rpar_real,0,0.0,.false.,' ')

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(spng_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise spng_otd module
!! @note This routine should be called in frame_usr_init
!! @remark This routine uses global scratch space \a SCRUZ
      subroutine spng_otd_init()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'FRAMELP'
      include 'SPNGOTD'

      ! local variables
      integer ierr, nhour, nmin
      integer itmp
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      integer ntot, il

      real xx,yy,ytmp,n1,n2,b,d

      ! functions
      real dnekclock, mth_stepf

#define DEBUG
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (spng_ifinit) then
         call mntr_warn(spng_id,
     $        'module ['//trim(spng_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()

      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_str_id,rpar_real)
      spng_str = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_ws_id,rpar_real)
      spng_ws = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_wt_id,rpar_real)
      spng_wt = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_pt_id(1),rpar_real)
      spng_pt(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_pt_id(2),rpar_real)
      spng_pt(2) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_a_id(1),rpar_real)
      spng_a(1) = rtmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,spng_a_id(2),rpar_real)
      spng_a(2) = rtmp

      ! initialise sponge variables
      ntot = lx1*ly1*lz1*lelv
      call rzero(spng_fun,ntot)

      if(spng_str.gt.0.0) then
         call mntr_log(spng_id,lp_inf,"Sponge turned on")

         ! save reference field (zero)
         call rzero(spng_vr(1,1),ntot)
         call rzero(spng_vr(1,2),ntot)
         if (IF3D) call rzero(spng_vr(1,NDIM),ntot)

         if (spng_ws.gt.0.0) then
           if (spng_wt.gt.spng_ws) then
             call mntr_abort(spng_id,"Wrong sponge parameters")
           endif

           ! get SPNG_FUN

           ! normal vector from top bdry
           rtmp = sqrt(1 + 1/spng_a(1)**2)
           n1 = 1/rtmp
           n2 = -1/spng_a(1)/rtmp
           ! b = (y + d*n2) - a*(x + d*n1) 
           b = (spng_pt(2) - n2*spng_ws) 
     $         - spng_a(1)*(spng_pt(1) - n1*spng_ws)

           do il=1,ntot
             xx = xm1(il,1,1,1)
             yy = ym1(il,1,1,1)
             ytmp = spng_a(1)*xx + b
             d = abs(spng_a(1)*xx - yy + b)/sqrt(1+spng_a(1)**2)
             if (yy.lt.ytmp) then
               rtmp = 0.0
             elseif (d.lt.spng_wt) then
               rtmp = mth_stepf(d/spng_wt)
             else 
               rtmp = 1.0
             endif
             spng_fun(il) = max(spng_fun(il),rtmp)
           enddo        ! jl=1,ntot

           ! normal vector from outflow bdry
           rtmp = sqrt(1 + 1/spng_a(2)**2)
             write(6,*) 'CHK',n1,n2,spng_a(2)
           n1 = 1/rtmp
           n2 = -1/spng_a(2)/rtmp
           ! b = (y + d*n2) - a*(x + d*n1) 
           b = (spng_pt(2) - n2*spng_ws) 
     $         - spng_a(2)*(spng_pt(1) - n1*spng_ws)

           do il=1,ntot
             xx = xm1(il,1,1,1)
             yy = ym1(il,1,1,1)
             ytmp = spng_a(2)*xx + b
             d = abs(spng_a(2)*xx - yy + b)/sqrt(1+spng_a(2)**2)
             if (ytmp.lt.yy) then
               rtmp = 0.0
             elseif (d.lt.spng_wt) then
               rtmp = mth_stepf(d/spng_wt)
             else 
               rtmp = 1.0
             endif
             spng_fun(il) = spng_str*max(spng_fun(il),rtmp)
           enddo        ! jl=1,ntot
         endif            ! spng_ws.gt.0.0
      endif               ! spng_str.gt.0.0

#ifdef DEBUG
      ! for debugging
      ltmp = ifto
      ifto = .TRUE.
      call outpost2(spng_vr(1,1),spng_vr(1,2),spng_vr(1,NDIM),spng_fun,
     $              spng_fun,1,'spg')
      ifto = ltmp
#endif

      ! is everything initialised
      spng_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(spng_tmr_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup sponge_box
!! @return spng_is_initialised
      logical function spng_is_initialised()
      implicit none

      include 'SIZE'
      include 'SPNGOTD'
!-----------------------------------------------------------------------
      spng_is_initialised = spng_ifinit

      return
      end function
!=======================================================================
!> @brief Get sponge forcing
!! @ingroup sponge_box
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    ieg             global element number
      subroutine spng_forcing(ffx,ffy,ffz,ix,iy,iz,ieg)
      implicit none

      include 'SIZE'            !
      include 'INPUT'           ! IF3D
      include 'PARALLEL'        ! GLLEL
      include 'SOLN'            ! JP
      include 'SPNGOTD'

      ! argument list
      real ffx, ffy, ffz
      integer ix,iy,iz,ieg

      ! local variables
      integer iel, ip
!-----------------------------------------------------------------------
      iel=GLLEL(ieg)
      if (SPNG_STR.gt.0.0) then
         ip=ix+NX1*(iy-1+NY1*(iz-1+NZ1*(iel-1)))

         if (JP.eq.0) then
            return ! dns: do nothing
         else
            ! perturbation
            ffx = ffx + SPNG_FUN(ip)*(spng_vr(ip,1) - VXP(ip,JP))
            ffy = ffy + SPNG_FUN(ip)*(spng_vr(ip,2) - VYP(ip,JP))
            if (IF3D) then
               ffz = ffz + SPNG_FUN(ip)*(spng_vr(ip,ndim) - VZP(ip,JP))
            endif
         endif

      endif

      return
      end subroutine
!=======================================================================
