!!> \file modsubgrid.f90
!!!  Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module modsubgrid
use modsubgriddata
implicit none
save
  public :: subgrid, initsubgrid, exitsubgrid, subgridnamelist

contains
  subroutine initsubgrid
    use modglobal, only : ih,i1,jh,j1,k1,delta,zf,fkar,pi
    use modmpi, only : myid

    implicit none

    integer   :: k

    real :: ceps, ch
    real :: mlen

    call subgridnamelist

    allocate(ekm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ekh(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(zlt(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbdiss(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbshr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbbuo(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(csz(k1))

    ! Initialize variables to avoid problems when not using subgrid scheme JvdD
    ekm=0.; ekh=0.; zlt=0.; sbdiss=0.; sbshr=0.; sbbuo=0.; csz=0.

    cm = cf / (2. * pi) * (1.5*alpha_kolm)**(-1.5)

!     ch   = 2. * alpha_kolm / beta_kolm
    ch   = 1.0/Prandtl
    ch2  = ch-ch1

    ceps = 2. * pi / cf * (1.5*alpha_kolm)**(-1.5)
    ce1  = (cn**2)* (cm/Rigc - ch1*cm )
    ce2  = ceps - ce1

    if(cs < 0.) then
      csz(:)  = (cm**3/ceps)**0.25   !< Smagorinsky constant
    else
      csz(:)  = cs
    end if

    if(lmason) then
      do k = 1,k1
        mlen   = (1. / (csz(k) * delta(k))**nmason + 1. / (fkar * zf(k))**nmason)**(-1./nmason)
        csz(k) = mlen / delta(k)
      end do
    end if

    if (myid==0) then
      write (6,*) 'cf    = ',cf
      write (6,*) 'cm    = ',cm
      write (6,*) 'ch    = ',ch
      write (6,*) 'ch1   = ',ch1
      write (6,*) 'ch2   = ',ch2
      write (6,*) 'ceps  = ',ceps
      write (6,*) 'ceps1 = ',ce1
      write (6,*) 'ceps2 = ',ce2
      write (6,*) 'cs    = ',cs
      write (6,*) 'Rigc  = ',Rigc
    endif

  end subroutine initsubgrid

  subroutine subgridnamelist
    use modglobal, only : ifnamopt,fname_options,checknamelisterror
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical

    implicit none

    integer :: ierr

    namelist/NAMSUBGRID/ &
        ldelta,lmason,cf,cn,Rigc,Prandtl,lsmagorinsky,cs,nmason,ch1

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSUBGRID,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMSUBGRID')
      write(6 ,NAMSUBGRID)
      close(ifnamopt)
    end if

    call MPI_BCAST(ldelta     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lmason     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(nmason     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lsmagorinsky,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(cs         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cf         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cn         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Rigc       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Prandtl    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(sgs_surface_fix ,1,MPI_LOGICAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ch1        ,1,MY_REAL   ,0,comm3d,mpierr)

  end subroutine subgridnamelist

  subroutine subgrid

 ! Diffusion subroutines
! Thijs Heus, Chiel van Heerwaarden, 15 June 2007

    use modglobal, only : nsv, lmoist
    use modfields, only : up,vp,wp,e12p,thl0,thlp,qt0,qtp,sv0,svp
    use modsurfdata,only : thlflux,qtflux,svflux
	use modmpi,     only : myid
    implicit none
    integer n

    call closure
    call diffu(up)
    call diffv(vp)
    call diffw(wp)
    if (.not. lsmagorinsky) call diffe(e12p)
    call diffc(thl0,thlp,thlflux)
    if (lmoist) call diffc( qt0, qtp, qtflux)
	!if(myid==0) write(6,*) 'Scalar before diffc in subgrid',svp(8,9,10,1)
    do n=1,nsv
      call diffc(sv0(:,:,:,n),svp(:,:,:,n),svflux(:,:,n))
    end do
	!if(myid==0) write(6,*) 'Scalar after diffc in subgrid',svp(8,9,10,1)
    if (.not. lsmagorinsky) call sources
  end subroutine

  subroutine exitsubgrid
    implicit none
    deallocate(ekm,ekh,zlt,sbdiss,sbbuo,sbshr,csz)
  end subroutine exitsubgrid

   subroutine closure

!-----------------------------------------------------------------|
!                                                                 |
!*** *closure*  calculates K-coefficients                         |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.   06/01/1995                      |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     All the K-closure factors are calculated.                   |
!                                                                 |
!     ekm(i,j,k) = k sub m : for velocity-closure                 |
!     ekh(i,j,k) = k sub h : for temperture-closure               |
!     ekh(i,j,k) = k sub h = k sub c : for concentration-closure  |
!                                                                 |
!     We will use the next model for these factors:               |
!                                                                 |
!     k sub m = 0.12 * l * sqrt(E)                                |
!                                                                 |
!     k sub h = k sub c = ( 1 + (2*l)/D ) * k sub m               |
!                                                                 |
!           where : l = mixing length  ( in model = z2 )          |
!                   E = subgrid energy                            |
!                   D = grid-size distance                        |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *closure* is called from *program*.                 |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,     only : i1,j1,kmax,k1,ih,jh,i2,j2,delta,ekmin,grav,zf,fkar,deltai, &
                          dxi,dyi,dzf,dzh,imax,jmax
  use modfields,     only : dthvdz,e120,u0,v0,w0,thvf
  use modsurfdata,   only : dudz,dvdz,z0m
  use modmpi,        only : excjs,myidx,myidy
  use modruraldata,  only : applydamping,bc_height
  implicit none

  real    :: strain2,mlen
  integer :: i,j,k,kp,km,jp,jm,kmin

  if(lsmagorinsky) then
    do k = 1,kmax
      mlen        = csz(k) * delta(k)
    end do

    do i = 2,i1
      do j = 2,j1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
        do k=kmin,kmax
          kp=k+1
          km=k-1
          jp=j+1
          jm=j-1

          if(k == kmin) then
            strain2 =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

            strain2 = strain2 + 0.5 * ( &
              ( 0.25*(w0(i+1,j,kp)-w0(i-1,j,kp))*dxi + &
              dudz(i,j)   )**2 )

            strain2 = strain2 + 0.125 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain2 = strain2 + 0.5 * ( &
              ( 0.25*(w0(i,jp,kp)-w0(i,jm,kp))*dyi + &
              dvdz(i,j)   )**2 )

          else

            strain2 =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

            strain2 = strain2 + 0.125 * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  *dxi     + &
              (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    *dxi     + &
              (u0(i,j,k)-u0(i,j,km))     / dzh(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    *dxi     + &
              (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi     + &
              (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )

            strain2 = strain2 + 0.125 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain2 = strain2 + 0.125 * ( &
              ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
              (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     / dzh(k)+ &
              (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))   / dzh(k)+ &
              (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
              (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )
          end if

          ekm(i,j,k)  = mlen ** 2. * sqrt(2. * strain2)
          call applydamping(ekm(i,j,k),i,j,k)   !< MK - apply damping function in modruraldata.f90 when it is used.
          ekh(i,j,k)  = ekm(i,j,k) / Prandtl

          ekm(i,j,k) = max(ekm(i,j,k),ekmin)
          ekh(i,j,k) = max(ekh(i,j,k),ekmin)
        end do
      end do
    end do

  ! do TKE scheme
  else
    do j=2,j1
      do i=2,i1
      kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
      do k=kmin,kmax
          if (ldelta .or. (dthvdz(i,j,k)<=0)) then
            zlt(i,j,k) = delta(k)
            if (lmason) zlt(i,j,k) = (1. / zlt(i,j,k) ** nmason + 1. / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)
            ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
            call applydamping(ekm(i,j,k),i,j,k)   !< MK - apply damping function in modruraldata.f90 when it is used.
            ekh(i,j,k) = (ch1 + ch2) * ekm(i,j,k)

            ekm(i,j,k) = max(ekm(i,j,k),ekmin)
            ekh(i,j,k) = max(ekh(i,j,k),ekmin)
          else
             ! zlt(i,j,k) = min(delta(k),cn*e120(i,j,k)/sqrt(grav/thvf(k)*abs(dthvdz(i,j,k))))
             ! faster calculation: evaluate sqrt only if the second argument is actually smaller
             zlt(i,j,k) = delta(k)
             if ( grav*abs(dthvdz(i,j,k)) * delta(k)**2 > (cn*e120(i,j,k))**2 * thvf(k) ) then
                zlt(i,j,k) = cn*e120(i,j,k)/sqrt(grav/thvf(k)*abs(dthvdz(i,j,k)))
             end if

            if (lmason) zlt(i,j,k) = (1. / zlt(i,j,k) ** nmason + 1. / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)

            ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
            call applydamping(ekm(i,j,k),i,j,k)   !< MK - apply damping function in modruraldata.f90 when it is used.
            ekh(i,j,k) = (ch1 + ch2 * zlt(i,j,k)*deltai(k)) * ekm(i,j,k)

            ekm(i,j,k) = max(ekm(i,j,k),ekmin)
            ekh(i,j,k) = max(ekh(i,j,k),ekmin)
          endif
        end do
      end do
    end do
  end if

!*************************************************************
!     Set cyclic boundary condition for K-closure factors.
!*************************************************************

  call excjs( ekm           , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( ekh           , 2,i1,2,j1,1,k1,ih,jh)

  do j=1,j2
    do i=1,i2
      ekm(i,j,k1)  = ekm(i,j,kmax)
      ekh(i,j,k1)  = ekh(i,j,kmax)
    end do
  end do

  return
  end subroutine closure
  subroutine sources


!-----------------------------------------------------------------|
!                                                                 |
!*** *sources*                                                    |
!      calculates various terms from the subgrid TKE equation     |
!                                                                 |
!     Hans Cuijpers   I.M.A.U.     06/01/1995                     |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Subroutine sources calculates all other terms in the       |
!      subgrid energy equation, except for the diffusion terms.   |
!      These terms are calculated in subroutine diff.             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *sources* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,   only : i1,j1,kmax,delta,dx,dy,dxi,dyi,dzf,dzh,grav,cu,cv,deltai,imax,jmax
  use modfields,   only : u0,v0,w0,e120,e12p,dthvdz,thvf,e12m
  use modsurfdata,  only : dudz,dvdz,ustar,thlflux
  use modsubgriddata, only: sgs_surface_fix
  use modruraldata,   only: bc_height
  use modmpi,         only: myidx,myidy

  implicit none

  real    tdef2, uwflux, vwflux, local_dudz, local_dvdz, local_dthvdz, horv
  integer i,j,k,jm,jp,km,kp,kmin


  do j=2,j1
  do i=2,i1
  kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
  do k=(kmin+1),kmax
    kp=k+1
    km=k-1
    jp=j+1
    jm=j-1

    tdef2 = 2. * ( &
             ((u0(i+1,j,k)-u0(i,j,k))   /dx         )**2    + &
             ((v0(i,jp,k)-v0(i,j,k))    /dy         )**2    + &
             ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

    tdef2 = tdef2 + 0.25 * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  / dx     + &
               (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    / dx     + &
               (u0(i,j,k)-u0(i,j,km))     / dzh(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    / dx     + &
               (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  / dx     + &
               (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )

    tdef2 = tdef2 + 0.25 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     / dy     + &
               (v0(i,jp,k)-v0(i-1,jp,k))  / dx        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     / dy     + &
               (v0(i,j,k)-v0(i-1,j,k))    / dx        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) / dy     + &
               (v0(i+1,j,k)-v0(i,j,k))    / dx        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) / dy     + &
               (v0(i+1,jp,k)-v0(i,jp,k))  / dx        )**2    )

    tdef2 = tdef2 + 0.25 * ( &
              ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
               (w0(i,j,kp)-w0(i,jm,kp))   / dy        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     / dzh(k)+ &
               (w0(i,j,k)-w0(i,jm,k))     / dy        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))   / dzh(k)+ &
               (w0(i,jp,k)-w0(i,j,k))     / dy        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
               (w0(i,jp,kp)-w0(i,j,kp))   / dy        )**2    )


!    sbshr(i,j,k)  = ekm(i,j,k)*tdef2/ ( 2*e120(i,j,k))
!    sbbuo(i,j,k)  = -ekh(i,j,k)*grav/thvf(k)*dthvdz(i,j,k)/ ( 2*e120(i,j,k))
!    sbdiss(i,j,k) = - (ce1 + ce2*zlt(i,j,k)/delta(k)) * e120(i,j,k)**2 /(2.*zlt(i,j,k))

!     e12p(2:i1,2:j1,1) = e12p(2:i1,2:j1,1)+ &
!            sbshr(2:i1,2:j1,1)+sbbuo(2:i1,2:j1,1)+sbdiss(2:i1,2:j1,1)
    e12p(i,j,k) = e12p(i,j,k) &
                + (ekm(i,j,k)*tdef2 - ekh(i,j,k)*grav/thvf(k)*dthvdz(i,j,k) ) / (2*e120(i,j,k)) &  !  sbshr and sbbuo
                - (ce1 + ce2*zlt(i,j,k)*deltai(k)) * e120(i,j,k)**2 /(2.*zlt(i,j,k))               !  sbdiss

  end do
  end do
  end do
!     ----------------------------------------------end i,j,k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=kmin
!     --------------------------------------------


  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1
    kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1



! **  Calculate "shear" production term: tdef2  ****************

    tdef2 =  2. * ( &
            ((u0(i+1,j,kmin)-u0(i,j,kmin))*dxi)**2 &
          + ((v0(i,jp,kmin)-v0(i,j,kmin))*dyi)**2 &
          + ((w0(i,j,kmin+1)-w0(i,j,kmin))/dzf(kmin))**2   )

    if (sgs_surface_fix) then
          ! Use known surface flux and exchange coefficient to derive
          ! consistent gradient (such that correct flux will occur in
          ! shear production term)
          ! Make sure that no division by zero occurs in determination of the
          ! directional component; ekm should already be >= ekmin
          ! Replace the dudz by surface flux -uw / ekm
          horv = max(sqrt((u0(i,j,kmin)+cu)**2+(v0(i,j,kmin)+cv)**2),  0.01)
          uwflux = -ustar(i,j)*ustar(i,j)* ((u0(i,j,kmin)+cu)/horv)
          local_dudz = -uwflux / ekm(i,j,kmin)
          tdef2 = tdef2 + ( 0.25*(w0(i+1,j,kmin+1)-w0(i-1,j,kmin+1))*dxi + &
               local_dudz )**2
    else
          tdef2 = tdef2 + ( 0.25*(w0(i+1,j,kmin+1)-w0(i-1,j,kmin+1))*dxi + &
                                  dudz(i,j)   )**2
    endif

    tdef2 = tdef2 +   0.25 *( &
          ((u0(i,jp,kmin)-u0(i,j,kmin))*dyi+(v0(i,jp,kmin)-v0(i-1,jp,kmin))*dxi)**2 &
         +((u0(i,j,kmin)-u0(i,jm,kmin))*dyi+(v0(i,j,kmin)-v0(i-1,j,kmin))*dxi)**2 &
         +((u0(i+1,j,kmin)-u0(i+1,jm,kmin))*dyi+(v0(i+1,j,kmin)-v0(i,j,kmin))*dxi)**2 &
         +((u0(i+1,jp,kmin)-u0(i+1,j,kmin))*dyi+ &
                                 (v0(i+1,jp,kmin)-v0(i,jp,kmin))*dxi)**2   )

    if (sgs_surface_fix) then
          ! Use known surface flux and exchange coefficient to derive
          ! consistent gradient (such that correct flux will occur in
          ! shear production term)
          ! Make sure that no division by zero occurs in determination of the
          ! directional component; ekm should already be >= ekmin
          ! Replace the dvdz by surface flux -vw / ekm
          horv = max(sqrt((u0(i,j,kmin)+cu)**2+(v0(i,j,kmin)+cv)**2),  0.01)
          vwflux = -ustar(i,j)*ustar(i,j)* ((v0(i,j,kmin)+cv)/horv)
          local_dvdz = -vwflux / ekm(i,j,kmin)
          tdef2 = tdef2 + ( 0.25*(w0(i,jp,kmin+1)-w0(i,jm,kmin+1))*dyi + &
                        local_dvdz  )**2
    else
         tdef2 = tdef2 + ( 0.25*(w0(i,jp,kmin+1)-w0(i,jm,kmin+1))*dyi + &
                                dvdz(i,j)   )**2
    endif

! **  Include shear and buoyancy production terms and dissipation **

    sbshr(i,j,kmin)  = ekm(i,j,kmin)*tdef2/ ( 2*e120(i,j,kmin))
    if (sgs_surface_fix) then
          ! Replace the -ekh *  dthvdz by the surface flux of thv
          ! (but we only have the thlflux , which seems at the surface to be
          ! equivalent
          local_dthvdz = -thlflux(i,j)/ekh(i,j,kmin)
          sbbuo(i,j,kmin)  = -ekh(i,j,kmin)*grav/thvf(kmin)*local_dthvdz/ ( 2*e120(i,j,kmin))
    else
          sbbuo(i,j,kmin)  = -ekh(i,j,kmin)*grav/thvf(kmin)*dthvdz(i,j,kmin)/ ( 2*e120(i,j,kmin))
    endif
    sbdiss(i,j,kmin) = - (ce1 + ce2*zlt(i,j,kmin)*deltai(kmin)) * e120(i,j,kmin)**2 /(2.*zlt(i,j,kmin))

  e12p(i,j,kmin) = e12p(i,j,kmin) + &
            sbshr(i,j,kmin)+sbbuo(i,j,kmin)+sbdiss(i,j,kmin)

  end do
  end do

!  e12p(2:i1,2:j1,1:kmax) = e12p(2:i1,2:j1,1:kmax)+ &
!            sbshr(2:i1,2:j1,1:kmax)+sbbuo(2:i1,2:j1,1:kmax)+sbdiss(2:i1,2:j1,1:kmax)
!  e12p(2:i1,2:j1,kmin) = e12p(2:i1,2:j1,kmin) + &
!            sbshr(2:i1,2:j1,kmin)+sbbuo(2:i1,2:j1,kmin)+sbdiss(2:i1,2:j1,kmin)

  return
  end subroutine sources

  subroutine diffc (putin,putout,flux)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dy2i,dzh,imax,jmax
    use modfields, only : rhobf,rhobh
    use modmpi,    only : myidx,myidy,myid
    use modruraldata, only : bc_height
    implicit none

    real, intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in)    :: flux (i2,j2)

    integer i,j,k,jm,jp,km,kp,kmin

    do j=2,j1
      jp=j+1
      jm=j-1

      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
        do k=(kmin+1),kmax
          kp=k+1
          km=k-1


          putout(i,j,k) = putout(i,j,k) &
                    +  0.5 * ( &
                  ( (ekh(i+1,j,k)+ekh(i,j,k))*(putin(i+1,j,k)-putin(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k)))*dx2i &
                    + &
                  ( (ekh(i,jp,k)+ekh(i,j,k)) *(putin(i,jp,k)-putin(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i,jm,k)) *(putin(i,j,k)-putin(i,jm,k)) )*dy2i &
                  + &
                  ( rhobh(kp)/rhobf(k) * (dzf(kp)*ekh(i,j,k) + dzf(k)*ekh(i,j,kp)) &
                    *  (putin(i,j,kp)-putin(i,j,k)) / dzh(kp)**2 &
                    - &
                    rhobh(k)/rhobf(k) * (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km)) &
                    *  (putin(i,j,k)-putin(i,j,km)) / dzh(k)**2           )/dzf(k) &
                            )

        end do
      end do
    end do

    do j=2,j1
      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
		!if(myid==0 .and.i==7 .and. j==3) write(6,*) 'kmin = ',kmin
        putout(i,j,kmin) = putout(i,j,kmin) &
                  + 0.5 * ( &
                ( (ekh(i+1,j,kmin)+ekh(i,j,kmin))*(putin(i+1,j,kmin)-putin(i,j,kmin)) &
                  -(ekh(i,j,kmin)+ekh(i-1,j,kmin))*(putin(i,j,kmin)-putin(i-1,j,kmin)) )*dx2i &
                  + &
                ( (ekh(i,j+1,kmin)+ekh(i,j,kmin))*(putin(i,j+1,kmin)-putin(i,j,kmin)) &
                  -(ekh(i,j,kmin)+ekh(i,j-1,kmin))*(putin(i,j,kmin)-putin(i,j-1,kmin)) )*dy2i &
                  + &
                ( rhobh(kmin+1)/rhobf(kmin) * (dzf(kmin+1)*ekh(i,j,kmin) + dzf(kmin)*ekh(i,j,kmin+1)) &
                  *  (putin(i,j,kmin+1)-putin(i,j,kmin)) / dzh(kmin+1)**2 &
                  + rhobh(kmin)/rhobf(kmin)*flux(i,j) *2.                        )/dzf(kmin) &
                          )

      end do
    end do

  end subroutine diffc



  subroutine diffe(putout)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dx2i,dzf,dy2i,dzh,imax,jmax
    use modfields, only : e120,rhobf,rhobh
    use modmpi,    only : myidx,myidy
    use modruraldata, only : bc_height
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    integer             :: i,j,k,jm,jp,km,kp,kmin

    do j=2,j1
      jp=j+1
      jm=j-1

      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
        do k=(kmin+1),kmax
          kp=k+1
          km=k-1

          putout(i,j,k) = putout(i,j,k) &
                  +  ( &
              ((ekm(i+1,j,k)+ekm(i,j,k))*(e120(i+1,j,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k)))*dx2i &
                  + &
              ((ekm(i,jp,k)+ekm(i,j,k)) *(e120(i,jp,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i,jm,k)) *(e120(i,j,k)-e120(i,jm,k)) )*dy2i &
                  + &
              (rhobh(kp)/rhobf(k) * (dzf(kp)*ekm(i,j,k) + dzf(k)*ekm(i,j,kp)) &
              *(e120(i,j,kp)-e120(i,j,k)) / dzh(kp)**2 &
              - rhobh(k)/rhobf(k) * (dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km)) &
              *(e120(i,j,k)-e120(i,j,km)) / dzh(k)**2        )/dzf(k) &
                            )

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=kmin
  !     --------------------------------------------

    do j=2,j1
      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1

        putout(i,j,kmin) = putout(i,j,kmin) + &
            ( (ekm(i+1,j,kmin)+ekm(i,j,kmin))*(e120(i+1,j,kmin)-e120(i,j,kmin)) &
              -(ekm(i,j,kmin)+ekm(i-1,j,kmin))*(e120(i,j,kmin)-e120(i-1,j,kmin)) )*dx2i &
            + &
            ( (ekm(i,j+1,kmin)+ekm(i,j,kmin))*(e120(i,j+1,kmin)-e120(i,j,kmin)) &
              -(ekm(i,j,kmin)+ekm(i,j-1,kmin))*(e120(i,j,kmin)-e120(i,j-1,kmin)) )*dy2i &
            + &
              ( rhobh(kmin+1)/rhobf(kmin) * (dzf(kmin+1)*ekm(i,j,kmin) + dzf(kmin)*ekm(i,j,kmin+1)) &
              *  (e120(i,j,kmin+1)-e120(i,j,kmin)) / dzh(kmin+1)**2              )/dzf(kmin)

      end do
    end do

  end subroutine diffe


  subroutine diffu (putout)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dx2i,dzf,dy,dyi,dzh, cu,cv,imax,jmax
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modsurfdata,only : ustar
    use modmpi,     only : myidx,myidy
    use modruraldata, only : bc_height
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emmo,emom,emop,empo
    real                :: fu
    real                :: ucu, upcu
    integer             :: i,j,k,jm,jp,km,kp,kmin

    do j=2,j1
      jp=j+1
      jm=j-1

      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
        do k=kmin+1,kmax
          kp=k+1
          km=k-1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )

          emop = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i-1,j,kp) ) ) / &
                    ( 4.   * dzh(kp) )

          empo = 0.25 * ( &
                  ekm(i,j,k)+ekm(i,jp,k)+ekm(i-1,jp,k)+ekm(i-1,j,k)  )

          emmo = 0.25 * ( &
                  ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )


          putout(i,j,k) = putout(i,j,k) &
                  + &
                  ( ekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k)) &
                    -ekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k)) ) * 2. * dx2i &
                  + &
                  ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxi) &
                    -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxi)   ) / dy &
                  + &
                  ( rhobh(kp)/rhobf(k) * emop * ( (u0(i,j,kp)-u0(i,j,k))   /dzh(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxi) &
                    - rhobh(k)/rhobf(k) * emom * ( (u0(i,j,k)-u0(i,j,km))   /dzh(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxi)   ) /dzf(k)

        end do
      end do
    end do

  !     ------------------------------------------------
  !     special treatment for lowest full level: k=kmin
  !     ------------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1

      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1

        empo = 0.25 * ( &
              ekm(i,j,kmin)+ekm(i,jp,kmin)+ekm(i-1,jp,kmin)+ekm(i-1,j,kmin)  )

        emmo = 0.25 * ( &
              ekm(i,j,kmin)+ekm(i,jm,kmin)+ekm(i-1,jm,kmin)+ekm(i-1,j,kmin)  )

        emop = ( dzf(kmin+1) * ( ekm(i,j,kmin) + ekm(i-1,j,kmin) )  + &
                    dzf(kmin) * ( ekm(i,j,kmin+1) + ekm(i-1,j,kmin+1) ) ) / &
                  ( 4.   * dzh(kmin+1) )


        ucu   = 0.5*(u0(i,j,kmin)+u0(i+1,j,kmin))+cu

        if(ucu >= 0.) then
          upcu  = max(ucu,1.e-10)
        else
          upcu  = min(ucu,-1.e-10)
        end if


        fu = ( 0.5*( ustar(i,j)+ustar(i-1,j) ) )**2  * &
                upcu/sqrt(upcu**2  + &
                ((v0(i,j,kmin)+v0(i-1,j,kmin)+v0(i,jp,kmin)+v0(i-1,jp,kmin))/4.+cv)**2)

        putout(i,j,kmin) = putout(i,j,kmin) &
                + &
              ( ekm(i,j,kmin)  * (u0(i+1,j,kmin)-u0(i,j,kmin)) &
              -ekm(i-1,j,kmin)* (u0(i,j,kmin)-u0(i-1,j,kmin)) ) * 2. * dx2i &
                + &
              ( empo * ( (u0(i,jp,kmin)-u0(i,j,kmin))   *dyi &
                        +(v0(i,jp,kmin)-v0(i-1,jp,kmin))*dxi) &
              -emmo * ( (u0(i,j,kmin)-u0(i,jm,kmin))   *dyi &
                        +(v0(i,j,kmin)-v0(i-1,j,kmin))  *dxi)   ) / dy &
               + &
              ( rhobh(kmin+1)/rhobf(kmin) * emop * ( (u0(i,j,kmin+1)-u0(i,j,kmin))    /dzh(kmin+1) &
                        +(w0(i,j,kmin+1)-w0(i-1,j,kmin+1))  *dxi) &
                -rhobh(kmin)/rhobf(kmin)*fu   ) / dzf(kmin)

      end do
    end do

  end subroutine diffu


  subroutine diffv (putout)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dx,dxi,dzf,dyi,dy2i,dzh, cu,cv,imax,jmax
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modsurfdata,only : ustar
    use modmpi,     only : myidx,myidy
    use modruraldata, only : bc_height

    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emmo, eomm,eomp,epmo
    real                :: fv, vcv,vpcv
    integer             :: i,j,k,jm,jp,km,kp,kmin

    do j=2,j1
      jp=j+1
      jm=j-1

      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
        do k=(kmin+1),kmax
          kp=k+1
          km=k-1

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
                    ( 4.   * dzh(k) )

          eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) / &
                    ( 4.   * dzh(kp) )

          emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )

          epmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,jm,k)+ekm(i+1,j,k)  )


        putout(i,j,k) = putout(i,j,k) &
                + &
              ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxi &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxi &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi)   ) / dx &
                + &
              (ekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
              -ekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &
                + &
              ( rhobh(kp)/rhobf(k) * eomp * ( (v0(i,j,kp)-v0(i,j,k))    /dzh(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                - rhobh(k)/rhobf(k) * eomm * ( (v0(i,j,k)-v0(i,j,km))    /dzh(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi)   ) / dzf(k)

        end do
      end do
    end do

  !     -----------------------------------------------
  !     special treatment for lowest full level: k=kmin
  !     -----------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1
      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1

        emmo = 0.25 * ( &
              ekm(i,j,kmin)+ekm(i,jm,kmin)+ekm(i-1,jm,kmin)+ekm(i-1,j,kmin)  )

        epmo = 0.25  * ( &
              ekm(i,j,kmin)+ekm(i,jm,kmin)+ekm(i+1,jm,kmin)+ekm(i+1,j,kmin)  )

        eomp = ( dzf(kmin+1) * ( ekm(i,j,kmin) + ekm(i,jm,kmin)  )  + &
                    dzf(kmin) * ( ekm(i,j,kmin+1) + ekm(i,jm,kmin+1) ) ) / &
                  ( 4.   * dzh(kmin+1) )

        vcv   = 0.5*(v0(i,j,kmin)+v0(i,j+1,kmin))+cv
        if(vcv >= 0.) then
          vpcv  = max(vcv,1.e-10)
        else
          vpcv  = min(vcv,-1.e-10)
        end if


        fv    = ( 0.5*( ustar(i,j)+ustar(i,j-1) ) )**2  * &
                    vpcv/sqrt(vpcv**2  + &
                ((u0(i,j,kmin)+u0(i+1,j,kmin)+u0(i,jm,kmin)+u0(i+1,jm,kmin))/4.+cu)**2)

        putout(i,j,kmin) = putout(i,j,kmin) &
                  + &
                  ( epmo * ( (v0(i+1,j,kmin)-v0(i,j,kmin))   *dxi &
                            +(u0(i+1,j,kmin)-u0(i+1,jm,kmin))*dyi) &
                    -emmo * ( (v0(i,j,kmin)-v0(i-1,j,kmin))   *dxi &
                            +(u0(i,j,kmin)-u0(i,jm,kmin))    *dyi)   ) / dx &
                  + &
                ( ekm(i,j,kmin) * (v0(i,jp,kmin)-v0(i,j,kmin)) &
                  -ekm(i,jm,kmin)* (v0(i,j,kmin)-v0(i,jm,kmin))  ) * 2. * dy2i &
                  + &
                ( rhobh(kmin+1)/rhobf(kmin) * eomp * ( (v0(i,j,kmin+1)-v0(i,j,kmin))     /dzh(kmin+1) &
                          +(w0(i,j,kmin+1)-w0(i,jm,kmin+1))    *dyi) &
                  -rhobh(kmin)/rhobf(kmin)*fv   ) / dzf(kmin)

      end do
    end do

  end subroutine diffv



  subroutine diffw(putout)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dx,dxi,dy,dyi,dzf,dzh,imax,jmax
    use modfields, only : u0,v0,w0,rhobh,rhobf
    use modmpi,     only : myidx,myidy
    use modruraldata, only : bc_height
    implicit none

  !*****************************************************************

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emom, eomm, eopm, epom
    integer             :: i,j,k,jm,jp,km,kp,kmin

    do j=2,j1
      jp=j+1
      jm=j-1
      do i=2,i1
        kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1

        do k=(kmin+1),kmax
          kp=k+1
          km=k-1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
                    ( 4.   * dzh(k) )

          eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) / &
                    ( 4.   * dzh(k) )

          epom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i+1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i+1,j,km) ) ) / &
                    ( 4.   * dzh(k) )


          putout(i,j,k) = putout(i,j,k) &
                + &
                  ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxi &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) /dzh(k) ) &
                    -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxi &
                            +(u0(i,j,k)-u0(i,j,km))     /dzh(k) ))/dx &
                + &
                  ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   /dzh(k) ) &
                    -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     /dzh(k) ))/dy &
                + (1./rhobh(k))*&
                  ( rhobf(k) * ekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) /dzf(k) &
                  - rhobf(km) * ekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) /dzf(km) ) * 2. &
                                                              / dzh(k)

        end do
      end do
    end do

  end subroutine diffw

end module
