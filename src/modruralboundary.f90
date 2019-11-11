!> \file modruralboundary.f90
!! By Michael Koene, TU Delft, section Atmospheric Physics, October 8, 2019
!! TODO: Write comments to output file
!!       clean up code (write statements)
!!       Test performance

module modruralboundary

  use modruraldata, only : lruralboundary, ldefrural, lnoslip, lwallfunc, damping, bc_height
  implicit none
  save
  public :: initruralboundary, exitruralboundary,&
            applyruralboundary

  ! Fields
  logical, allocatable :: limmersed_boundary(:,:,:)        !< Boolean array where .true. is the immersed boundary

  !< Shear layers in x,y,z-directions
  logical, allocatable :: lshear_x(:,:,:), lshear_y(:,:,:), lshear_z(:,:,:)
  !< Normal immersed boundary layers for incoming x,y,z-velocities
  logical, allocatable :: lnorm_x(:,:,:), lnorm_y(:,:,:), lnorm_z(:,:,:)

  !< Field for the minimal distance between a cellcenter and a wall
  real, allocatable    :: mindist(:,:,:)

  !< Field for the wall shear stress
  real, allocatable    :: shear(:,:,:,:)

  !< Fields for velocities at the ghost positions (k=0)
  real, allocatable    :: u0g(:,:), umg(:,:), v0g(:,:), vmg(:,:)

  !< Field for the ekm at the ghost positions (k=0)
  real, allocatable    :: ekmg(:,:)

contains
  subroutine initruralboundary
    use modglobal,  only : itot, jtot, ih, i1, jh, j1, imax, jmax, kmax, cexpnr, ifnamopt, ifinput, fname_options
    use modmpi,     only : myid, mpi_logical, comm3d, mpierr, MPI_INTEGER
    implicit none

    integer       :: i, j, k, ierr
    character(600) :: readstring

    namelist/NAMRURALBOUNDARY/ lruralboundary, ldefrural, &
                               lwallfunc, lnoslip


    if(myid==0) then    !first myid
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMRURALBOUNDARY,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMRURALBOUNDARY'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMRURALBOUNDARY'
      endif
      write(6 ,NAMRURALBOUNDARY)
      close(ifnamopt)
    endif

    !if(.not.(myid==0)) return

    call MPI_BCAST(lruralboundary   ,    1, mpi_logical , 0, comm3d, mpierr)

    if (.not. (lruralboundary)) return

    call MPI_BCAST(lwallfunc        ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(lnoslip          ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(ldefrural        ,    1, mpi_logical , 0, comm3d, mpierr)

    if (lnoslip .and. lwallfunc) then
      print *, 'Problem in namoptions NAMRURALBOUNDARY'
      print *, 'Cannot use both no slip conditions and wall functions for the shear'
      print *, 'Either set lnoslip to true or lwallfunc to true but not both.'
      stop 'ERROR: Problem in namoptions NAMRURALBOUNDARY'
    endif

    if (.not. (lnoslip .or. lwallfunc)) then
      print *, 'Problem in namoptions NAMRURALBOUNDARY'
      print *, 'Cannot go without no slip conditions or wall functions for the shear'
      print *, 'Either set lnoslip to true or lwallfunc to true by the use of lnoslip = .true. or lwallfunc = .true.'
      stop 'ERROR: Problem in namoptions NAMRURALBOUNDARY'
    endif

    write(6,*) 'allocating fields in ruralboundary'

    allocate(limmersed_boundary (itot+1,jtot+1,kmax))
    if (lnoslip) then
      allocate(lshear_x (2-ih:imax+ih,2-jh:jmax+jh,kmax))
      allocate(lshear_y (2-ih:imax+ih,2-jh:jmax+jh,kmax))
      allocate(lshear_z (2-ih:imax+ih,2-jh:jmax+jh,kmax))
    endif
    allocate(lnorm_x (2-ih:imax+ih,2-jh:jmax+jh,kmax))
    allocate(lnorm_y (2-ih:imax+ih,2-jh:jmax+jh,kmax))
    allocate(lnorm_z (2-ih:imax+ih,2-jh:jmax+jh,kmax))
    if (lwallfunc) then
      allocate(shear(2-ih:i1+ih,2-jh:j1+jh,kmax,12))
      allocate(damping(2-ih:i1+ih,2-jh:j1+jh,kmax))
      allocate(u0g(2-ih:i1+ih,2-jh:j1+jh))
      allocate(umg(2-ih:i1+ih,2-jh:j1+jh))
      allocate(v0g(2-ih:i1+ih,2-jh:j1+jh))
      allocate(vmg(2-ih:i1+ih,2-jh:j1+jh))
      allocate(ekmg(2-ih:i1+ih,2-jh:j1+jh))
    endif

    write(6,*) 'succesfully allocated fields in ruralboundary'

    limmersed_boundary(:,:,:) = .false.
    if (lwallfunc) then
      damping(:,:,:)=1.
      shear(:,:,:,:)=0.
      u0g(:,:)=0.
      umg(:,:)=0.
      v0g(:,:)=0.
      vmg(:,:)=0.
    endif


    ! Definition of ruralboundary
    if (myid==0) then
      if (ldefrural) then  !< Profile prescribed by use in the file rural_bc.inp.<expnr>
        write(6,*) 'Reading inputfile in ruralboundary'
        open (ifinput,file='rural_bc.inp.'//cexpnr)
          do i=1,itot
            read (ifinput,'(a200)') readstring

            do while (readstring(1:1)=='#')  ! Skip the lines that are commented (like headers)
              read (ifinput,'(a200)') readstring
            end do
            read(readstring,*) (bc_height(i,j),j=1,jtot)
          end do
        close(ifinput)

        do i=1,itot
          do j=1,jtot
            do k=1,kmax
              if (k.LE.bc_height(i,j)) then
                limmersed_boundary(i+1,j+1,k)=.true.
              endif
            end do
          end do
        end do
        write(6,*) 'Succesfully read inputfile in ruralboundary'
      else           !< Simple block in the middle of the domain
      write(6,*) 'Generating standard boundary in ruralboundary'
        bc_height(NINT(itot*0.5):(NINT(itot*0.5)+1),NINT(jtot*0.5):(NINT(jtot*0.5)+1))=NINT(kmax*0.5)
        do i=1,itot
          do j=1,jtot
            do k=1,kmax
              if (k.LE.bc_height(i,j)) then
                limmersed_boundary(i+1,j+1,k)=.true.
              endif
            end do
          end do
        end do
        write(6,*) 'Succesfully generated immersed boundary in ruralboundary'
      endif
    endif
    limmersed_boundary(1,:,:)=limmersed_boundary(itot+1,:,:)
    limmersed_boundary(:,1,:)=limmersed_boundary(:,jtot+1,:)
    bc_height(itot+1,:)=bc_height(1,:)
    bc_height(:,jtot+1)=bc_height(:,1)
    call MPI_BCAST(bc_height,(itot+1)*(jtot+1),MPI_INTEGER ,0,comm3d,mpierr)
    call MPI_BCAST(limmersed_boundary,itot*jtot*kmax,MPI_LOGICAL ,0,comm3d,mpierr)

    call constructboundarytypes

    !if (lwallfunc) call mindistance
    return
  end subroutine initruralboundary

  subroutine constructboundarytypes   !< Calculate the positions of the different boundary layers in multiple directions

    use modglobal,  only : imax, i1, ih, jmax, j1, jh, kmax, k1
    use modmpi,     only : myid, myidx, myidy, boolexcjs
    implicit none
    integer i,j,k,ipos,jpos


    write(6,*) 'Starting constructboundarytypes in ruralboundary, myid =',myid

    !< Fill the layer types with .false.
    if(lnoslip) then
      lshear_x(:,:,:)=.false.
      lshear_y(:,:,:)=.false.
      lshear_z(:,:,:)=.false.
    endif
    lnorm_x(:,:,:)=.false.
    lnorm_y(:,:,:)=.false.
    lnorm_z(:,:,:)=.false.

    !< Find Shear layers perpendicular to the x-direction and normal layers in x-direction
    if(lnoslip) then
      do k=1,kmax
        do j=2,jmax
          do i=2,imax-1
            ipos=i+myidx*imax
            jpos=j+myidy*jmax
            if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos+1,jpos,k))) then
              lshear_x(i,j,k)=limmersed_boundary(ipos+1,jpos,k)
              lshear_x(i+1,j,k)=limmersed_boundary(ipos,jpos,k)
              lnorm_x(i+1,j,k)=.true.
            endif
          end do
        end do
      end do
      write(6,*) 'Succesfully found shear and normal layers in x-direction'

      !< Find Shear layers perpendicular to the y-direction and normal layers in y-direction
      do k=1,kmax
        do i=2,imax
          do j=2,jmax-1
            ipos=i+myidx*imax
            jpos=j+myidy*jmax
            if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos+1,k))) then
              lshear_y(i,j,k)=limmersed_boundary(ipos,jpos+1,k)
              lshear_y(i,j+1,k)=limmersed_boundary(ipos,jpos,k)
              lnorm_y(i,j+1,k)=.true.
            endif
          end do
        end do
      end do

      !< Find Shear layers perpendicular to the z-direction and normal layers in z-direction
      do i=2,imax
        do j=2,jmax
          do k=1,kmax-1
            ipos=i+myidx*imax
            jpos=j+myidy*jmax
            if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos,k+1))) then
              !lshear_z(i,j,k)=limmersed_boundary(ipos,jpos,k+1)
              !lshear_z(i,j,k+1)=limmersed_boundary(ipos,jpos,k)
              !lnorm_z(i,j,k+1)=.true.
            endif
          end do
        end do
      end do
      call boolexcjs( lnorm_x  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lnorm_y  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lnorm_z  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lshear_x  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lshear_y  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lshear_z  , 2,imax,2,jmax,1,kmax,ih,jh)

    write(6,*) 'Succesfully found shear and normal layers in all directions'
    elseif(lwallfunc) then
      !< Find normal layers in x-direction
      do k=1,kmax
        do j=1,jmax
          do i=2,imax
            ipos=i+myidx*imax
            jpos=j+myidy*jmax
            if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos-1,jpos,k))) then
              lnorm_x(i,j,k)=.true.
            endif
          end do
        end do
      end do
      write(6,*) 'Succesfully found normal layers in x-direction'
      !write(6,*) 'jmax = ',jmax

      !< Find normal layers in y-direction
      do k=1,kmax
        do i=1,imax
          do j=2,jmax
            ipos=i+myidx*imax
            jpos=j+myidy*jmax
            if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos-1,k))) then
              lnorm_y(i,j,k)=.true.
            endif
          end do
        end do
      end do

      !< Find normal layers in z-direction
      do i=1,imax
        do j=1,jmax
          do k=2,kmax
            ipos=i+myidx*imax
            jpos=j+myidy*jmax
            if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos,k-1))) then
              !lnorm_z(i,j,k)=.true.
            endif
          end do
        end do
      end do

      write(6,*) 'before exjs: , myid=',myid
      write(6,*) 'lnorm_x(:,1,5)=',lnorm_x(:,1,5)
      write(6,*) 'lnorm_x(:,2,5)=',lnorm_x(:,2,5)
      write(6,*) 'lnorm_x(1,:,5)=',lnorm_x(1,:,5)

      call boolexcjs( lnorm_x  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lnorm_y  , 2,imax,2,jmax,1,kmax,ih,jh)
      call boolexcjs( lnorm_z  , 2,imax,2,jmax,1,kmax,ih,jh)

      !lnorm_x(itot+1,:,:)=lnorm_x(1,:,:)
      !lnorm_x(:,jtot+1,:)=lnorm_x(:,1,:)
      !lnorm_y(itot+1,:,:)=lnorm_y(1,:,:)
      !lnorm_y(:,jtot+1,:)=lnorm_y(:,1,:)
      !lnorm_z(itot+1,:,:)=lnorm_z(1,:,:)
      !lnorm_z(:,jtot+1,:)=lnorm_z(:,1,:)

      write(6,*) 'Succesfully found normal layers in all directions'
    endif
    write(6,*) 'finished constructboundarytypes'

    write(6,*) 'after exjs:, myid=',myid
    write(6,*) 'lnorm_x(:,1,5)=',lnorm_x(:,1,5)
    write(6,*) 'lnorm_x(:,2,5)=',lnorm_x(:,2,5)
    write(6,*) 'lnorm_x(1,:,5)=',lnorm_x(1,:,5)
    !write(6,*) 'lnorm_y(6,:,2)=',lnorm_y(6,:,2)
    !write(6,*) 'lnorm_x(:,5,2)=',lnorm_x(:,5,2)
    !write(6,*) 'lnorm_y(:,5,2)=',lnorm_y(:,5,2)
    !write(6,*) 'lnorm_z(9,14,:)=',lnorm_z(9,14,:)
  end subroutine constructboundarytypes

!  subroutine mindistance !< Determine the distance between the cellcenter and nearest wall for each cell (i,j,k)
!
!    use modglobal, only : itot, jtot, kmax, zh, zf, dx, dy
!    implicit none
!    integer :: i,j,k,iw,jw,kw
!    real    :: dist
!
!    mindist(:,:,:)=10*(itot+jtot+kmax)
!    write(6,*) 'Calculating minimal distances between each point in space and the walls'
!    do i=1,i1
!      do j=1,j1
!        do k=1,kmax
!          mindist(i,j,k)=zh(k)
!          do iw=1,itot+1
!            do jw=1,jtot+1
!              do kw=1,kmax
!                if(lnorm_x(iw,jw,kw)) then
!                  dist=sqrt(((iw-0.5-i)*dx)**2+((jw-j)*dy)**2+(zh(iw)-zh(i))**2)
!                  if (dist < mindist(i,j,k)) then
!                    mindist(i,j,k) = dist
!                  endif
!                elseif(lnorm_y(iw,jw,kw)) then
!                  dist=sqrt(((iw-i)*dx)**2+((jw-0.5-j)*dy)**2+(zh(iw)-zh(i))**2)
!                  if (dist < mindist(i,j,k)) then
!                    mindist(i,j,k) = dist
!                  endif
!                elseif(lnorm_z(iw,jw,kw))then
!                  dist=sqrt(((iw-i)*dx)**2+((jw-j)*dy)**2+(zf(iw)-zh(i))**2)
!                  if (dist < mindist(i,j,k)) then
!                    mindist(i,j,k) = dist
!                  endif
!                endif
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
!
!    return
!  end subroutine mindistance

  subroutine exitruralboundary
    use modmpi, only : myid
    implicit none

    deallocate(bc_height)
    !if (.not. (myid==0)) return
    if (.not. (lruralboundary)) return
    if(myid==0) write(6,*) 'Starting with exitruralboundary'
    if(lnoslip) then
      deallocate(lshear_x)
      deallocate(lshear_y)
      deallocate(lshear_z)
    elseif(lwallfunc) then
      !deallocate(mindist)
      !write(6,*) 'deallocating shear'
      deallocate(shear)
      !write(6,*) 'deallocating damping'
      deallocate(damping)
      !write(6,*) 'deallocating ghost points'
      deallocate(u0g)
      deallocate(umg)
      deallocate(v0g)
      deallocate(vmg)
      deallocate(ekmg)
    endif
    !write(6,*) 'deallocating lnorm_x'
    deallocate(lnorm_x)
    !write(6,*) 'deallocating lnorm_y'
    deallocate(lnorm_y)
    !write(6,*) 'deallocating lnorm_z'
    deallocate(lnorm_z)
    !write(6,*) 'deallocating bc_height'
    deallocate(limmersed_boundary)
    if(myid==0) write(6,*) 'Finished with exitruralboundary'
    return
  end subroutine exitruralboundary

  subroutine applyruralboundary   !< apply rural boundary conditions using the immersed boundary method
    use modfields,      only : um, vm, wm, u0, v0, w0, up, vp, wp
    use modglobal,      only : rk3step, itot, imax, jmax, jtot, kmax, i1, j1, k1, ih, jh, dt, rdt, timee, dx, dy, dzh, dzf, ekmin, nsv
    use modfields,      only : rhobf, rhobh, thl0, thlp, sv0, svp, e12p
    use modsubgriddata, only : ekm
    use modmicrodata,   only : nu_a
    use modmpi,         only : myid, excjs,myidx

    implicit none
    integer  :: i, j, k, m, nc
    real     :: rk3coef,rk3coefi
    real     :: emmo, emom, eomm, empo, epmo, epom, emop, eopm, eomp
    real     :: yplus

    integer  :: maxlocx(3)

    !if (.not. (myid==0)) return
    if (.not. (lruralboundary)) return

    if (myid==0) then
      !write(6,*) 'Starting with applyruralboundary, t = ',timee
      !write(6,*) 'The timestep dt = ',dt
    endif
    !write(6,*) 'um,vm,wm = ',um(2,2,2),',',vm(2,2,2),',',wm(2,2,2)

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    damping(:,:,:)=1.

    do i=1,i1  !1+myid*imax,myid*imax+1+imax
      do j=1,j1              !1+myid*jmax,myid*jmax+1+jmax
      !if(myid==1) write(6,*) 'myid=1: j=',j
        do k=2,kmax
          if(lnoslip) then
            write(6,*) 'lnoslip is used'
            !< Shear layer first
            if (lshear_x(i,j,k)) then
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
              wp(i,j,k)=-wm(i,j,k)*rk3coefi
            endif
            if (lshear_y(i,j,k)) then
              up(i,j,k)=-um(i,j,k)*rk3coefi
              wp(i,j,k)=-wm(i,j,k)*rk3coefi
            endif
            if (lshear_z(i,j,k)) then
              up(i,j,k)=-um(i,j,k)*rk3coefi
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
            endif
          elseif(lwallfunc) then         !< wallfunctions
            !write(6,*) 'start wallfunctions'
            !if(j==14) write(6,*) 'j==14 is reached'
            !if(j==14 .and. myid ==0) write(6,*) 'myid=0,j=14'
            !if(j==3 .and. myid ==1) write(6,*) 'myid=1,j=3'
            !if(j==3 .and. myid ==0) write(6,*) 'myid=0,j=3'
            !if(i==1 .and. myidx==1 .and. j==4 .and. k==5) write(6,*) 'i==1,myidx==1,lnorms are:',lnorm_x(i,j,k),lnorm_y(i,j,k),lnorm_z(i,j,k)
            !if(i==2 .and. myidx==1 .and. j==4 .and. k==5) write(6,*) 'i==2,myidx==1,lnorms are:',lnorm_x(i,j,k),lnorm_y(i,j,k),lnorm_z(i,j,k)
            !if(i==3 .and. myidx==1 .and. j==4 .and. k==5) write(6,*) 'i==3,myidx==1,lnorms are:',lnorm_x(i,j,k),lnorm_y(i,j,k),lnorm_z(i,j,k)
            !if(i==16 .and. myidx==0 .and. j==4 .and. k==5) write(6,*) 'i==16,myidx==0,lnorms are:',lnorm_x(i,j,k),lnorm_y(i,j,k),lnorm_z(i,j,k)
            !if(i==15 .and. myidx==0 .and. j==4 .and. k==5) write(6,*) 'i==15,myidx==0,lnorms are:',lnorm_x(i,j,k),lnorm_y(i,j,k),lnorm_z(i,j,k)
            !if(i==17 .and. myidx==0 .and. j==4 .and. k==5) write(6,*) 'i==17,myidx==0,lnorms are:',lnorm_x(i,j,k),lnorm_y(i,j,k),lnorm_z(i,j,k)
            if (lnorm_x(i,j,k)) then     !< Wall in x-direction
              !if(myidx==0 .and. j==4 .and. k==5) write(6,*) 'myidx=0 lnorm_x=T,i=',i
              !if(myidx==1 .and. j==4 .and. k==5) write(6,*) 'myidx=1 lnorm_x=T,i=',i
              emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,j-1,k)+ekm(i-1,j-1,k)+ekm(i-1,j,k)  )

              !write(6,*) 'norm_x at (i,j,k) = (',i,j,k,') and myid =',myid

            !if(mindist(i,j,k)-0.5*dx<1e-6) write(6,*) 'mindist equal to dx/2'
              call wallaw(v0(i-1,j,k),0.5*dx,nu_a,shear(i-1,j,k,1))
              call wallaw(v0(i,j,k),  0.5*dx,nu_a,shear(i,j,k,2))

              vp(i-1,j,k) = vp(i-1,j,k) - 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * shear(i-1,j,k,1)/dx
              vp(i,j,k)   = vp(i,j,k)   + 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * shear(i,j,k,2)  /dx

              empo = 0.25  * ( &
                ekm(i,j+1,k)+ekm(i,j,k)+ekm(i-1,j,k)+ekm(i-1,j+1,k)  )

              call wallaw(v0(i-1,j+1,k),0.5*dx,nu_a,shear(i-1,j+1,k,1))
              call wallaw(v0(i,j+1,k),  0.5*dx,nu_a,shear(i,j+1,k,2))

              vp(i-1,j+1,k) = vp(i-1,j+1,k) - 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * shear(i-1,j+1,k,1)/dx
              vp(i,j+1,k)   = vp(i,j+1,k)   + 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * shear(i,j+1,k,2)  /dx

              emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                    ( 4.   * dzh(k) )

              call wallaw(w0(i-1,j,k),0.5*dx,nu_a,shear(i-1,j,k,3))
              call wallaw(w0(i,j,k),  0.5*dx,nu_a,shear(i,j,k,4))

              wp(i-1,j,k) = wp(i-1,j,k) - 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * shear(i-1,j,k,3)/dx/dx
              wp(i,j,k)   = wp(i,j,k)   + 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * shear(i,j,k,4)  /dx/dx

              emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i-1,j,k+1)  )  + &
                      dzf(k+1)  * ( ekm(i,j,k) + ekm(i-1,j,k) ) ) / &
                    ( 4.   * dzh(k+1) )

              call wallaw(w0(i-1,j,k+1),0.5*dx,nu_a,shear(i-1,j,k+1,3))
              call wallaw(w0(i,j,k+1),  0.5*dx,nu_a,shear(i,j,k+1,4))

              wp(i-1,j,k+1) = wp(i-1,j,k+1) - 0.5 * emop * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * shear(i-1,j,k+1,3)/dx/dx
              wp(i,j,k+1)   = wp(i,j,k+1)   + 0.5 * emop * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * shear(i,j,k+1,4)  /dx/dx

              call xwallscalar(i,j,k,thl0,thlp)
              do nc=1,nsv
                call xwallscalar(i,j,k,sv0(:,:,:,nc),svp(:,:,:,nc))
              end do
              call xwalle12(i,j,k)

            endif
            if (lnorm_y(i,j,k)) then     !< Wall in y-direction
              emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,j-1,k)+ekm(i-1,j-1,k)+ekm(i-1,j,k)  )

              call wallaw(u0(i,j-1,k),0.5*dy,nu_a,shear(i,j-1,k,5))
              call wallaw(u0(i,j,k)  ,0.5*dy,nu_a,shear(i,j,k,6))

              up(i,j-1,k) = up(i,j-1,k) - 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j-1,k,5)/dy
              up(i,j,k)   = up(i,j,k)   + 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j,k,6)  /dy

              epmo = 0.25  * ( &
                ekm(i+1,j,k)+ekm(i+1,j-1,k)+ekm(i,j-1,k)+ekm(i,j,k)  )

              call wallaw(u0(i+1,j-1,k),0.5*dy,nu_a,shear(i+1,j-1,k,5))
              call wallaw(u0(i+1,j,k)  ,0.5*dy,nu_a,shear(i+1,j,k,6))

              up(i+1,j-1,k) = up(i+1,j-1,k) - 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * shear(i+1,j-1,k,5)/dy
              up(i+1,j,k)   = up(i+1,j,k)   + 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * shear(i+1,j,k,6)  /dy

              eomm = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / ( 4.  * dzh(k) )

              call wallaw(w0(i,j-1,k),0.5*dy,nu_a,shear(i,j-1,k,7))
              call wallaw(w0(i,j,k)  ,0.5*dy,nu_a,shear(i,j,k,8))

              wp(i,j-1,k) = wp(i,j-1,k) - 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j-1,k,7)/dy
              wp(i,j,k)   = wp(i,j,k)   + 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j,k,8)  /dy

              eomp = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i,j-1,k+1)  )  + &
                dzf(k+1) * ( ekm(i,j,k) + ekm(i,j-1,k) ) ) / ( 4.  * dzh(k+1) )

              call wallaw(w0(i,j-1,k+1),0.5*dy,nu_a,shear(i,j-1,k+1,7))
              call wallaw(w0(i,j,k+1)  ,0.5*dy,nu_a,shear(i,j,k+1,8))

              wp(i,j-1,k+1) = wp(i,j-1,k+1) - 0.5 * eomp * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * shear(i,j-1,k+1,7)/dy
              wp(i,j,k+1)   = wp(i,j,k+1)   + 0.5 * eomp * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * shear(i,j,k+1,8)  /dy

              call ywallscalar(i,j,k,thl0,thlp)
              do nc=1,nsv
                call ywallscalar(i,j,k,sv0(:,:,:,nc),svp(:,:,:,nc))
              end do
              call ywalle12(i,j,k)

            endif
            if (lnorm_z(i,j,k)) then     !< Wall in z-direction
              emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                    ( 4.   * dzh(k) )

              call wallaw(u0(i,j,k-1),0.5*dzf(k-1),nu_a,shear(i,j,k-1,9))
              call wallaw(u0(i,j,k)  ,0.5*dzf(k)  ,nu_a,shear(i,j,k,10))

              up(i,j,k-1) = up(i,j,k-1) - 0.5 * emom * rhobh(k)/rhobf(k-1) *((u0(i,j,k)-u0(i,j,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i,j,k-1,9)/dzf(k-1)
              up(i,j,k)   = up(i,j,k)   + 0.5 * emom * rhobh(k)/rhobf(k) *((u0(i,j,k)-u0(i,j,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i,j,k,10) /dzf(k)

              epom = ( dzf(k-1) * ( ekm(i+1,j,k)  + ekm(i,j,k)  )  + &
                      dzf(k)  * ( ekm(i+1,j,k-1) + ekm(i,j,k-1) ) ) / &
                    ( 4.   * dzh(k) )

              call wallaw(u0(i+1,j,k-1),0.5*dzf(k-1),nu_a,shear(i+1,j,k-1,9))
              call wallaw(u0(i+1,j,k)  ,0.5*dzf(k)  ,nu_a,shear(i+1,j,k,10))

              up(i+1,j,k-1) = up(i+1,j,k-1) - 0.5 * epom * rhobh(k)/rhobf(k-1) *((u0(i+1,j,k)-u0(i+1,j,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i+1,j,k-1,9)/dzf(k-1)
              up(i+1,j,k)   = up(i+1,j,k)   + 0.5 * epom * rhobh(k)/rhobf(k) *((u0(i+1,j,k)-u0(i+1,j,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i+1,j,k,10) /dzf(k)

              eomm = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / ( 4.  * dzh(k) )

              call wallaw(v0(i,j,k-1),0.5*dzf(k-1),nu_a,shear(i,j,k-1,11))
              call wallaw(v0(i,j,k)  ,0.5*dzf(k)  ,nu_a,shear(i,j,k,12))

              vp(i,j,k-1) = vp(i,j,k-1) - 0.5 * eomm * rhobh(k)/rhobf(k-1) *((v0(i,j,k)-v0(i,j,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i,j,k-1,11)/dzf(k-1)
              vp(i,j,k)   = vp(i,j,k)   + 0.5 * eomm * rhobh(k)/rhobf(k) *((v0(i,j,k)-v0(i,j,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i,j,k,12)  /dzf(k)

              eopm = ( dzf(k-1) * ( ekm(i,j+1,k)  + ekm(i,j,k)  )  + &
                dzf(k) * ( ekm(i,j+1,k-1) + ekm(i,j,k-1) ) ) / ( 4.  * dzh(k) )

              call wallaw(v0(i,j+1,k-1),0.5*dzf(k-1),nu_a,shear(i,j+1,k-1,11))
              call wallaw(v0(i,j+1,k)  ,0.5*dzf(k)  ,nu_a,shear(i,j+1,k,12))

              vp(i,j+1,k-1) = vp(i,j+1,k-1) - 0.5 * eopm * rhobh(k)/rhobf(k-1) *((v0(i,j+1,k)-v0(i,j+1,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i,j+1,k-1,11)/dzf(k-1)
              vp(i,j+1,k)   = vp(i,j+1,k)   + 0.5 * eopm * rhobh(k)/rhobf(k) *((v0(i,j+1,k)-v0(i,j+1,k-1))/dzh(k))/dzf(k-1) - 0.5 * shear(i,j+1,k,12)  /dzf(k)


              call zwallscalar(i,j,k,thl0,thlp)
              do nc=1,nsv
                call zwallscalar(i,j,k,sv0(:,:,:,nc),svp(:,:,:,nc))
              end do
              call zwalle12(i,j,k)

            endif
            !TODO: Maybe remove mindist and put in 0.5*dx,0.5*dy,0.5*dz depending on direction
            if (lnorm_x(i,j,k) .or. lnorm_y(i,j,k) .or. lnorm_z(i,j,k)) then !< Calculate and apply Piomelli wall damping
              !yplus= mindist(i,j,k)*sqrt(sum(abs(shear(i,j,k,:))))/nu_a
              if(lnorm_x(i,j,k)) then
                yplus = 0.5 * dx * sqrt(sum(abs(shear(i,j,k,1:4))))/nu_a
                damping(i,j,k)   = min(damping(i,j,k),  1.-exp(-(yplus*0.04)**3.))
                yplus = 0.5 * dx * sqrt(sum(abs(shear(i-1,j,k,1:4))))/nu_a
                damping(i-1,j,k) = min(damping(i-1,j,k),1.-exp(-(yplus*0.04)**3.))
              endif
              if(lnorm_y(i,j,k)) then
                yplus = 0.5 * dy * sqrt(sum(abs(shear(i,j,k,5:8))))/nu_a
                damping(i,j,k)   = min(damping(i,j,k),  1.-exp(-(yplus*0.04)**3.))
                yplus = 0.5 * dy * sqrt(sum(abs(shear(i,j-1,k,5:8))))/nu_a
                damping(i,j-1,k) = min(damping(i,j-1,k),1.-exp(-(yplus*0.04)**3.))
              endif
              if(lnorm_z(i,j,k)) then
                yplus = 0.5*dzf(k)*sqrt(sum(abs(shear(i,j,k,9:12))))/nu_a
                damping(i,j,k)   = min(damping(i,j,k),  1.-exp(-(yplus*0.04)**3.))
                yplus = 0.5*dzf(k)*sqrt(sum(abs(shear(i,j,k-1,9:12))))/nu_a
                damping(i,j,k-1) = min(damping(i,j,k-1),1.-exp(-(yplus*0.04)**3.))
              endif
              !write(6,*) 'yplus = ',yplus,', dx = ',dx, ', sum(abs(shear))=',sum(abs(shear(i,j,k,1:12)))
              !write(6,*) 'damping =',damping(i,j,k)
            endif
          endif
        end do
      end do
    end do
    if(lwallfunc) then !< special treatment for lowest full level: k=1
      do i=1,i1
        do j=1,j1

          call wallaw(u0(i,j,1),0.5*dzh(1),nu_a,shear(i,j,1,10))
          u0g(i,j)  = u0(i,j,1)  - shear(i,j,1,10)*dzf(1)/nu_a
          umg(i,j)  = um(i,j,1)  - shear(i,j,1,10)*dzf(1)/nu_a

          call wallaw(v0(i,j,1),0.5*dzh(1),nu_a,shear(i,j,1,12))
          v0g(i,j)  = v0(i,j,1)  - shear(i,j,1,12)*dzf(1)/nu_a
          vmg(i,j)  = -vm(i,j,1) - shear(i,j,1,12)*dzf(1)/nu_a

          ekmg(i,j) = -ekm(i,j,1) + 2.*nu_a

        end do
      end do
      do i=2,i1
        do j=2,j1

          if (lnorm_x(i,j,1)) then     !< Wall in x-direction
            emmo = 0.25  * ( &
                ekm(i,j,1)+ekm(i,j-1,1)+ekm(i-1,j-1,1)+ekm(i-1,j,1)  )

            !if(mindist(i,j,1)-0.5*dx<1e-6) write(6,*) 'mindist equal to dx/2'
            call wallaw(v0(i-1,j,1),0.5*dx,nu_a,shear(i-1,j,1,1))
            call wallaw(v0(i,j,1),  0.5*dx,nu_a,shear(i,j,1,2))

            vp(i-1,j,1) = vp(i-1,j,1) - 0.5 * emmo*((v0(i,j,1)-v0(i-1,j,1))/dx) / dx - 0.5 * shear(i-1,j,1,1)/dx
            vp(i,j,1)   = vp(i,j,1)   + 0.5 * emmo*((v0(i,j,1)-v0(i-1,j,1))/dx) / dx - 0.5 * shear(i,j,1,2)  /dx

            empo = 0.25  * ( &
                ekm(i,j+1,1)+ekm(i,j,1)+ekm(i-1,j,1)+ekm(i-1,j+1,1)  )

            !if(mindist(i,j,1)-0.5*dx<1e-6) write(6,*) 'mindist equal to dx/2'
            call wallaw(v0(i-1,j+1,1),0.5*dx,nu_a,shear(i-1,j+1,1,1))
            call wallaw(v0(i,j+1,1),  0.5*dx,nu_a,shear(i,j+1,1,2))

            vp(i-1,j+1,1) = vp(i-1,j+1,1) - 0.5 * emmo*((v0(i,j+1,1)-v0(i-1,j+1,1))/dx) / dx - 0.5 * shear(i-1,j+1,1,1)/dx
            vp(i,j+1,1)   = vp(i,j+1,1)   + 0.5 * emmo*((v0(i,j+1,1)-v0(i-1,j+1,1))/dx) / dx - 0.5 * shear(i,j+1,1,2)  /dx


            call xwallscalar(i,j,1,thl0,thlp)
            do nc=1,nsv
              call xwallscalar(i,j,1,sv0(:,:,:,nc),svp(:,:,:,nc))
            end do
            call xwalle12(i,j,1)
          endif

          if (lnorm_y(i,j,1)) then     !< Wall in y-direction
            emmo = 0.25  * ( &
              ekm(i,j,1)+ekm(i,j-1,1)+ekm(i-1,j-1,1)+ekm(i-1,j,1)  )

            call wallaw(u0(i,j-1,1),0.5*dy,nu_a,shear(i,j-1,1,5))
            call wallaw(u0(i,j,1)  ,0.5*dy,nu_a,shear(i,j,1,6))

            up(i,j-1,1) = up(i,j-1,1) - 0.5 * emmo * ((u0(i,j,1)-u0(i,j-1,1))/dy)/dy - 0.5 * shear(i,j-1,1,5)/dy
            up(i,j,1)   = up(i,j,1)   + 0.5 * emmo * ((u0(i,j,1)-u0(i,j-1,1))/dy)/dy - 0.5 * shear(i,j,1,6)  /dy

            epmo = 0.25  * ( &
              ekm(i+1,j,1)+ekm(i+1,j-1,1)+ekm(i,j-1,1)+ekm(i,j,1)  )

            call wallaw(u0(i+1,j-1,1),0.5*dy,nu_a,shear(i+1,j-1,1,5))
            call wallaw(u0(i+1,j,1)  ,0.5*dy,nu_a,shear(i+1,j,1,6))

            up(i+1,j-1,1) = up(i+1,j-1,1) - 0.5 * emmo * ((u0(i+1,j,1)-u0(i+1,j-1,1))/dy)/dy - 0.5 * shear(i+1,j-1,1,5)/dy
            up(i+1,j,1)   = up(i+1,j,1)   + 0.5 * emmo * ((u0(i+1,j,1)-u0(i+1,j-1,1))/dy)/dy - 0.5 * shear(i+1,j,1,6)  /dy


            call ywallscalar(i,j,1,thl0,thlp)
            do nc=1,nsv
              call ywallscalar(i,j,1,sv0(:,:,:,nc),svp(:,:,:,nc))
            end do
            call ywalle12(i,j,1)

          endif

          ! Wall shear with the ground
          ! write(6,*) 'i=',i
          !emom = ( dzh(1) * ( ekm(i,j,1)  + ekm(i-1,j,1)  )  + &
          !         dzf(1) * ( ekmg(i,j) + ekmg(i-1,j)     ) ) / &
          !       ( 4.   * dzh(1) )

          !up(i,j,1)   = up(i,j,1)  - shear(i,j,1,10) /dzf(1) + emom * rhobh(1)/rhobf(1) *((u0(i,j,1)-u0g(i,j))/dzh(1))/dzh(1)

          !eomm = ( dzh(1) * ( ekm(i,j,1)  + ekm(i,j-1,1)  )  + &
          !         dzf(1) * ( ekmg(i,j)   + ekmg(i,j-1)   ) ) / ( 4.  * dzh(1) )

          !vp(i,j,1)   = vp(i,j,1) - shear(i,j,1,12)  /dzf(1)   + eomm * rhobh(1)/rhobf(1) *((v0(i,j,1)-v0g(i,j))/dzh(1))/dzh(1)

          !yplus = 0.5*dzh(1)*sqrt(sum(abs(shear(i,j,1,9:12))))/nu_a
          !damping(i,j,1)   = min(damping(i,j,1),  1.-exp(-(yplus*0.04)**3.))
        end do
      end do
      do i=1,i1
        do j=1,j1
          do k=1,kmax
            !< Now the normal velocities
            if (lnorm_x(i,j,k)) then
              !write(6,*) 'lnorm_x found, um = ',um(i,j,k)
              !if((i==2 .and. j==4) .and. k==1) write(6,*) 'Before: u0 = ',u0(i,j,k),' um = ',um(i,j,k),' up = ',up(i,j,k)
              up(i,j,k)=-um(i,j,k)*rk3coefi
              !if((i==2 .and. j==4) .and. k==1) write(6,*) 'After: u0 = ',u0(i,j,k),' um = ',um(i,j,k),' up = ',up(i,j,k),' rk3coef = ',rk3coef
            end if
            if (lnorm_y(i,j,k)) then
              !write(6,*) 'lnorm_y found, vm = ',vm(i,j,k)
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
            end if
            if (lnorm_z(i,j,k)) then
              !write(6,*) 'lnorm_z found, wm = ',wm(i,j,k)
              wp(i,j,k)=-wm(i,j,k)*rk3coefi
            end if
          end do
        end do
      end do
    endif
    !write(6,*) 'Finishing applyruralboundary, max up = ',maxval(up),', vp = ',maxval(vp),', wp = ',maxval(wp)
    !write(6,*) 'Reached the excjs part'
    call excjs( up  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( vp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( wp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( damping  , 2,i1,2,j1,1,kmax,ih,jh)
    call excjs( e12p  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( thlp  , 2,i1,2,j1,1,k1,ih,jh)
    do nc=1,nsv
      call excjs( svp(:,:,:,nc)  , 2,i1,2,j1,1,k1,ih,jh)
    enddo

    !write(6,*) 'Values at 2,2,10: u0=',up(2,2,10),' vp=',vp(2,2,10),'wp=',wp(2,2,10),'damping=',damping(2,2,10),'thlp=',thlp(2,2,10),'rhobh=',rhobh(10),'thl0=',thl0(2,2,10)

    if(maxval(sqrt(u0**2.+v0**2.+w0**2.))>16) then
      maxlocx=maxloc(u0**2.+v0**2.+w0**2.)
      write(6,*) 'maxlocx = ',maxlocx(1),maxlocx(2),maxlocx(3)
      !maxlocy=maxloc(u0**2.+v0**2.+w0**2.,dim=2)
      !maxlocz=maxloc(u0**2.+v0**2.+w0**2.,dim=3)
      write(6,*) 'ERROR: vel>16, maxloc = ',maxloc(sqrt(u0**2.+v0**2.+w0**2.)), 'u0 and um and up are here:',u0(maxlocx(1),maxlocx(2),maxlocx(3)),um(maxlocx(1),maxlocx(2),maxlocx(3)),up(maxlocx(1),maxlocx(2),maxlocx(3))
      write(6,*) 'v0 and vm and vp are here:',v0(maxlocx(1),maxlocx(2),maxlocx(3)),vm(maxlocx(1),maxlocx(2),maxlocx(3)),vp(maxlocx(1),maxlocx(2),maxlocx(3))
      write(6,*) 'w0 and wm and wp are here:',w0(maxlocx(1),maxlocx(2),maxlocx(3)),wm(maxlocx(1),maxlocx(2),maxlocx(3)),wp(maxlocx(1),maxlocx(2),maxlocx(3))
      write(6,*) 'shear(maxloc) = ',shear(maxlocx(1),maxlocx(2),maxlocx(3),:)
    endif
    !write(6,*) 'Max shear for x= ',maxval(shear(:,:,:,1:4)),' y= ',maxval(shear(:,:,:,5:8)),' z= ',maxval(shear(:,:,:,9:12))
    !write(6,*) 'Min damping = ',minval(damping)
    !if(myid==0) write(6,*) 'Finished apply ruralboundary'
    return
  end subroutine applyruralboundary

  subroutine xwallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dx2i
    use modsubgriddata, only : ekh

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * (ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k))*dx2i
    putout(i-1,j,k) = putout(i-1,j,k) &
                      - 0.5 * (ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k))*dx2i

    return
  end subroutine xwallscalar

  subroutine ywallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dy2i
    use modsubgriddata, only : ekh

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * (ekh(i,j,k)+ekh(i,j-1,k)) *(putin(i,j,k)-putin(i,j-1,k))*dy2i
    putout(i,j-1,k) = putout(i,j-1,k) &
                      - 0.5 * (ekh(i,j,k)+ekh(i,j-1,k)) *(putin(i,j,k)-putin(i,j-1,k))*dy2i

    return
  end subroutine ywallscalar

  subroutine zwallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dzf, dzh
    use modsubgriddata, only : ekh
    use modfields,      only : rhobh, rhobf

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * rhobh(k)/rhobf(k-1) * (dzf(k-1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k-1)) &
                      *  (putin(i,j,k)-putin(i,j,k-1)) / dzh(k)**2 /dzf(k)
    putout(i,j,k-1) = putout(i,j,k-1) &
                      - 0.5 * rhobh(k)/rhobf(k) * (dzf(k-1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k-1)) &
                      *  (putin(i,j,k)-putin(i,j,k-1)) / dzh(k)**2 /dzf(k)

    return
  end subroutine zwallscalar

  subroutine xwalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dx2i, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not.(k==1)) then
      e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
                                     -((w0(i,j,k+1)-w0(i-1,j,k+1))  / dx             + &
                                       (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &
                                     +(2.*(w0(i,j,k+1))             / dx             + &
                                       (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &

                                     -((w0(i,j,k)-w0(i-1,j,k))      / dx             + &
                                       (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &
                                     +(2.*(w0(i,j,k))               / dx             + &
                                       (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &

                                     -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
                                     +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (2.*v0(i,j+1,k))             / dx       )**2  + &

                                     -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
                                     +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (2.*v0(i,j,k))               / dx       )**2    &
                                    )

      e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms
                                       -((w0(i,j,k)-w0(i-1,j,k))    / dx             + &
                                         (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &
                                       +((-2.*w0(i-1,j,k))          / dx             + &
                                         (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &

                                       -((w0(i,j,k+1)-w0(i-1,j,k+1))/ dx             + &
                                         (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &
                                       +((-2.*w0(i-1,j,k+1))        / dx             + &
                                         (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &

                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (-2.*v0(i-1,j,k))          / dx       )**2  + &

                                       -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
                                       +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (-2.*v0(i-1,j+1,k))        / dx       )**2    &
                                    )
    else !Special treatment for the lowest full level: k=1
      e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
                                     -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
                                     +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (2.*v0(i,j+1,k))             / dx       )**2  + &

                                     -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
                                     +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (2.*v0(i,j,k))               / dx       )**2    &
                                    )

      e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms

                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (-2.*v0(i-1,j,k))          / dx       )**2  + &

                                       -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
                                       +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (-2.*v0(i-1,j+1,k))        / dx       )**2    &
                                    )
    endif
  end subroutine xwalle12


  subroutine ywalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dy2i, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not.(k==1)) then
      e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((2.*u0(i,j,k))             / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((2.*u0(i+1,j,k))           / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

                                       -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
                                       +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (2.*w0(i,j,k+1))           / dy       )**2  + &

                                       -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
                                       +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (2.*w0(i,j,k))             / dy       )**2    &
                                    )

      e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i,j-1,k))          / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i+1,j-1,k))        / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

                                       -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
                                       +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (-2.*w0(i,j-1,k))          / dy       )**2  + &

                                       -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
                                       +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (-2.*w0(i,j-1,k+1))        / dy       )**2    &
                                    )
    else !Special treatment for the lowest full level: k=1
      e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((2.*u0(i,j,k))             / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((2.*u0(i+1,j,k))           / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
                                    )

      e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i,j-1,k))          / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i+1,j-1,k))        / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
                                    )
    endif
  end subroutine ywalle12

   subroutine zwalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dzf, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0, rhobh, rhobf

    implicit none

    integer, intent(in)    :: i,j,k

    e12p(i,j,k)   = e12p(i,j,k) - rhobh(k)/rhobf(k-1) * (dzf(k)*ekm(i,j,k-1) + dzf(k-1)*ekm(i,j,k)) &
                                  *(e120(i,j,k)-e120(i,j,k-1)) / dzh(k)**2 /dzf(k)     &
                                + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                   -((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (u0(i,j,k)-u0(i,j,k-1))        / dzh(k)   )**2  + &
                                   +((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (2.*u0(i,j,k))                 / dzh(k)   )**2  + &

                                   -((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (u0(i+1,j,k)-u0(i+1,j,k-1))    / dzh(k)   )**2  + &
                                   +((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (2.*u0(i+1,j,k))               / dzh(k)   )**2  + &

                                   -((v0(i,j,k)-v0(i,j,k-1))        / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &
                                   +((2.*v0(i,j,k))                 / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &

                                   -((v0(i,j+1,k)-v0(i,j+1,k-1))    / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2  + &
                                   +(2.*v0(i,j+1,k)                 / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2    &
                                )

    e12p(i,j,k-1) = e12p(i,j,k-1) + rhobh(k)/rhobf(k) * (dzf(k-1)*ekm(i,j,k) + dzf(k)*ekm(i,j,k-1)) &
                                    *(e120(i,j,k)-e120(i,j,k-1)) / dzh(k)**2  /dzf(k)  &
                                  + ekm(i,j,k-1)/(2.*e120(i,j,k-1))* (&  !source terms
                                   -((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (u0(i,j,k)-u0(i,j,k-1))        / dzh(k)   )**2  + &
                                   +((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (-2.*u0(i,j,k-1))              / dzh(k)   )**2  + &

                                   -((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (u0(i+1,j,k)-u0(i+1,j,k-1))    / dzh(k)   )**2  + &
                                   +((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (-2.*u0(i+1,j,k-1))            / dzh(k)   )**2  + &
                                   -((v0(i,j,k)-v0(i,j,k-1))        / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &
                                   +((-2.*v0(i,j,k-1))              / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &

                                   -((v0(i,j+1,k)-v0(i,j+1,k-1))    / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2  + &
                                   +((-2.*v0(i,j+1,k-1))            / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2    &
                                )
  end subroutine zwalle12

  subroutine wallaw(utan,dx,visc,tau)  !< Copied from work by Jasper Tomas
    !use modglobal,         only :

    implicit none

    real, intent(in)  :: utan  ! tangential velocity component
    real, intent(in)  :: dx    ! distance to the wall
    real, intent(in)  :: visc  ! viscosity near the wall
    real, intent(out) :: tau   !

    real dxi,dx5
    real utanabs, utankr, dutan, sub
    real tausub, taupow

    real const1, const2, const3, const4
    real aaa, bbb

    parameter(aaa = 8.3)
    parameter(bbb = 0.1428571429)

    !write(6,*) 'Starting with wallaw'

    dxi = 1./dx
    dx5 = 0.5*dx

    const1 = 0.5 * (1. - bbb) * aaa ** ((1. + bbb) / (1. - bbb))
    const2 = (1. + bbb) / aaa
    const3 = aaa ** (2. / (1. - bbb))
    const4 = 2. / (1. + bbb)

    utanabs = abs(utan)
    utankr  = 0.5 * visc * dxi * const3
    dutan   = utankr - utanabs
    sub     = max (sign(1.,dutan),0.) ! sub = 1 for viscous sublayer and 0 for outside.

    tausub  = 2. * visc * utanabs * dxi
    taupow  = (const1 * (visc * dxi) ** (1. + bbb) + const2 * ((visc * dxi) ** bbb) * utanabs) ** const4

    tau = sub * tausub + (1. - sub) * taupow ! use linear profile inside viscous sublayer and use power law outside
    tau = sign(tau,utan)   ! give tau the same sign as utan
    return
  end subroutine wallaw

end module modruralboundary
