
subroutine qinit(meqn,mbc,mx,my,xlow_patch,ylow_patch,dx,dy,q,maux,aux)

    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level, coordinate_system
    use topo_module
    use vp2_module

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow_patch,ylow_patch,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i,j,m,mx2,my2,mcapa
    real(kind=8) :: xlow2,ylow2,xhi2,yhi2,dx2,dy2
    real(kind=8) :: xim,xip,yjm,yjp,xipc,ximc,yjmc,yjpc,dq,x,y,topointegral
    real(kind=8), allocatable :: init_remold(:)

    ! Set flat state based on sea_level
    q = 0.d0

    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlow_patch,ylow_patch,dx,dy,q,maux,aux)
    endif

    if (fname_remold .ne. " ") then
        call read_topo_header(fname_remold, 3 ,mx2,my2, &
                 xlow2,ylow2,xhi2,yhi2,dx2,dy2)

        if (.not.allocated(init_remold)) allocate(init_remold(mx2*my2))

        call read_topo_file(mx2,my2, 3 ,fname_remold,xlow2,ylow2,init_remold)
        
        do i=1,mx
            x = xlow_patch + (i-0.5d0)*dx
            xim = x - 0.5d0*dx
            xip = x + 0.5d0*dx
            do j=1,my
                y = ylow_patch + (j-0.5d0)*dy
                yjm = y - 0.5d0*dy
                yjp = y + 0.5d0*dy

                ! Check to see if we are in the qinit region at this grid point
                if ((xip > xlow2).and.(xim < xhi2).and.  &
                    (yjp > ylow2).and.(yjm < yhi2)) then

                    xipc=min(xip,xhi2)
                    ximc=max(xim,xlow2)

                    yjpc=min(yjp,yhi2)
                    yjmc=max(yjm,ylow2)

                    dq = topointegral(ximc,xipc,yjmc,yjpc,xlow2, &
                                      ylow2,dx2,dy2,mx2, &
                                      my2,init_remold,1)
                    if (coordinate_system == 2) then
                        dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(mcapa,i,j))
                    else
                        dq = dq / ((xipc-ximc)*(yjpc-yjmc))
                    endif

                    q(6,i,j) = q(6,i,j) + dq

                endif
            enddo
        enddo
    endif

end subroutine qinit
