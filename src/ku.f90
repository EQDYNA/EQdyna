! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine ku
    use globalvar
    implicit none
    include 'mpif.h'
    
    integer (kind = 4) :: nel,m,i,j,ntemp,k,k1, non,itag,eqn,label
    real (kind = dp) :: pstrinc,efPML(96), elresf(nee), dl(ned,nen), vl(ned,nen), al(ned,nen)
        
    time1 = MPI_WTIME()
    
    do nel = 1, numel
        !dl(1:ned,1:nen) = d(1:ned,ien(1:nen,nel))
        !vl(1:ned,1:nen) = v(1:ned,ien(1:nen,nel))
        al(1:ned,1:nen) = 0.0d0
        al(1:ned,1:nen) = al(1:ned,1:nen) + rdampm*v(1:ned,ien(1:nen,nel))
        al(3,1:nen)     = al(3,1:nen) + (1.0d0-C_elastic)*grav*(roumax-(gamar+1.0d0)*rhow)/roumax

        elresf = 0.0d0 
        call contma(elemass(1:nee,nel),al,elresf)

        if (et(nel)==1 .or. et(nel)>10) then
            call qdckd(eleshp(1,1,nel), mat(nel,1:5), v(1:ned,ien(1:nen,nel)), &
                        d(1:ned,ien(1:nen,nel)), s1(ids(nel)+1:ids(nel)+12), &
                        elresf, -eledet(nel), eleporep(nel), pstrinc, &
                        x(1:3,ien(1:8,nel)))
            pstrain(nel) = pstrain(nel) + pstrinc
           
            do i=1,nen
                !non=ien(i,nel)
                do j=1,ned
                    !itag=locid(ien(i,nel))+j
                    k=id1(locid(ien(i,nel))+j)
                    if(k > 0) then
                        brhs(k) = brhs(k) + elresf((i-1)*ned+j)
                    endif
                enddo
            enddo
        elseif (et(nel)==2) then ! PML element. 
            efPML=0.0d0 
            do i=1,8
                do j=1,3
                    efPML((i-1)*12+9+j)=elresf((i-1)*3+j)
                enddo
            enddo

            call PMLwhg(v(1:ned,ien(1:nen,nel)), efPML, s1(ids(nel)+1:ids(nel)+21), &
                        x(1:3,ien(1:8,nel)), mat(nel,1:5), eleshp(1,1,nel), &
                        eledet(nel), nel)
            do i = 1, 8
                !non=ien(i,nel)
                if (dof1(ien(i,nel))==12) then
                    do j = 1, 12
                        !itag=locid(ien(i,nel))+j
                        eqn=id1(locid(ien(i,nel))+j)
                        if (eqn > 0) then
                            brhs(eqn)=brhs(eqn)+efPML((i-1)*12+j)
                        endif
                    enddo
                elseif (dof1(ien(i,nel)) == ndof) then
                    !itag=locid(ien(i,nel))+1
                    eqn=id1(locid(ien(i,nel))+1)
                    brhs(eqn)=brhs(eqn)+efPML((i-1)*12+1)+efPML((i-1)*12+2)+efPML((i-1)*12+3)+efPML((i-1)*12+10)
                    !itag=locid(non)+2
                    eqn=id1(locid(ien(i,nel))+2)
                    brhs(eqn)=brhs(eqn)+efPML((i-1)*12+4)+efPML((i-1)*12+5)+efPML((i-1)*12+6)+efPML((i-1)*12+11)
                    !itag=locid(non)+3
                    eqn=id1(locid(ien(i,nel))+3)
                    brhs(eqn)=brhs(eqn)+efPML((i-1)*12+7)+efPML((i-1)*12+8)+efPML((i-1)*12+9)+efPML((i-1)*12+12)				
                endif
            enddo
        endif
    enddo
    timeused(4) = timeused(4) + MPI_WTIME() - time1
end subroutine ku
