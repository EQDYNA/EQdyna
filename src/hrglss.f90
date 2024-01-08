! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine hrglss
 
    use globalvar
    implicit none
    include 'mpif.h'

    integer (kind = 4) :: nel, i , j, k, m, itmp, non, itag, fi(4,8)
    real (kind = dp) :: phid(ned), dl(ned,nen), vl(ned,nen), fhr(ned,nen), f(24), det, coef, q(3,4)
    
    time1 = MPI_WTIME()
    do nel = 1, numel

        m = 1
        
        do i = 1,nen
            k = ien(i,nel)
            do j = 1,ned
                dl(j,i) = d(j,k)
                vl(j,i) = v(j,k)
                dl(j,i) = dl(j,i) + rdampk* vl(j,i)
            enddo		
        enddo
        if (C_hg == 1) then
            do itmp = 1, 4
                fhr = 0.0d0
                !... calculate sum(phi*dl)
                do i = 1,ned
                    phid(i) = 0.0d0
                    do j = 1,nen
                        phid(i) = phid(i) + phi(j,itmp,nel) * dl(i,j)
                    enddo
                enddo
                !... calculate hourglass resistence
                do i = 1,nen
                    fhr(1,i) = phi(i,itmp,nel) * (ss(1,nel)*phid(1)  &
                                + ss(2,nel)*phid(2) + ss(3,nel)*phid(3))
                    
                    fhr(2,i) = phi(i,itmp,nel) * (ss(2,nel)*phid(1)  &
                                + ss(4,nel)*phid(2) + ss(5,nel)*phid(3))
                    
                    fhr(3,i) = phi(i,itmp,nel) * (ss(3,nel)*phid(1)  &
                                + ss(5,nel)*phid(2) + ss(6,nel)*phid(3))
                enddo
                !... assemble to global right-hand-side force vector
                do i = 1,nen
                    do j = 1,ned						
                        non = ien(i,nel) 
                    
                        if (dof1(non) == 3) then
                            itag = locid(non)+j
                        elseif (dof1(non) == 12) then
                            itag = locid(non)+j+9
                        endif
                        
                        k = id1(itag)
                        
                        if(k > 0) then
                            brhs(k) = brhs(k) - fhr(j,i)
                        endif
                    enddo
                enddo
            enddo 

        elseif (C_hg==2) then
            !viscous hourglass control
            det=eledet(nel)
            coef=0.25d0*kapa_hg*mat(nel,3)*mat(nel,1)*(det*w)**(2.0d0/3.0d0)
            fi(1,1)=1;fi(1,2)=1;fi(1,3)=-1;fi(1,4)=-1;fi(1,5)=-1;fi(1,6)=-1;fi(1,7)=1;fi(1,8)=1
            fi(2,1)=1;fi(2,2)=-1;fi(2,3)=-1;fi(2,4)=1;fi(2,5)=-1;fi(2,6)=1;fi(2,7)=1;fi(2,8)=-1
            fi(3,1)=1;fi(3,2)=-1;fi(3,3)=1;fi(3,4)=-1;fi(3,5)=1;fi(3,6)=-1;fi(3,7)=1;fi(3,8)=-1
            fi(4,1)=1;fi(4,2)=-1;fi(4,3)=1;fi(4,4)=-1;fi(4,5)=-1;fi(4,6)=1;fi(4,7)=-1;fi(4,8)=1
            !fi1=[1,1,-1,-1,-1,-1,1,1]
            !fi2=[1,-1,-1,1,-1,1,1,-1]
            !fi3=[1,-1,1,-1,1,-1,1,-1]
            !fi4=[1,-1,1,-1,-1,1,-1,1]
            q=0.0d0!q(3,4)
            f=0.0d0
            do i=1,3
                do j=1,4
                    do k=1,8
                        q(i,j)=q(i,j)+vl(i,k)*fi(j,k)
                    enddo
                enddo
            enddo
            do k=1,8
                do i=1,3
                    do j=1,4
                        f((k-1)*3+i)=f((k-1)*3+i)-coef*q(i,j)*fi(j,k)
                    enddo
                enddo
            enddo	
            do i=1,nen
                do j=1,ned
                    non=ien(i,nel) 
                    if (dof1(non).eq.3) then
                        itag=locid(non)+j
                    elseif (dof1(non).eq.12) then
                        itag=locid(non)+j+9
                    endif
                    k=id1(itag)
                    if(k > 0) then
                        brhs(k) = brhs(k)+f((i-1)*3+j)
                    endif
                enddo
            enddo			
            !
        endif!C_HG(choose hg type)
    enddo
    timeused(5) = timeused(5) + MPI_WTIME() - time1
end subroutine hrglss