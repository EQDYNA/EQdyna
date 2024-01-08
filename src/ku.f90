! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
subroutine ku
    use globalvar
    implicit none
    include 'mpif.h'
    
    integer (kind = 4) :: nel,m,i,j,ntemp,k,k1, non,itag,eqn,label
    real (kind = dp) :: det,constk,pstrinc,xc(3),matelement(5),esPML(21),es(12),ex(3,8),efPML(96),&
        elresf(nee), eleffm(nee), dl(ned,nen), vl(ned,nen), al(ned,nen)
        
    time1 = MPI_WTIME()
    
    do nel = 1, numel
        dl(1:ned,1:nen) = d(1:ned,ien(1:nen,nel))
        vl(1:ned,1:nen) = v(1:ned,ien(1:nen,nel))
        al(1:ned,1:nen) = 0.0d0
        al(1:ned,1:nen) = al(1:ned,1:nen) + rdampm*vl(1:ned,1:nen)
        
        do j = 1, nen
           do i = 1, ned
                if(i==3.and.C_elastic==0) then  !for inelastic off-fault, gravity included
                    al(i,j) = al(i,j) + grav*(roumax-(gamar+1.0d0)*rhow)/roumax
                endif
            enddo
        enddo

        elresf = 0.0d0 
        det = eledet(nel)
        eleffm(1:nee)=elemass(1:nee,nel)
        
        call contma(eleffm,al,elresf)

        if (et(nel)==1.or.et(nel)>10) then
            constk=-det
            es(1:12)    = s1(ids(nel)+1:ids(nel)+12)
            ex(1:3,1:8) = x(1:3,ien(1:8,nel))
            matelement(1:5) = mat(nel,1:5)
            call qdckd(eleshp(1,1,nel),matelement,vl,dl,es,elresf,constk,&
                    eleporep(nel),pstrinc,ex)

            pstrain(nel) = pstrain(nel) + pstrinc
            s1(ids(nel)+1:ids(nel)+12)=es(1:12)
            
            do i=1,nen					
                non=ien(i,nel)
                do j=1,ned
                    itag=locid(non)+j
                    k=id1(itag)					
                    if(k > 0) then
                        brhs(k) = brhs(k) + elresf((i-1)*ned+j)
                    endif
                enddo
            enddo
        elseif (et(nel)==2) then ! PML element. 
            efPML=0.0d0
            esPML=0.0d0 
            do i=1,8
                do j=1,3
                efPML((i-1)*12+9+j)=elresf((i-1)*3+j)
                enddo
            enddo
            do i=1,8
                do j=1,3
                    ex(j,i)=x(j,ien(i,nel))
                enddo
            enddo
            !do i=1,21
                esPML(1:21)=s1(ids(nel)+1:ids(nel)+21)
            !enddo
            !do i=1,5 
                matelement(1:5)=mat(nel,1:5)
            !enddo		
            call PMLwhg(vl,efPML,esPML,ex,matelement,eleshp(1,1,nel),det,nel)
            do i=1,8
                non=ien(i,nel)
                if (dof1(non)==12) then
                    do j=1,12
                        itag=locid(non)+j
                        eqn=id1(itag)
                        if (eqn.gt.0) then
                            brhs(eqn)=brhs(eqn)+efPML((i-1)*12+j)
                        endif
                    enddo
                elseif (dof1(non)==3) then
                    itag=locid(non)+1
                    eqn=id1(itag)
                    brhs(eqn)=brhs(eqn)+efPML((i-1)*12+1)+efPML((i-1)*12+2)+efPML((i-1)*12+3)+efPML((i-1)*12+10)
                    itag=locid(non)+2
                    eqn=id1(itag)
                    brhs(eqn)=brhs(eqn)+efPML((i-1)*12+4)+efPML((i-1)*12+5)+efPML((i-1)*12+6)+efPML((i-1)*12+11)
                    itag=locid(non)+3
                    eqn=id1(itag)
                    brhs(eqn)=brhs(eqn)+efPML((i-1)*12+7)+efPML((i-1)*12+8)+efPML((i-1)*12+9)+efPML((i-1)*12+12)				
                endif
            enddo
            !update the stress components 
            !do i=1,21
                s1(ids(nel)+1:ids(nel)+21)=esPML(1:21)
            !enddo
        endif!et(nel)=1/2
    enddo!nel loop
    timeused(4) = timeused(4) + MPI_WTIME() - time1
end subroutine ku
