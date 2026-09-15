subroutine checkInputConsistency

    use globalvar
    use errorCodes
    implicit none

    if (C_elastic==0.and.C_Q==1) then
        call abortRun(ERR_CFG_Q_NEEDS_ELASTIC, &
            'Q model (C_Q=1) can only work with the elastic code (C_elastic=1). Set C_Q=0 or C_elastic=1.')
    endif
    if (C_Q==1.and.rat>1) then
        call abortRun(ERR_CFG_Q_NEEDS_UNIFORM, &
            'Q model (C_Q=1) can only work with uniform element size; rat must be 1.0.')
    endif
    if (output_plastic == 1 .and. C_elastic/=0) then
        write(*,*) 'Now, C_elastic = ', C_elastic
        call abortRun(ERR_CFG_PLASTIC_OUTPUT, &
            'Plastic strains are only output for C_elastic=0. Set output_plastic=0 or C_elastic=0.')
    endif
end subroutine checkInputConsistency
