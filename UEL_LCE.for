      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=4000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=20000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global


***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      USE GLOBAL
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      uvar(1) = globalSdv(noel-ElemOffset,npt,1)

      return
      END SUBROUTINE uvarm
      
 
      !==================== UEL SUBROUTINE ====================     
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     &     PROPS,NPROPS,coords,MCRD,nNODE,Uall,DUall,Vel,Accn,JTYPE,
     &     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &     NJPROP,PERIOD)   
      
      USE GLOBAL
      
      IMPLICIT NONE
       
!     VARIABLES DEFINED IN UEL,PASSED BACK TO ABAQUS
      REAL(KIND=8) :: RHS(MLVARX,*),AMATRX(nDofEl,nDofEl),SVARS(*),
     &   ENERGY(8),PNEWDT
!     VARIABLES PASSED INTO UEL 
      INTEGER :: nDofEl,NRHS,NSVARS,NPROPS,MCRD,nNODE,JTYPE,
     &   KSTEP,KINC
      INTEGER :: JELEM,NDLOAD,NPREDF,MLVARX,MDLOAD,NJPROP
      INTEGER :: JDLTYP(MDLOAD,*),LFLAGS(*),JPROPS(*)
      REAL(KIND=8) :: PROPS(*),coords(MCRD,nNODE),UALL(nDofEl),
     &   DUALL(MLVARX,*)
      REAL(KIND=8) :: VEL(nDofEl),ACCN(nDofEl),TIME(2),DTIME,
     &   PARAMS(*)
      REAL(KIND=8) :: ADLMAG(MDLOAD,*),PREDEF(2,NPREDF,nNODE),
     &   DDLMAG(MDLOAD,*),PERIOD
      
      integer lenJobName,lenOutDir,nDim,nInt,nIntS
      character*256 jobName,outDir,fileName

      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     &     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')
      

      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (64,65,72,73)
      !
      if((lflags(1).eq.1).or.(lflags(1).eq.2)) then
         !
         ! all is good
         !
      else
        write(*,*) 'Abaqus does not have the right procedure'
        write(*,*) 'go back and chekc the procedure type'
        write(*,*) 'lflags(1)=',lflags(1)
        write(80,*) 'Abaqus does not have the right procedure'
        write(80,*) 'go back and chekc the procedure type'
        write(80,*) 'lflags(1)=',lflags(1)
        call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      
      if(lflags(2).eq.0) then
        !
        ! lflags(2)=0 -> small disp.
        ! lflags(2)=1 -> large disp.
        !
        write(*,*) 'Abaqus thinks you are doing'
        write(*,*) 'a small displacement analysis'
        write(*,*) 'go in and set nlgeom=yes'
        write(80,*) 'Abaqus thinks you are doing'
        write(80,*) 'a small displacement analysis'
        write(80,*) 'go in and set nlgeom=yes'
        call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
        !
        ! lflags(4)=0 -> general step
        ! lflags(4)=1 -> linear purturbation step
        !
        write(*,*) 'Abaqus thinks you are doing'
        write(*,*) 'a linear purturbation step'
        write(80,*) 'Abaqus thinks you are doing'
        write(80,*) 'a linear purturbation step'
        call xit         
      endif
      !
      ! DO nothing IF a ``dummy'' step
      !
      IF(dtime.EQ.0.0) RETURN
      !----------------------------------------------------------------
      ! 
      ! CALL the paricular element to perform the analysis
      !
      IF(jtype.EQ.1) then
         !
         ! This is a plane strain analysis
         !
         nINT=4
         nINTS=1
         nDIM=2
         CALL UPE4(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     &        PROPS,NPROPS,coords,MCRD,nNODE,Uall,DUall,Vel,Accn,JTYPE,
     &        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &        NJPROP,PERIOD,
     &        nDIM,nINT,nINTS)
         !
         !
         ELSEIF(jtype.EQ.2) then
         !
         ! This is an axisymmetric analysis
         !
         nINT=4
         nINTS=1
         nDIM=2
         CALL UAX4(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     &        PROPS,NPROPS,coords,MCRD,nNODE,Uall,DUall,Vel,Accn,JTYPE,
     &        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &        NJPROP,PERIOD,
     &        nDIM,nINT,nINTS)
         !
         !
      ELSEIF(jtype.EQ.3) then
         !
         ! This is a 3D analysis
         !
         nINT=8
         nINTS=1
         nDIM=3
         CALL U3D8(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,nNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDIM,nINT,nINTS)
         !
         !
      ELSE
         !
         ! We have a problem...
         !
         print *, 'Element type not supported, jtype=',jtype
         call xit
         !
      ENDIF
      !
      ! DOne with this element, RHS and AMATRX already returned
      !  as output from the specIFic element routine CALLed
      !
      !----------------------------------------------------------------

      RETURN
      END SUBROUTINE UEL
      
!==================== UPE4 SUBROUTINE ====================     
      SUBROUTINE UPE4(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     &    PROPS,NPROPS,coords,MCRD,nNODE,UALL,DUALL,VEL,ACCN,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,
     &    JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,
     &    MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD,
     &    nDIM,nINT,nINTS)   
      
      USE GLOBAL
      
      IMPLICIT NONE
      !
      ! variables defined in uel,passed back to Abaqus
      !
      REAL(KIND=8) :: RHS(MLVARX,*),AMATRX(nDofEl,nDofEl),SVARS(*),
     &   ENERGY(8),PNEWDT
      !
      ! variables passed into UEL
      !
      INTEGER :: nDofEl,NRHS,NSVARS,NPROPS,MCRD,nNODE,JTYPE,
     &   KSTEP,KINC
      INTEGER :: JELEM,NDLOAD,NPREDF,MLVARX,MDLOAD,NJPROP
      INTEGER :: JDLTYP(MDLOAD,*),LFLAGS(*),JPROPS(*)
      REAL(KIND=8) :: PROPS(*),coords(MCRD,nNODE),UALL(nDofEl),
     &   DUALL(MLVARX,*)
      REAL(KIND=8) :: VEL(nDofEl),ACCN(nDofEl),TIME(2),DTIME,
     &   PARAMS(*)
      REAL(KIND=8) :: ADLMAG(MDLOAD,*),PREDEF(2,NPREDF,nNODE),
     &   DDLMAG(MDLOAD,*),PERIOD
      !
      ! variables defined and used in the UEL
      !
      INTEGER :: INDI,INDJ,INDK,INDJJ,INDA,ROW,COL,INDB,INDL
      INTEGER :: INTPT,STAT
      REAL(KIND=8),PARAMETER :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0
      REAL(KIND=8),PARAMETER :: HALF=0.5D0,PI=3.141592653D0
      REAL(KIND=8),PARAMETER :: THREE=3.0D0,THIRD=1.0D0/3.0D0
      INTEGER :: nDIM 
      INTEGER :: nSDV,nlSDV,ngSDV 
      INTEGER :: nINT, nINTS 
      REAL(KIND=8) :: U(nNODE,nDIM),DU(nNODE,nDIM),
     &   Uold(nNODE,nDofEl),BODY(3)
      REAL(KIND=8) :: coordsC(MCRD,nNODE),RU(2*nNODE,1),
     &   KUU(2*nNODE,2*nNODE)
      REAL(KIND=8) :: KROND(3,3),XI(nINT,nDIM),W(nINT),
     &   SH0(nNODE),SH(nNODE)
      REAL(KIND=8) :: DSH0(nNODE,nDIM),DSHC0(nNODE,nDIM),
     &   DSH(nNODE,nDIM),DSHC(nNODE,nDIM),DSHXI(nNODE,nDIM)
      REAL(KIND=8) :: DetMAPJ0,DetMAPJ0C,DetMAPJ,DetMAPJC,
     &   Fc_Tau(3,3),DetFc_Tau,Fc_T(3,3),DetFc_T
      REAL(KIND=8) :: D_T(3),D0(3),D_Tau(3)
      REAL(KIND=8) :: F_Tau(3,3),F_T(3,3),DetF_Tau,
     &   T_Tau(3,3),FV_T(3,3),FV_Tau(3,3),DetF_T,DetF
      REAL(KIND=8) :: Q0,TangentStiffness(3,3,3,3),PKSTRESS(3,3)
      REAL(KIND=8) :: ZERO_MATRIX(1,2),Le      
      !
      ! Get element parameters
      !
      nlSDV = JPROPS(1) ! number of local sdv's per integ point
      ngSDV = JPROPS(2) ! number of global sdv's per integ point
      
      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
        !
        ! allocate memory for the globalSdv's
        !
        ! numElem needs to be set in the MODULE
        ! nInt needs to be set in the UEL
        !
        stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
        allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
        if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
        endif
        write(*,*) '-------------------------------------------------'
        write(*,*) '----------- globalSDV ALLOCATED -----------------'
        write(*,*) '-------------------------------------------------'
        write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
        write(*,*) '---------- numElem=',numElem
        write(*,*) '---------- UPE4 ELEMENTS ------------------------'
        write(*,*) '-------------------------------------------------'
        write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
        write(*,*) '---------- nInt =',nInt
        write(*,*) '---------- nIntS=',nIntS
        write(*,*) '-------------------------------------------------'
        write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
        write(*,*) '---------- ngSdv=',ngSdv
        write(*,*) '-------------------------------------------------'
      endif
      
      
      ! obtain initial conditions
      Q0=PROPS(9)
      D0(1)=COS(Q0)
      D0(2)=SIN(Q0)
      D0(3)=0   

      ! Initialize the residual and tangent matrices to zero
      RU=ZERO
      KUU=ZERO
      ENERGY=ZERO
      BODY=ZERO
      
      ! Identity tensor
      CALL CREATE_IDENTITY_MATRIX(3,KROND)

      ! Obtain nodal displacements
      INDK=0
      DO INDI=1,nNODE
        DO INDJ=1,nDIM
        INDK=INDK + 1
        U(INDI,INDJ)=UALL(INDK)
        DU(INDI,INDJ)=DUALL(INDK,1)
        Uold(INDI,INDJ)=U(INDI,INDJ) - DU(INDI,INDJ)
        ENDDO
      ENDDO

      ! Obtain current nodal coordinates
      DO INDI=1,nNODE
        DO INDJ=1,nDIM
            coordsC(INDJ,INDI)=coords(INDJ,INDI) + U(INDI,INDJ)
        ENDDO
      ENDDO
      
      ! displacement increment, based on element diagonal
      Le=dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     &     ((coordsC(2,1)-coordsC(2,3))**two))
      
      DO INDI=1,nNODE
        DO INDJ=1,nDIM
          IF(dabs(DU(INDI,INDJ)).GT.10.0*Le) THEN
            PNEWDT=0.5
            RETURN
          ENDIF
         ENDDO
      ENDDO
      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto,E.A.,Peric,D.,Dutko,M.,Owen,D.R.J.,1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33,3277-3296.
      !
      !
      ZERO_MATRIX=0.0D0
      IF (nNODE.EQ.4) THEN
          CALL calcShape2DLinear(1,ZERO_MATRIX,1,SH0,DSHXI)
      ELSE
          print *, 'Incorrect number of nodes: nNODE.NE.4'
          call xit
      ENDIF

      ! Map shape functions from local to global reference coordinate system
      CALL MAPShape2D(nNODE,DSHXI,coords,DSH0,DetMAPJ0,STAT)
      IF (STAT.EQ.0) THEN
          PNEWDT=0.5
          RETURN
      ENDIF     
      
      ! Map shape functions from local to global current coordinate system
      CALL MAPShape2D(nNODE,DSHXI,coordsC,DSHC0,DetMAPJ0C,STAT)
      IF (STAT.EQ.0) THEN
          PNEWDT=0.5
          RETURN
      ENDIF
      

      ! Calculate the deformation gradient at the element centroid
      ! at the end of the increment for use in the `F-bar' method
      ! The subscript tau denotes the time at the end of the increment.

      Fc_Tau=KROND
      Fc_T=KROND
      DO INDI=1,nDIM
        DO INDJ=1,nDIM
          DO INDK=1,nNODE
              Fc_Tau(INDI,INDJ)=Fc_Tau(INDI,INDJ) + 
     &            DSH0(INDK,INDJ) * U(INDK,INDI)
              Fc_T(INDI,INDJ)=Fc_T(INDI,INDJ) +
     &            DSH0(INDK, INDJ)*Uold(INDK,INDI)
          ENDDO
        ENDDO
      ENDDO
      !
      ! Modify for plane-strain
      !
      Fc_Tau(3,3)=ONE
      Fc_T(3,3)=ONE      
      !
      ! 2D plane-strain implementation DetF
      !
      DetFc_T=Fc_T(1,1)*Fc_T(2,2) - Fc_T(1,2)*Fc_T(2,1)
      DetFc_Tau=Fc_Tau(1,1)*Fc_Tau(2,2) - Fc_Tau(1,2)*Fc_Tau(2,1)
      
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      ! Begin the loop over integration points
       
      ! Obtain integration point local coordinates and weights
      IF(nNODE.EQ.4) then
         !
         ! gauss integration for a rectangular element
         !
         IF(nINT.EQ.4) then    
            CALL xint2D4pt(xi,w,nINT) ! 4-pt integration, nINT=4 above
         ELSEIF(nINT.EQ.1) then
            CALL xint2D1pt(xi,w,nINT) ! 1-pt integration, nINT=1 above
         ELSE
            print *, 'Invalid number of int points, nINT=',nINT
            call xit
         ENDIF
      ELSE
         print *, 'Incorrect number of nodes: nNODE.ne.4'
         call xit
      ENDIF
      
      
      ! Loop over integration points
      INDJJ=0  ! JJ is used for tracking the state variables

      DO INTPT=1,nINT
        ! Obtain state variables from previous increment
        IF ((KINC.LE.1) .AND. (KSTEP.EQ.1)) THEN
            ! This is the first increment of the first step.
            ! Give initial conditions.
            FV_T=KROND
            F_T=KROND
            D_T=D0
        ELSE
            ! This is not the first increment; read old values.
            FV_T(1,1)=SVARS(1 + INDJJ)
            FV_T(2,2)=SVARS(2 + INDJJ)
            FV_T(3,3)=SVARS(3 + INDJJ)
            FV_T(2,3)=SVARS(4 + INDJJ)
            FV_T(3,1)=SVARS(5 + INDJJ)
            FV_T(1,2)=SVARS(6 + INDJJ)
            FV_T(3,2)=SVARS(7 + INDJJ)
            FV_T(1,3)=SVARS(8 + INDJJ)
            FV_T(2,1)=SVARS(9 + INDJJ)
            F_T(1,1)=SVARS(10 + INDJJ)
            F_T(2,2)=SVARS(11 + INDJJ)
            F_T(3,3)=SVARS(12 + INDJJ)
            F_T(2,3)=SVARS(13 + INDJJ)
            F_T(3,1)=SVARS(14 + INDJJ)
            F_T(1,2)=SVARS(15 + INDJJ)
            F_T(3,2)=SVARS(16 + INDJJ)
            F_T(1,3)=SVARS(17 + INDJJ)
            F_T(2,1)=SVARS(18 + INDJJ)
            D_T(1)=SVARS(19 + INDJJ)
            D_T(2)=SVARS(20 + INDJJ)
            D_T(3)=SVARS(21 + INDJJ)
        END IF
    
        ! Obtain shape functions and their local gradients
        IF (nNODE.EQ.4) THEN
            CALL calcShape2DLinear(nINT,XI,INTPT,SH,DSHXI)
        ELSE
            print *, 'Incorrect number of nodes: nNODE /= 4'
            call xit
        ENDIF

        ! Map shape functions from local to global reference coordinate system
        CALL MAPShape2D(nNODE,DSHXI,coords,DSH,DetMAPJ,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Map shape functions from local to global current coordinate system
        CALL MAPShape2D(nNODE,DSHXI,coordsC,DSHC,DetMAPJC,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Obtain the deformation gradient at this integration point.
        ! The subscript tau denotes the time at the end of the increment.

        F_Tau=KROND
        F_t=KROND
        DO INDI=1,nDIM
            DO INDJ=1,nDIM
                DO INDK=1,nNODE
                    F_Tau(INDI,INDJ)=F_Tau(INDI,INDJ) +
     &               DSH(INDK,INDJ) * U(INDK,INDI)
                    F_T(INDI,INDJ)=F_T(INDI,INDJ) +
     &               DSH(INDK,INDJ) * Uold(INDK,INDI)
                ENDDO
            ENDDO
        ENDDO
        !
        ! Modify F(3,3) for plane-strain 
        F_Tau(3,3)=one
        F_T(3,3)=one
        
        !
        ! Modify the deformation gradient for the `F-bar' method
        !  only when using the 4 node fully integrated linear
        !  element, do not use the `F-bar' method for any other element
        !
        if((nNode.EQ.4).and.(nInt.EQ.4)) then
            !
            !  2D plane-strain implementation
            !
            DetF_t = F_t(1,1)*F_t(2,2) - F_t(1,2)*F_t(2,1)
            DetF_Tau = F_Tau(1,1)*F_Tau(2,2) - F_Tau(1,2)*F_Tau(2,1)
            do INDI=1,nDim
                do INDJ=1,nDim
                  F_Tau(INDI,INDJ)=((DetFc_tau/DetF_Tau)**half)*
     &        F_Tau(INDI,INDJ)
                  F_t(INDI,INDJ) = ((DetFc_t/DetF_t)**half)*
     &        F_t(INDI,INDJ)
                enddo
            enddo
        endif
        call mDet(F_Tau,DetF_Tau)
        !
        ! Perfom the constitutive time integration at integration point
        !
        CALL MATERIAL(PROPS,F_T,FV_T,D_T,F_Tau,FV_Tau,D_Tau,
     &    DTIME,D0,PKSTRESS,TangentStiffness,STAT)   
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Save the state variables at this integ point
        !  at the end of the increment.
        !
        SVARS(1 + INDJJ)=FV_Tau(1,1)
        SVARS(2 + INDJJ)=FV_Tau(2,2)
        SVARS(3 + INDJJ)=FV_Tau(3,3)
        SVARS(4 + INDJJ)=FV_Tau(2,3)
        SVARS(5 + INDJJ)=FV_Tau(3,1)
        SVARS(6 + INDJJ)=FV_Tau(1,2)
        SVARS(7 + INDJJ)=FV_Tau(3,2)
        SVARS(8 + INDJJ)=FV_Tau(1,3)
        SVARS(9 + INDJJ)=FV_Tau(2,1)
        SVARS(10 + INDJJ)=F_Tau(1,1)
        SVARS(11 + INDJJ)=F_Tau(2,2)
        SVARS(12 + INDJJ)=F_Tau(3,3)
        SVARS(13 + INDJJ)=F_Tau(2,3)
        SVARS(14 + INDJJ)=F_Tau(3,1)
        SVARS(15 + INDJJ)=F_Tau(1,2)
        SVARS(16 + INDJJ)=F_Tau(3,2)
        SVARS(17 + INDJJ)=F_Tau(1,3)
        SVARS(18 + INDJJ)=F_Tau(2,1)
        SVARS(19 + INDJJ)=D_Tau(1)
        SVARS(20 + INDJJ)=D_Tau(2)
        SVARS(21 + INDJJ)=D_Tau(3)
        INDJJ=INDJJ + nlSDV
             
        ! Save the state variables at this integ point in the
        !  global array used for plotting field output
        !
        globalSdv(jelem,intPt,1) = D_Tau(1)   ! polymer volume fraction
        
        ! Compute/update the displacement residual vector
        DO INDA=1,nNODE
            DO INDI=1,nDIM
            ROW=nDIM * (INDA - 1) + INDI
            DO INDJ=1,nDIM
            RU(ROW,1)=RU(ROW,1) - PKSTRESS(INDI,INDJ) *
     &         DSH(INDA,INDJ) * W(INTPT) * DetMAPJ
            ENDDO
            ENDDO
        ENDDO

        ! Compute/update the displacement tangent matrix
        DO INDA=1,nNODE
            DO INDI=1,nDIM
                DO INDB=1,nNODE
                    DO INDK=1,nDIM
                        ROW=nDIM * (INDA - 1) + INDI
                        COL=nDIM * (INDB - 1) + INDK
                        DO INDJ=1,nDIM
                            DO INDL=1,nDIM
                            KUU(ROW,COL)=KUU(ROW,COL) + 
     &      TangentStiffness(INDI,INDJ,INDK,INDL) * DSH(INDB,INDL) 
     &           * DSH(INDA,INDJ) * W(INTPT) * DetMAPJ
                            ENDDO
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
      ENDDO
            !
            ! End the loop over integration points
            !----------------------------------------------------------------
            ! RETURN to Abaqus the RHS vector and the StIFfness matrix.  This
            !  is essentially giving Abaqus the residual and the tangent matrix.
            !
            ! RETURN to Abaqus the right hand side vector.
            AMATRX=KUU
            DO INDI=1,nNODE
                RHS(nDIM * (INDI - 1) + 1,1)=RU(nDIM * (INDI - 1) + 1,1)
                RHS(nDIM * (INDI - 1) + 2,1)=RU(nDIM * (INDI - 1) + 2,1)
                RHS(nDIM * (INDI - 1) + 3,1)=RU(nDIM * (INDI - 1) + 3,1)
            ENDDO
      
      END SUBROUTINE UPE4

      
!==================== UAX4 SUBROUTINE ====================     
      SUBROUTINE UAX4(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     &    PROPS,NPROPS,coords,MCRD,nNODE,UALL,DUALL,VEL,ACCN,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,
     &    JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,
     &    MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD,
     &    nDIM,nINT,nINTS)   
      
        IMPLICIT NONE
        !
        ! variables defined in uel,passed back to Abaqus
        !
        REAL(KIND=8) :: RHS(MLVARX,*),AMATRX(nDofEl,nDofEl),SVARS(*),
     &   ENERGY(8),PNEWDT
        !
        ! variables passed into UEL
        !
        INTEGER :: nDofEl,NRHS,NSVARS,NPROPS,MCRD,nNODE,JTYPE,
     &   KSTEP,KINC
        INTEGER :: JELEM,NDLOAD,NPREDF,MLVARX,MDLOAD,NJPROP
        INTEGER :: JDLTYP(MDLOAD,*),LFLAGS(*),JPROPS(*)
        REAL(KIND=8) :: PROPS(*),coords(MCRD,nNODE),UALL(nDofEl),
     &   DUALL(MLVARX,*)
        REAL(KIND=8) :: VEL(nDofEl),ACCN(nDofEl),TIME(2),DTIME,
     &   PARAMS(*)
        REAL(KIND=8) :: ADLMAG(MDLOAD,*),PREDEF(2,NPREDF,nNODE),
     &   DDLMAG(MDLOAD,*),PERIOD
        !
        ! variables defined and used in the UEL
        !
        INTEGER :: INDI,INDJ,INDK,INDJJ,INDA,ROW,COL,INDB,INDL
        INTEGER :: INTPT,STAT
        INTEGER :: nINTSURF,INTPTSURF
        REAL(KIND=8),PARAMETER :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0
        REAL(KIND=8),PARAMETER :: HALF=0.5D0,PI=3.141592653D0
        REAL(KIND=8),PARAMETER :: THREE=3.0D0,THIRD=1.0D0 / 3.0D0
        INTEGER :: nDOF
        INTEGER :: nDIM 
        INTEGER :: nSDV,nlSDV,ngSDV
        INTEGER :: nINT,nINTS  
        REAL(KIND=8) :: U(nNODE,3),DU(nNODE,nDIM),
     &   Uold(nNODE,nDofEl)
        REAL(KIND=8) :: coordsC(MCRD,nNODE),RU(3*nNODE,1),
     &   KUU(3*nNODE,3*nNODE)
        REAL(KIND=8) :: KROND(3,3),XI(nINT,3),W(nINT),
     &   SH0(nNODE),SH(nNODE)
        REAL(KIND=8) :: DSH0(nNODE,3),DSHC0(nNODE,3),DSH(nNODE,3),
     &   DSHC(nNODE,3),DSHXI(nNODE,3)
        REAL(KIND=8) :: DetMAPJ0,DetMAPJ0C,DetMAPJ,DetMAPJC,
     &   Fc_Tau(3,3),DetFc_Tau
        REAL(KIND=8) :: D_T(3),D0(3),D_Tau(3)
        REAL(KIND=8) :: F_Tau(3,3),F_T(3,3),DetF_Tau,
     &   T_Tau(3,3),FV_T(3,3),FV_Tau(3,3)
        REAL(KIND=8) :: Q0,TangentStiffness(3,3,3,3),PKSTRESS(3,3)
        REAL(KIND=8) :: ZERO_MATRIX(1,3)

        REAL(KIND=8) :: SMAT(6,1),BMAT(6,3*nNODE),GMAT(9,3*nNODE)
        REAL(KIND=8) :: G0MAT(9,3*nNODE),AMAT(9,9),QMAT(9,9),
     &   BODY(3)
        REAL(KIND=8) :: BODYFORCERES(3*nNODE,1),GAMMA
        ! Check the procedure type; this should be a 
        !  *Static step,which is either 1 or 2
        !
        IF ((LFLAGS(1).EQ.1) .OR. (LFLAGS(1).EQ.2)) THEN
        ! Correct procedure specIFied
        ELSE
        print *, 'Abaqus DOes not have the right procedure'
        print *, 'Go back and check the procedure type'
        print *, 'LFLAGS(1)=',LFLAGS(1)
        call xit
        ENDIF
        ! Make sure Abaqus knows you are DOing a large
        !  deformation problem
        !
        IF (LFLAGS(2).EQ.0) THEN
        ! LFLAGS(2)=0 -> Small displacement
        ! LFLAGS(2)=1 -> Large displacement
        print *, 'Abaqus thinks you are DOing'
        print *, 'a small displacement analysis'
        print *, 'Go in and set nlgeom=yes'
        call xit
        ENDIF
        ! Check to see IF you are DOing a general
        !  step or a linear perturbation step
        !
        IF (LFLAGS(4).EQ.1) THEN
        ! LFLAGS(4)=0 -> General step
        ! LFLAGS(4)=1 -> Linear perturbation step
        print *, 'Abaqus thinks you are DOing'
        print *, 'a linear perturbation step'
        call xit
        ENDIF
        ! DO nothing IF a ``dummy'' step
        !
        IF(dtime.EQ.zero) RETURN
      
        ! Initialize the residual and tangent matrices to zero
        RU=0.0D0
        KUU=0.0D0
        ENERGY=0.0D0
      
        CALL CREATE_IDENTITY_MATRIX(3,KROND)

        ! Obtain nodal displacements
        INDK=0
        DO INDI=1,nNODE
            DO INDJ=1,nDIM
                INDK=INDK + 1
                U(INDI,INDJ)=UALL(INDK)
                DU(INDI,INDJ)=DUALL(INDK,1)
                Uold(INDI,INDJ)=U(INDI,INDJ) - DU(INDI,INDJ)
            ENDDO
        ENDDO
        ! Obtain current nodal coordinates
        DO INDI=1,nNODE
            DO INDJ=1,nDIM
                coordsC(INDJ,INDI)=coords(INDJ,INDI) + U(INDI,INDJ)
            ENDDO
        ENDDO
        !----------------------------------------------------------------
        ! 
        ! Take this opportunity to perform calculations at the element
        !  centroid.  Get the deformation gradient for use in the `F-bar' method.
        !
        ! Reference for the F-bar method:
        !  de Souza Neto,E.A.,Peric,D.,Dutko,M.,Owen,D.R.J.,1996.
        !  Design of simple low order finite elements for large strain
        !  analysis of nearly incompressible solids. International Journal
        !  of Solids and Structures 33,3277-3296.
        !
        ! Obtain shape functions and their local gradients at the element
        !  centroid,that means xi_1=xi_2=x_3=0.0,and nINT=1
        !
        ZERO_MATRIX=0.0D0
        IF (nNODE.EQ.8) THEN
            CALL CALCSHAPE3DLINEAR(1,ZERO_MATRIX,1,SH0,DSHXI)
        ELSE
            print *, 'Incorrect number of nodes: nNODE /= 8'
            call xit
        ENDIF

        ! Map shape functions from local to global reference coordinate system
        CALL MAPSHAPE3D(nNODE,DSHXI,coords,DSH0,DetMAPJ0,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF     
        ! Map shape functions from local to global current coordinate system
        CALL MAPSHAPE3D(nNODE,DSHXI,coordsC,DSHC0,DetMAPJ0C,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Calculate the deformation gradient at the element centroid
        ! at the end of the increment for use in the `F-bar' method
        ! The subscript tau denotes the time at the end of the increment.

        Fc_Tau=KROND
        DO INDI=1,nDIM
            DO INDJ=1,nDIM
                DO INDK=1,nNODE
                    Fc_Tau(INDI,INDJ)=Fc_Tau(INDI,INDJ) + 
     &               DSH0(INDK,INDJ) * U(INDK,INDI)
                ENDDO
            ENDDO
        ENDDO
        CALL MDet(Fc_Tau,DetFc_Tau)

        !----------------------------------------------------------------
        ! Begin the loop over integration points

        ! Obtain integration point local coordinates and weights

        IF (nNODE.EQ.8) THEN
            IF (nINT.EQ.8) THEN
                CALL XINT3D8PT(XI,W,nINT)  ! 8-pt integration,nINT=8 above
            ELSE
                print *, 'Incorrect number of int points: nINT /= 8'
                call xit
            END IF
        ELSE
            print *, 'Incorrect number of nodes: nNODE /= 8'
            call xit
        ENDIF

        ! Loop over integration points
        Q0=PROPS(9)
        D0(1)=COS(Q0)
        D0(2)=SIN(Q0)
        D0(3)=0
        INDJJ=0  ! JJ is used for tracking the state variables
        DO INTPT=1,nINT
            ! Obtain state variables from previous increment
            IF ((KINC.LE.1) .AND. (KSTEP.EQ.1)) THEN
                ! This is the first increment of the first step.
                ! Give initial conditions.
                FV_T=KROND
                F_T=KROND
                D_T=D0
            ELSE
                ! This is not the first increment; read old values.
                FV_T(1,1)=SVARS(1 + INDJJ)
                FV_T(2,2)=SVARS(2 + INDJJ)
                FV_T(3,3)=SVARS(3 + INDJJ)
                FV_T(2,3)=SVARS(4 + INDJJ)
                FV_T(3,1)=SVARS(5 + INDJJ)
                FV_T(1,2)=SVARS(6 + INDJJ)
                FV_T(3,2)=SVARS(7 + INDJJ)
                FV_T(1,3)=SVARS(8 + INDJJ)
                FV_T(2,1)=SVARS(9 + INDJJ)
                F_T(1,1)=SVARS(10 + INDJJ)
                F_T(2,2)=SVARS(11 + INDJJ)
                F_T(3,3)=SVARS(12 + INDJJ)
                F_T(2,3)=SVARS(13 + INDJJ)
                F_T(3,1)=SVARS(14 + INDJJ)
                F_T(1,2)=SVARS(15 + INDJJ)
                F_T(3,2)=SVARS(16 + INDJJ)
                F_T(1,3)=SVARS(17 + INDJJ)
                F_T(2,1)=SVARS(18 + INDJJ)
                D_T(1)=SVARS(19 + INDJJ)
                D_T(2)=SVARS(20 + INDJJ)
                D_T(3)=SVARS(21 + INDJJ)
            END IF
    
        ! Obtain shape functions and their local gradients
        IF (nNODE.EQ.8) THEN
            CALL CALCSHAPE3DLINEAR(nINT,XI,INTPT,SH,DSHXI)
        ELSE
            print *, 'Incorrect number of nodes: nNODE /= 8'
            call xit
        ENDIF

        ! Map shape functions from local to global reference coordinate system
        CALL MAPSHAPE3D(nNODE,DSHXI,coords,DSH,DetMAPJ,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Map shape functions from local to global current coordinate system
        CALL MAPSHAPE3D(nNODE,DSHXI,coordsC,DSHC,DetMAPJC,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Obtain the deformation gradient at this integration point.
        ! The subscript tau denotes the time at the end of the increment.

        F_Tau=KROND
        DO INDI=1,nDIM
            DO INDJ=1,nDIM
                DO INDK=1,nNODE
                    F_Tau(INDI,INDJ)=F_Tau(INDI,INDJ) +
     &               DSH(INDK,INDJ) * U(INDK,INDI)
                ENDDO
            ENDDO
        ENDDO

        ! Modify the deformation gradient for the `F-bar' method.

        !CALL MDet(F_Tau,DetF_Tau)
        !F_Tau=((DetFCTAU / DetF_Tau) ** THIRD) * F_Tau

        !
        ! Perform the constitutive update at this integ. point
        !
        CALL MATERIAL(PROPS,F_T,FV_T,D_T,F_Tau,FV_Tau,D_Tau,
     &    DTIME,D0,PKSTRESS,TangentStiffness,STAT)   
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF

        ! Save the state variables at this integ point
        !  at the end of the increment.
        !
        SVARS(1 + INDJJ)=FV_Tau(1,1)
        SVARS(2 + INDJJ)=FV_Tau(2,2)
        SVARS(3 + INDJJ)=FV_Tau(3,3)
        SVARS(4 + INDJJ)=FV_Tau(2,3)
        SVARS(5 + INDJJ)=FV_Tau(3,1)
        SVARS(6 + INDJJ)=FV_Tau(1,2)
        SVARS(7 + INDJJ)=FV_Tau(3,2)
        SVARS(8 + INDJJ)=FV_Tau(1,3)
        SVARS(9 + INDJJ)=FV_Tau(2,1)
        SVARS(10 + INDJJ)=F_Tau(1,1)
        SVARS(11 + INDJJ)=F_Tau(2,2)
        SVARS(12 + INDJJ)=F_Tau(3,3)
        SVARS(13 + INDJJ)=F_Tau(2,3)
        SVARS(14 + INDJJ)=F_Tau(3,1)
        SVARS(15 + INDJJ)=F_Tau(1,2)
        SVARS(16 + INDJJ)=F_Tau(3,2)
        SVARS(17 + INDJJ)=F_Tau(1,3)
        SVARS(18 + INDJJ)=F_Tau(2,1)
        SVARS(19 + INDJJ)=D_Tau(1)
        SVARS(20 + INDJJ)=D_Tau(2)
        SVARS(21 + INDJJ)=D_Tau(3)
        INDJJ=INDJJ + nSDV
        ! Compute/update the displacement residual vector
        DO INDA=1,nNODE
        DO INDI=1,nDOF
        ROW=nDOF * (INDA - 1) + INDI
        DO INDJ=1,nDIM
            RU(ROW,1)=RU(ROW,1) - PKSTRESS(INDI,INDJ) *
     &         DSH(INDA,INDJ) * W(INTPT) * DetMAPJ
        ENDDO
        ENDDO
        ENDDO
        ! Compute/update the displacement tangent matrix
        DO INDA=1,nNODE
          DO INDI=1,nDOF
            DO INDB=1,nNODE
              DO INDK=1,nDOF
                ROW=nDOF * (INDA - 1) + INDI
                COL=nDOF * (INDB - 1) + INDK
                DO INDJ=1,nDIM
                DO INDL=1,nDIM
                    KUU(ROW,COL)=KUU(ROW,COL) + 
     &      TangentStiffness(INDI,INDJ,INDK,INDL) * DSH(INDB,INDL) 
     &           * DSH(INDA,INDJ) * W(INTPT) * DetMAPJ
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
      ENDDO
        !
        ! End the loop over integration points
        !----------------------------------------------------------------
        ! RETURN to Abaqus the RHS vector and the StIFfness matrix.  This
        !  is essentially giving Abaqus the residual and the tangent matrix.
        !
        ! RETURN to Abaqus the right hand side vector.
        AMATRX=KUU
        DO INDI=1,nNODE
            RHS(nDIM * (INDI - 1) + 1,1)=RU(nDIM * (INDI - 1) + 1,1)
            RHS(nDIM * (INDI - 1) + 2,1)=RU(nDIM * (INDI - 1) + 2,1)
            RHS(nDIM * (INDI - 1) + 3,1)=RU(nDIM * (INDI - 1) + 3,1)
        ENDDO
      END SUBROUTINE UAX4      

!==================== U3D8 SUBROUTINE ====================     
      SUBROUTINE U3D8(RHS,AMATRX,SVARS,ENERGY,nDofEl,NRHS,NSVARS,
     &    PROPS,NPROPS,coords,MCRD,nNODE,UALL,DUALL,VEL,ACCN,
     &    JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,
     &    JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,
     &    MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD,
     &    nDIM,nINT,nINTS)   
      
      
      USE GLOBAL
      
      IMPLICIT NONE
      !
      ! variables defined in uel,passed back to Abaqus
      !
      REAL(KIND=8) :: RHS(MLVARX,*),AMATRX(nDofEl,nDofEl),SVARS(*),
     &   ENERGY(8),PNEWDT
      !
      ! variables passed into UEL
      !
      INTEGER :: nDofEl,NRHS,NSVARS,NPROPS,MCRD,nNODE,JTYPE,
     &   KSTEP,KINC
      INTEGER :: JELEM,NDLOAD,NPREDF,MLVARX,MDLOAD,NJPROP
      INTEGER :: JDLTYP(MDLOAD,*),LFLAGS(*),JPROPS(*)
      REAL(KIND=8) :: PROPS(*),coords(MCRD,nNODE),UALL(nDofEl),
     &   DUALL(MLVARX,*)
      REAL(KIND=8) :: VEL(nDofEl),ACCN(nDofEl),TIME(2),DTIME,
     &   PARAMS(*)
      REAL(KIND=8) :: ADLMAG(MDLOAD,*),PREDEF(2,NPREDF,nNODE),
     &   DDLMAG(MDLOAD,*),PERIOD
      !
      ! variables defined and used in the UEL
      !
      INTEGER :: INDI,INDJ,INDK,INDJJ,INDA,ROW,COL,INDB,INDL
      INTEGER :: INTPT,STAT
      REAL(KIND=8),PARAMETER :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0
      REAL(KIND=8),PARAMETER :: HALF=0.5D0,PI=3.141592653D0
      REAL(KIND=8),PARAMETER :: THREE=3.0D0,THIRD=1.0D0/3.0D0
      INTEGER :: nDIM 
      INTEGER :: nSDV,nlSDV,ngSDV 
      INTEGER :: nINT, nINTS 
      REAL(KIND=8) :: U(nNODE,nDIM),DU(nNODE,nDIM),
     &   Uold(nNODE,nDofEl),BODY(3)
      REAL(KIND=8) :: coordsC(MCRD,nNODE),RU(nDIM*nNODE,1),
     &   KUU(nDIM*nNODE,nDIM*nNODE),RUtemp(nDIM*nNODE,1),
     &   KUUtemp(nDIM*nNODE,nDIM*nNODE)
      REAL(KIND=8) :: KROND(3,3),XI(nINT,nDIM),W(nINT),
     &   SH0(nNODE),SH(nNODE)
      REAL(KIND=8) :: DSH0(nNODE,nDIM),DSHC0(nNODE,nDIM),
     &   DSH(nNODE,nDIM),DSHC(nNODE,nDIM),DSHXI(nNODE,nDIM)
      REAL(KIND=8) :: DetMAPJ0,DetMAPJ0C,DetMAPJ,DetMAPJC,
     &   Fc_Tau(3,3),DetFc_Tau,Fc_T(3,3),DetFc_T
      REAL(KIND=8) :: D_T(3),D0(3),D_Tau(3)
      REAL(KIND=8) :: F_Tau(3,3),F_T(3,3),DetF_Tau,
     &   T_Tau(3,3),FV_T(3,3),FV_Tau(3,3),DetF_T,DetF
      REAL(KIND=8) :: Q0,TangentStiffness(3,3,3,3),PKSTRESS(3,3)
      REAL(KIND=8) :: ZERO_MATRIX(1,3),Le
      REAL(KIND=8) :: SMAT(6,1),BMAT(6,nNode*nDim),Gmat(9,3*nNode),
     &   AMAT(9,9),QMAT(9,9)
      integer :: kk

      !
      ! Get element parameters
      !
      nlSDV = JPROPS(1) ! number of local sdv's per integ point
      ngSDV = JPROPS(2) ! number of global sdv's per integ point

      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
                !
                ! allocate memory for the globalSdv's
                !
                ! numElem needs to be set in the MODULE
                ! nInt needs to be set in the UEL
                !
                stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
          allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- U3D8 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif
      
      ! obtain initial conditions
      Q0=PROPS(9)
      D0(1)=COS(Q0)
      D0(2)=SIN(Q0)
      D0(3)=0   

      ! Initialize the residual and tangent matrices to zero
      RU=ZERO
      RUtemp = zero
      KUU=ZERO
      kuutemp = zero
      ENERGY=ZERO
      BODY=ZERO
      
      ! Identity tensor
      CALL CREATE_IDENTITY_MATRIX(3,KROND)

      ! Obtain nodal displacements
      INDK=0
      DO INDI=1,nNODE
       DO INDJ=1,nDIM
          INDK=INDK + 1
          U(INDI,INDJ)=UALL(INDK)
          DU(INDI,INDJ)=DUALL(INDK,1)
          Uold(INDI,INDJ)=U(INDI,INDJ) - DU(INDI,INDJ)
       ENDDO
      ENDDO
      
      ! Obtain current nodal coordinates
      DO INDI=1,nNODE
       DO INDJ=1,nDIM
          coordsC(INDJ,INDI)=coords(INDJ,INDI) + U(INDI,INDJ)
       ENDDO
      ENDDO
      
      ! displacement increment, based on element diagonal
      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) + 
     &     ((coordsC(2,1)-coordsC(2,7))**two) +
     &     ((coordsC(3,1)-coordsC(3,7))**two))
      
      DO INDI=1,nNODE
        DO INDJ=1,nDIM
            IF(dabs(DU(INDI,INDJ)).GT.10.0*Le) THEN
             PNEWDT=0.5
             RETURN
            ENDIF
        ENDDO
      ENDDO
      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto,E.A.,Peric,D.,Dutko,M.,Owen,D.R.J.,1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33,3277-3296.
      !
      !
      ZERO_MATRIX=0.0D0
      IF (nNODE.EQ.8) THEN
          CALL calcShape3DLinear(1,ZERO_MATRIX,1,SH0,DSHXI)
      ELSE
          print *, 'Incorrect number of nodes: nNODE.NE.8'
          call xit
      ENDIF
      

      ! Map shape functions from local to global reference coordinate system
      CALL MAPShape3D(nNODE,DSHXI,coords,DSH0,DetMAPJ0,STAT)
      IF (STAT.EQ.0) THEN
          PNEWDT=0.5
          RETURN
      ENDIF     
      
      ! Map shape functions from local to global current coordinate system
      CALL MAPShape3D(nNODE,DSHXI,coordsC,DSHC0,DetMAPJ0C,STAT)
      IF (STAT.EQ.0) THEN
          PNEWDT=0.5
          RETURN
      ENDIF
      

      ! Calculate the deformation gradient at the element centroid
      ! at the end of the increment for use in the `F-bar' method
      ! The subscript tau denotes the time at the end of the increment.

      Fc_Tau=KROND
      Fc_T=KROND
      DO INDI=1,nDIM
        DO INDJ=1,nDIM
            DO INDK=1,nNODE
              Fc_Tau(INDI,INDJ)=Fc_Tau(INDI,INDJ) + 
     &            DSH0(INDK,INDJ) * U(INDK,INDI)
              Fc_T(INDI,INDJ)=Fc_T(INDI,INDJ) +
     &            DSH0(INDK, INDJ)*Uold(INDK,INDI)
            ENDDO
        ENDDO
      ENDDO
      !
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      ! Begin the loop over integration points

      ! Obtain integration point local coordinates and weights
      IF(nNODE.EQ.8) then
         !
         ! gauss integration for a rectangular element
         !
         IF(nINT.EQ.8) then    
            CALL xint3D8pt(xi,w,nINT) ! 8-pt integration, nINT=8 above
         ELSEIF(nINT.EQ.1) then
            CALL xint3D1pt(xi,w,nINT) ! 1-pt integration, nINT=1 above
         ELSE
            print *, 'Invalid number of int points, nINT=',nINT
            call xit
         ENDIF
      ELSE
         print *, 'Incorrect number of nodes: nNODE.ne.8'
         call xit
      ENDIF

      
      ! Loop over integration points
      INDJJ=0  ! JJ is used for tracking the state variables
      DO INTPT=1,nINT
        ! Obtain state variables from previous increment
        IF ((KINC.LE.1) .AND. (KSTEP.EQ.1)) THEN
                ! This is the first increment of the first step.
                ! Give initial conditions.
                FV_T=KROND
                F_T=KROND
                D_T=D0
        ELSE
                ! This is not the first increment; read old values.
                FV_T(1,1)=SVARS(1 + INDJJ)
                FV_T(2,2)=SVARS(2 + INDJJ)
                FV_T(3,3)=SVARS(3 + INDJJ)
                FV_T(2,3)=SVARS(4 + INDJJ)
                FV_T(3,1)=SVARS(5 + INDJJ)
                FV_T(1,2)=SVARS(6 + INDJJ)
                FV_T(3,2)=SVARS(7 + INDJJ)
                FV_T(1,3)=SVARS(8 + INDJJ)
                FV_T(2,1)=SVARS(9 + INDJJ)
                F_T(1,1)=SVARS(10 + INDJJ)
                F_T(2,2)=SVARS(11 + INDJJ)
                F_T(3,3)=SVARS(12 + INDJJ)
                F_T(2,3)=SVARS(13 + INDJJ)
                F_T(3,1)=SVARS(14 + INDJJ)
                F_T(1,2)=SVARS(15 + INDJJ)
                F_T(3,2)=SVARS(16 + INDJJ)
                F_T(1,3)=SVARS(17 + INDJJ)
                F_T(2,1)=SVARS(18 + INDJJ)
                D_T(1)=SVARS(19 + INDJJ)
                D_T(2)=SVARS(20 + INDJJ)
                D_T(3)=SVARS(21 + INDJJ)
        END IF
        
        
        ! Obtain shape functions and their local gradients
        IF (nNODE.EQ.8) THEN
            CALL calcShape3DLinear(nINT,XI,INTPT,SH,DSHXI)
        ELSE
            print *, 'Incorrect number of nodes: nNODE /= 8'
            call xit
        ENDIF
         

        ! Map shape functions from local to global reference coordinate system
        CALL MAPShape3D(nNODE,DSHXI,coords,DSH,DetMAPJ,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF
        
        
        ! Map shape functions from local to global current coordinate system
        CALL MAPShape3D(nNODE,DSHXI,coordsC,DSHC,DetMAPJC,STAT)
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF     

        ! Obtain the deformation gradient at this integration point.
        ! The subscript tau denotes the time at the end of the increment.

        F_Tau=KROND
        F_t=KROND
        DO INDI=1,nDIM
            DO INDJ=1,nDIM
                DO INDK=1,nNODE
                    F_Tau(INDI,INDJ)=F_Tau(INDI,INDJ) +
     &               DSH(INDK,INDJ) * U(INDK,INDI)
                    F_T(INDI,INDJ)=F_T(INDI,INDJ) +
     &               DSH(INDK,INDJ) * Uold(INDK,INDI)
                ENDDO
            ENDDO
        ENDDO       
        !
        ! Modify the deformation gradient for the `F-bar' method
        !  only when using the 8 node fully integrated linear
        !  element, do not use the `F-bar' method for any other element
        !
        if((nNode.EQ.8).and.(nInt.EQ.8)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
        !    F_tau = ((detFc_tau/detF_tau)**third)*F_tau
        !    F_t = ((detFc_tau/detF_tau)**third)*F_t
            F_tau = F_tau
            F_t = F_t
        endif
        call mDet(F_Tau,DetF_Tau)
        !
        ! Perfom the constitutive time integration at integration point
        !
        CALL MATERIAL(PROPS,F_T,FV_T,D_T,F_Tau,FV_Tau,D_Tau,
     &    DTIME,D0,PKSTRESS,TangentStiffness,STAT) 
        IF (STAT.EQ.0) THEN
            PNEWDT=0.5
            RETURN
        ENDIF
        
        ! Save the state variables at this integ point
        !  at the end of the increment.
        !
        SVARS(1 + INDJJ)=FV_Tau(1,1)
        SVARS(2 + INDJJ)=FV_Tau(2,2)
        SVARS(3 + INDJJ)=FV_Tau(3,3)
        SVARS(4 + INDJJ)=FV_Tau(2,3)
        SVARS(5 + INDJJ)=FV_Tau(3,1)
        SVARS(6 + INDJJ)=FV_Tau(1,2)
        SVARS(7 + INDJJ)=FV_Tau(3,2)
        SVARS(8 + INDJJ)=FV_Tau(1,3)
        SVARS(9 + INDJJ)=FV_Tau(2,1)
        SVARS(10 + INDJJ)=F_Tau(1,1)
        SVARS(11 + INDJJ)=F_Tau(2,2)
        SVARS(12 + INDJJ)=F_Tau(3,3)
        SVARS(13 + INDJJ)=F_Tau(2,3)
        SVARS(14 + INDJJ)=F_Tau(3,1)
        SVARS(15 + INDJJ)=F_Tau(1,2)
        SVARS(16 + INDJJ)=F_Tau(3,2)
        SVARS(17 + INDJJ)=F_Tau(1,3)
        SVARS(18 + INDJJ)=F_Tau(2,1)
        SVARS(19 + INDJJ)=D_Tau(1)
        SVARS(20 + INDJJ)=D_Tau(2)
        SVARS(21 + INDJJ)=D_Tau(3)
        INDJJ=INDJJ + nlSDV
        
        ! Save the state variables at this integ point in the
        !  global array used for plotting field output
        !
        globalSdv(jelem,intPt,1) = D_Tau(1)

        ! Compute/update the displacement residual vector
        DO INDA=1,nNODE
            DO INDI=1,nDIM
            ROW=nDIM * (INDA - 1) + INDI
              DO INDJ=1,nDIM
            RU(ROW,1)=RU(ROW,1) - PKSTRESS(INDI,INDJ) *
     &         DSH(INDA,INDJ) * W(INTPT) * DetMAPJ
              ENDDO
            ENDDO
        ENDDO

         Smat(1,1) = PKstress(1,1)
         Smat(2,1) = PKstress(2,2)
         Smat(3,1) = PKstress(3,3)
         Smat(4,1) = PKstress(1,2)
         Smat(5,1) = PKstress(2,3)
         Smat(6,1) = PKstress(1,3)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dsh(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dsh(kk,2)
            Bmat(3,3+nDim*(kk-1)) = dsh(kk,3)
            Bmat(4,1+nDim*(kk-1)) = dsh(kk,2)
            Bmat(4,2+nDim*(kk-1)) = dsh(kk,1)
            Bmat(5,2+nDim*(kk-1)) = dsh(kk,3)
            Bmat(5,3+nDim*(kk-1)) = dsh(kk,2)
            Bmat(6,1+nDim*(kk-1)) = dsh(kk,3)
            Bmat(6,3+nDim*(kk-1)) = dsh(kk,1)
         enddo
       Rutemp = Rutemp + detmapJ*w(intpt)*
     &        (
     &        -matmul(transpose(Bmat),Smat)
     &        )
         
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dsh(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dsh(kk,1)
            Gmat(3,3+nDim*(kk-1)) = dsh(kk,1)
            Gmat(4,1+nDim*(kk-1)) = dsh(kk,2)
            Gmat(5,2+nDim*(kk-1)) = dsh(kk,2)
            Gmat(6,3+nDim*(kk-1)) = dsh(kk,2)
            Gmat(7,1+nDim*(kk-1)) = dsh(kk,3)
            Gmat(8,2+nDim*(kk-1)) = dsh(kk,3)
            Gmat(9,3+nDim*(kk-1)) = dsh(kk,3)
         enddo

         Amat = zero
         Amat(1,1) = TangentStiffness(1,1,1,1)
         Amat(1,2) = TangentStiffness(1,1,2,1)
         Amat(1,3) = TangentStiffness(1,1,3,1)
         Amat(1,4) = TangentStiffness(1,1,1,2)
         Amat(1,5) = TangentStiffness(1,1,2,2)
         Amat(1,6) = TangentStiffness(1,1,3,2)
         Amat(1,7) = TangentStiffness(1,1,1,3)
         Amat(1,8) = TangentStiffness(1,1,2,3)
         Amat(1,9) = TangentStiffness(1,1,3,3)
         Amat(2,1) = TangentStiffness(2,1,1,1)
         Amat(2,2) = TangentStiffness(2,1,2,1)
         Amat(2,3) = TangentStiffness(2,1,3,1)
         Amat(2,4) = TangentStiffness(2,1,1,2)
         Amat(2,5) = TangentStiffness(2,1,2,2)
         Amat(2,6) = TangentStiffness(2,1,3,2)
         Amat(2,7) = TangentStiffness(2,1,1,3)
         Amat(2,8) = TangentStiffness(2,1,2,3)
         Amat(2,9) = TangentStiffness(2,1,3,3)
         Amat(3,1) = TangentStiffness(3,1,1,1)
         Amat(3,2) = TangentStiffness(3,1,2,1)
         Amat(3,3) = TangentStiffness(3,1,3,1)
         Amat(3,4) = TangentStiffness(3,1,1,2)
         Amat(3,5) = TangentStiffness(3,1,2,2)
         Amat(3,6) = TangentStiffness(3,1,3,2)
         Amat(3,7) = TangentStiffness(3,1,1,3)
         Amat(3,8) = TangentStiffness(3,1,2,3)
         Amat(3,9) = TangentStiffness(3,1,3,3)
         Amat(4,1) = TangentStiffness(1,2,1,1)
         Amat(4,2) = TangentStiffness(1,2,2,1)
         Amat(4,3) = TangentStiffness(1,2,3,1)
         Amat(4,4) = TangentStiffness(1,2,1,2)
         Amat(4,5) = TangentStiffness(1,2,2,2)
         Amat(4,6) = TangentStiffness(1,2,3,2)
         Amat(4,7) = TangentStiffness(1,2,1,3)
         Amat(4,8) = TangentStiffness(1,2,2,3)
         Amat(4,9) = TangentStiffness(1,2,3,3)
         Amat(5,1) = TangentStiffness(2,2,1,1)
         Amat(5,2) = TangentStiffness(2,2,2,1)
         Amat(5,3) = TangentStiffness(2,2,3,1)
         Amat(5,4) = TangentStiffness(2,2,1,2)
         Amat(5,5) = TangentStiffness(2,2,2,2)
         Amat(5,6) = TangentStiffness(2,2,3,2)
         Amat(5,7) = TangentStiffness(2,2,1,3)
         Amat(5,8) = TangentStiffness(2,2,2,3)
         Amat(5,9) = TangentStiffness(2,2,3,3)
         Amat(6,1) = TangentStiffness(3,2,1,1)
         Amat(6,2) = TangentStiffness(3,2,2,1)
         Amat(6,3) = TangentStiffness(3,2,3,1)
         Amat(6,4) = TangentStiffness(3,2,1,2)
         Amat(6,5) = TangentStiffness(3,2,2,2)
         Amat(6,6) = TangentStiffness(3,2,3,2)
         Amat(6,7) = TangentStiffness(3,2,1,3)
         Amat(6,8) = TangentStiffness(3,2,2,3)
         Amat(6,9) = TangentStiffness(3,2,3,3)
         Amat(7,1) = TangentStiffness(1,3,1,1)
         Amat(7,2) = TangentStiffness(1,3,2,1)
         Amat(7,3) = TangentStiffness(1,3,3,1)
         Amat(7,4) = TangentStiffness(1,3,1,2)
         Amat(7,5) = TangentStiffness(1,3,2,2)
         Amat(7,6) = TangentStiffness(1,3,3,2)
         Amat(7,7) = TangentStiffness(1,3,1,3)
         Amat(7,8) = TangentStiffness(1,3,2,3)
         Amat(7,9) = TangentStiffness(1,3,3,3)
         Amat(8,1) = TangentStiffness(2,3,1,1)
         Amat(8,2) = TangentStiffness(2,3,2,1)
         Amat(8,3) = TangentStiffness(2,3,3,1)
         Amat(8,4) = TangentStiffness(2,3,1,2)
         Amat(8,5) = TangentStiffness(2,3,2,2)
         Amat(8,6) = TangentStiffness(2,3,3,2)
         Amat(8,7) = TangentStiffness(2,3,1,3)
         Amat(8,8) = TangentStiffness(2,3,2,3)
         Amat(8,9) = TangentStiffness(2,3,3,3)
         Amat(9,1) = TangentStiffness(3,3,1,1)
         Amat(9,2) = TangentStiffness(3,3,2,1)
         Amat(9,3) = TangentStiffness(3,3,3,1)
         Amat(9,4) = TangentStiffness(3,3,1,2)
         Amat(9,5) = TangentStiffness(3,3,2,2)
         Amat(9,6) = TangentStiffness(3,3,3,2)
         Amat(9,7) = TangentStiffness(3,3,1,3)
         Amat(9,8) = TangentStiffness(3,3,2,3)
         Amat(9,9) = TangentStiffness(3,3,3,3)
         
          Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         
         Kuutemp = Kuutemp + detMapJ*w(intpt)*
     &           (
     &           matmul(matmul(transpose(Gmat),Amat),Gmat)
     &           )
        ! Compute/update the displacement tangent matrix
        DO INDA=1,nNODE
          DO INDI=1,nDIM
            DO INDB=1,nNODE
              DO INDK=1,nDIM
                ROW=nDIM * (INDA - 1) + INDI
                COL=nDIM * (INDB - 1) + INDK
                  DO INDJ=1,nDIM
                    DO INDL=1,nDIM
                    KUU(ROW,COL)=KUU(ROW,COL) + 
     &      TangentStiffness(INDI,INDJ,INDK,INDL) * DSH(INDB,INDL) 
     &           * DSH(INDA,INDJ) * W(INTPT) * DetMAPJ
                    ENDDO
                  ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      
      
      ENDDO

        !
        ! End the loop over integration points
        !----------------------------------------------------------------
        ! RETURN to Abaqus the RHS vector and the StIFfness matrix.  This
        !  is essentially giving Abaqus the residual and the tangent matrix.
        !
        ! RETURN to Abaqus the right hand side vector.
        AMATRX=KUU
        DO INDI=1,nNODE
        RHS(nDIM * (INDI - 1) + 1,1)=RU(nDIM * (INDI - 1) + 1,1)
        RHS(nDIM * (INDI - 1) + 2,1)=RU(nDIM * (INDI - 1) + 2,1)
        RHS(nDIM * (INDI - 1) + 3,1)=RU(nDIM * (INDI - 1) + 3,1)
        ENDDO
        
    
      END SUBROUTINE U3D8             
                  
!==================== Material Model ====================        
      SUBROUTINE MATERIAL(PROPS,FN,FVN,DN,FITER,FVITER,DITER,
     &  DTIME,D0,PKSTRESS,TangentStiffness,STAT)  
        IMPLICIT NONE
        
        INTEGER :: STAT,INFO,N,NRHS,FLAG
        REAL(KIND=8),PARAMETER :: PI=3.1415927D0,TOL=1D-10
        REAL(KIND=8) :: Det,DTIME
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,FVITER,FITER,
     &   DELTAFV,CAUSTRESS,INVFITER,PKSTRESS
        REAL(KIND=8),DIMENSION(3) :: D0,DN,DITER,DELTAD,D
        REAL(KIND=8),DIMENSION(12) :: QITER,RR,DELTAQ
        REAL(KIND=8),DIMENSION(9) :: TEMP91,RG,TDELTAFVECTOR
        REAL(KIND=8),DIMENSION(8) :: RES8,DELTAFVECTOR
        REAL(KIND=8),DIMENSION(3,3) :: RES2,STIFFFD,FV,F
        REAL(KIND=8),DIMENSION(3,3,3) :: RES3
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES4,TangentStiffness, NumK
        REAL(KIND=8),DIMENSION(3) :: RES1,RF
        REAL(KIND=8),DIMENSION(3,9) :: RES39,STIFFFFV,
     &        STIFFFF,STIFFDF
        REAL(KIND=8),DIMENSION(9,3) :: RES93,STIFFGD,
     &        STIFFSIGMAD,STIFFPD
        REAL(KIND=8),DIMENSION(9,9) :: RES99,STIFFGFV,STIFFGF,
     &        STIFFSIGMAF,STIFFSIGMAFV,STIFFPF,
     &        STIFFPFV,STIFFFVF,TSTIFFPF
        REAL(KIND=8),DIMENSION(12,12) :: KQ,INVKQ
        REAL(KIND=8),DIMENSION(12,9) :: KF,M
    
        CALL INTERNALVARIABLES(PROPS,FN,DN,FVN,FITER,DTIME,D0,
     &   FVITER,DITER) 
        
        CALL SIGMA(PROPS,FITER,FVITER,DITER,D0,CAUSTRESS)

        CALL MATRIX_INVERSE(FITER,INVFITER,3,INFO)
        CALL MDet(FITER,Det)
        PKSTRESS= Det*MATMUL(CAUSTRESS,TRANSPOSE(INVFITER))
 
        CALL KFF(PROPS,FN,FVN,DN,FITER,FVITER,DITER,DTIME,
     &        D0,RES3)
        CALL TENSOR2MATRIX12(RES3,RES39)
        STIFFFF=RES39        
        CALL KFD(PROPS,FN,FVN,DN,FITER,FVITER,DITER,DTIME,
     &        D0,RES2)
        STIFFFD=RES2
        CALL KFV(PROPS,FN,FVN,DN,FITER,FVITER,DITER,DTIME,
     &        D0,RES3) 
        CALL TENSOR2MATRIX12(RES3,RES39)
        STIFFFFV=RES39

        CALL KGF(PROPS,FITER,DITER,FVITER,DTIME,D0,RES4)
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFGF=RES99
        
        CALL KGV(PROPS,FVN,FITER,DITER,FVITER,DTIME,D0,RES4)
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFGFV=RES99
                
        CALL KGD(PROPS,FITER,DITER,FVITER,DTIME,D0,RES3)
        CALL TENSOR2MATRIX21(RES3,RES93)
        STIFFGD=RES93
        CALL PARSIGMAPARF(PROPS,FITER,FVITER,DITER,D0,RES4)
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFSIGMAF=RES99

        CALL PARSIGMAPARD(PROPS,FITER,FVITER,DITER,D0,RES3)
        CALL TENSOR2MATRIX21(RES3,RES93)
        STIFFSIGMAD=RES93

        CALL PARSIGMAPARFV(PROPS,FITER,FVITER,DITER,D0,RES4)
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFSIGMAFV=RES99

        CALL PARPPARD(PROPS,FITER,FVITER,DITER,D0,RES3)
        CALL TENSOR2MATRIX21(RES3,RES93)
        STIFFPD=RES93
        
        CALL PARPPARF(PROPS,FITER,FVITER,DITER,D0,RES4) 
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFPF=RES99

        CALL PARPPARFV(PROPS,FITER,FVITER,DITER,D0,RES4)
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFPFV=RES99

        KQ=0.0D0
        KQ(1:3,1:3)=STIFFFD
        KQ(1:3,4:12)=STIFFFFV
        KQ(4:12,1:3)=STIFFGD
        KQ(4:12,4:12)=STIFFGFV
        KF=0.0D0
        KF(1:3,1:9)=STIFFFF
        KF(4:12,1:9)=STIFFGF
        N=12
        NRHS=9

        CALL SOLVE_LINEAR_SYSTEM(KQ,KF,M,N,NRHS,INFO)
        M=-M
       
        CALL MATRIX2VECTOR(DELTAFV,TEMP91)
        DELTAQ(1:3)=DELTAD
        DELTAQ(4:12)=TEMP91
        STIFFDF=M(1:3,:)
        STIFFFVF=M(4:12,:)
        
        TSTIFFPF=STIFFPF + MATMUL(STIFFPD,STIFFDF) 
     &    + MATMUL(STIFFPFV,STIFFFVF)
        CALL MATRIX2TENSOR(TSTIFFPF,TangentStiffness)
        CALL MDet(FVITER,Det)
        IF (Det <= 0.0D0) THEN
            WRITE(*,*) 'Det(Fviter) <= zero in INTEG'
            STAT=0
            RETURN
        END IF
      END SUBROUTINE MATERIAL    
      
!==================== Numerical tangent ====================  
      SUBROUTINE NKPF(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
      REAL(KIND=8),DIMENSION(12) :: PROPS
      REAL(KIND=8),DIMENSION(3,3) :: F,FV
      REAL(KIND=8),DIMENSION(3) :: D,D0
      COMPLEX,DIMENSION(3,3) :: pertF, caustress, invpertF,pertP
      REAL(KIND=8),DIMENSION(3,3,3,3) :: RES
      REAL :: h
      COMPLEX :: eps, Det_pertF
      INTEGER :: INDI,INDJ,INDK,INDL,ISTAT
      
      h = 1E-10
      eps = CMPLX(0.0, h, kind=kind(h))
      res = 0.0
      
      DO INDI = 1,3
        DO INDJ = 1,3
          DO INDK = 1,3
            DO INDL = 1,3
                pertF = F
                pertF(INDK, INDL) = pertF(INDK, INDL) + eps
                CALL ComplexSIGMA(PROPS,pertF,FV,D,D0,caustress)
                call ComplexMATINV3D(pertF,invpertF,Det_pertF,ISTAT)
                pertP = Det_pertF*MATMUL(caustress, invpertF)
                RES(INDI,INDJ,INDK,INDL) = aimag(pertP(INDI,INDJ))/h
            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE NKPF
      
!==================== Internal Variables ====================     
      SUBROUTINE INTERNALVARIABLES(PROPS,FN,DN,FVN,FITER,DTIME,
     & D0,FVITER,DITER)
        IMPLICIT NONE
        INTEGER :: INNERITER,LCNT
        INTEGER :: INDI,INDJ,INFO
        INTEGER,PARAMETER :: MAXITER=200,LCNTMAX=200
        REAL(KIND=8),PARAMETER :: PI=3.1415927D0,RTOL=1D-10,ATOL=1D-10
        REAL(KIND=8) :: INNERRES,Det,DTIME,REST
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: INNERERROR,INIRES
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,FVITER,FITER,
     &    DELTAFV,CAUSTRESS,INVFITER,PKSTRESS,FVITERtemp
        REAL(KIND=8),DIMENSION(3) :: D0,DN,DITER,DELTAD,D,DITERtemp
        REAL(KIND=8),DIMENSION(12) :: QITER,RR,DELTAQ,QITERtemp,DQITER
        REAL(KIND=8),DIMENSION(9) :: TEMP91,RG,TDELTAFVECTOR
        REAL(KIND=8),DIMENSION(8) :: RES8,DELTAFVECTOR
        REAL(KIND=8),DIMENSION(3,3) :: RES2,STIFFFD,FV
        REAL(KIND=8),DIMENSION(3,3,3) :: RES3
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES4
        REAL(KIND=8),DIMENSION(3) :: RES1,RF
        REAL(KIND=8),DIMENSION(3,9) :: RES39,STIFFFFV,
     &   STIFFFF,STIFFDF
        REAL(KIND=8),DIMENSION(9,3) :: RES93,STIFFGD,STIFFSIGMAD,
     &   STIFFPD
        REAL(KIND=8),DIMENSION(9,9) :: RES99,STIFFGFV,STIFFGF,
     &   STIFFSIGMAF,STIFFSIGMAFV,STIFFPF,STIFFPFV,
     &   STIFFFVF,TSTIFFPF
        REAL(KIND=8),DIMENSION(12,12) :: KQ,INVKQ
        REAL(KIND=8),DIMENSION(12,9) :: KF,M
        INTEGER :: io_status

        FVITER=FVN
        DITER=DN

        INNERRES=1.0D0
        INNERITER=1
        QITER(1:3)=DITER
        CALL MATRIX2VECTOR(FVITER,TEMP91)
        QITER(4:12)=TEMP91
        ! Compute residual at 0 iteration
        CALL FUNCF(PROPS,FN,DN,FVN,FITER,FVITER,DITER,DTIME,
     &  D0,RF)
        CALL FUNCG(PROPS,FVN,FITER,FVITER,DITER,DTIME,D0,RES2)
        CALL MATRIX2VECTOR(RES2,RG)
        RR(1:3)=RF
        RR(4:12)=RG
        CALL EUNORM(RR,INNERRES,12)         
        INIRES=ABS(INNERRES)
        INNERRES=INIRES
        REST=1.0D0
      
        DO WHILE (INNERITER.LT.MAXITER)
        CALL FUNCF(PROPS,FN,DN,FVN,FITER,FVITER,DITER,DTIME,
     &  D0,RF)
        CALL FUNCG(PROPS,FVN,FITER,FVITER,DITER,DTIME,D0,RES2)
        CALL MATRIX2VECTOR(RES2,RG)
        RR(1:3)=RF
        RR(4:12)=RG
        CALL KFD(PROPS,FN,FVN,DN,FITER,FVITER,DITER,DTIME,
     &    D0,RES2)
        STIFFFD=RES2
        CALL KFV(PROPS,FN,FVN,DN,FITER,FVITER,DITER,DTIME,
     &   D0,RES3) 
        CALL TENSOR2MATRIX12(RES3,RES39)
        STIFFFFV=RES39
        CALL KGV(PROPS,FVN,FITER,DITER,FVITER,DTIME,D0,RES4)
        CALL TENSOR2MATRIX22(RES4,RES99)
        STIFFGFV=RES99
        CALL KGD(PROPS,FITER,DITER,FVITER,DTIME,D0,RES3)
        CALL TENSOR2MATRIX21(RES3,RES93)
        STIFFGD=RES93
        KQ=0.0D0
        KQ(1:3,1:3)=STIFFFD
        KQ(1:3,4:12)=STIFFFFV
        KQ(4:12,1:3)=STIFFGD
        KQ(4:12,4:12)=STIFFGFV
        CALL MATRIX_INVERSE(KQ,INVKQ,12,INFO)
        DQITER=-MATMUL(INVKQ,RR)
        QITERtemp=QITER+DQITER
        DITERtemp=QITERtemp(1:3)
        CALL VECTOR2MATRIX(QITERtemp(4:12),FVITERtemp)                
        CALL FUNCF(PROPS,FN,DN,FVN,FITER,FVITERtemp,DITERtemp,DTIME,
     &  D0,RF)
        CALL FUNCG(PROPS,FVN,FITER,FVITERtemp,DITERtemp,DTIME,D0,RES2)
        CALL MATRIX2VECTOR(RES2,RG)
        RR(1:3)=RF
        RR(4:12)=RG                
        CALL EUNORM(RR,INNERRES,12)         
        INNERRES=ABS(INNERRES)   
        LCNT = 1
        DO WHILE ((INNERRES.GT.1.0*REST).AND.(INNERRES.GT.ATOL).AND.
     &    (LCNT.LT.LCNTMAX))
            REST=INNERRES
            DQITER = 1.0D0/3.0D0*DQITER
            QITERtemp=QITER+DQITER
            DITERtemp=QITERtemp(1:3)
            CALL VECTOR2MATRIX(QITERtemp(4:12),FVITERtemp)            
            CALL FUNCF(PROPS,FN,DN,FVN,FITER,FVITERtemp,DITERtemp,DTIME,
     &       D0,RF)
            CALL FUNCG(PROPS,FVN,FITER,FVITERtemp,DITERtemp,DTIME,
     &       D0,RES2)
            CALL MATRIX2VECTOR(RES2,RG)
            RR(1:3)=RF
            RR(4:12)=RG                
            CALL EUNORM(RR,INNERRES,12)         
            INNERRES=ABS(INNERRES)   
            QITER=QITERtemp
        ENDDO   
        INNERITER=INNERITER +1
        REST=INNERRES
        QITER=QITERtemp
        FVITER=FVITERtemp
        DITER=DITERtemp
        IF ((INNERRES.LT.(RTOL*INIRES)).OR.(INNERRES.LT.aTOL)) THEN
                EXIT
        ENDIF
        IF (INNERITER.EQ.MAXITER) THEN
            print *, innerres
            PRINT *, 'diverge locally'
            CALL XIT
        END IF    
        ENDDO   
      END SUBROUTINE INTERNALVARIABLES     

           
!==================== 3D INVERSE AND Det ====================    
      SUBROUTINE MATINV3D(A,A_INV,Det_A,ISTAT)
        ! RETURNs A_INV,the inverse and Det_A,the Determinant
        ! Note that the Determinant is of the original matrix,not the inverse

        IMPLICIT NONE

        INTEGER :: ISTAT
        REAL(KIND=8) :: A(3,3),A_INV(3,3),Det_A,Det_A_INV

        ISTAT=1

        Det_A=A(1,1) * (A(2,2) * A(3,3) - A(3,2) * A(2,3)) - 
     &        A(2,1) * (A(1,2) * A(3,3) - A(3,2) * A(1,3)) + 
     &        A(3,1) * (A(1,2) * A(2,3) - A(2,2) * A(1,3))

        IF (Det_A.LE.0.0D0) THEN
            print *, 'WARNING: SUBROUTINE MATINV3D:'
            print *, 'WARNING: Det of MAT=',Det_A
            ISTAT=0
            RETURN
        END IF

        Det_A_INV=1.0D0 / Det_A

        A_INV(1,1)=Det_A_INV*(A(2,2) * A(3,3) - A(3,2) * A(2,3))
        A_INV(1,2)=Det_A_INV*(A(3,2) * A(1,3) - A(1,2) * A(3,3))
        A_INV(1,3)=Det_A_INV*(A(1,2) * A(2,3) - A(2,2) * A(1,3))
        A_INV(2,1)=Det_A_INV*(A(3,1) * A(2,3) - A(2,1) * A(3,3))
        A_INV(2,2)=Det_A_INV*(A(1,1) * A(3,3) - A(3,1) * A(1,3))
        A_INV(2,3)=Det_A_INV*(A(2,1) * A(1,3) - A(1,1) * A(2,3))
        A_INV(3,1)=Det_A_INV*(A(2,1) * A(3,2) - A(3,1) * A(2,2))
        A_INV(3,2)=Det_A_INV*(A(3,1) * A(1,2) - A(1,1) * A(3,2))
        A_INV(3,3)=Det_A_INV*(A(1,1) * A(2,2) - A(2,1) * A(1,2))
      END SUBROUTINE MATINV3D 

!==================== 3D INVERSE AND Det ====================    
      SUBROUTINE ComplexMATINV3D(A,A_INV,Det_A,ISTAT)
        ! RETURNs A_INV,the inverse and Det_A,the Determinant
        ! Note that the Determinant is of the original matrix,not the inverse

        IMPLICIT NONE

        INTEGER :: ISTAT
        COMPLEX :: A(3,3),A_INV(3,3),Det_A,Det_A_INV

        ISTAT=1

        Det_A=A(1,1) * (A(2,2) * A(3,3) - A(3,2) * A(2,3)) - 
     &        A(2,1) * (A(1,2) * A(3,3) - A(3,2) * A(1,3)) + 
     &        A(3,1) * (A(1,2) * A(2,3) - A(2,2) * A(1,3))


        Det_A_INV=1.0D0 / Det_A

        A_INV(1,1)=Det_A_INV*(A(2,2) * A(3,3) - A(3,2) * A(2,3))
        A_INV(1,2)=Det_A_INV*(A(3,2) * A(1,3) - A(1,2) * A(3,3))
        A_INV(1,3)=Det_A_INV*(A(1,2) * A(2,3) - A(2,2) * A(1,3))
        A_INV(2,1)=Det_A_INV*(A(3,1) * A(2,3) - A(2,1) * A(3,3))
        A_INV(2,2)=Det_A_INV*(A(1,1) * A(3,3) - A(3,1) * A(1,3))
        A_INV(2,3)=Det_A_INV*(A(2,1) * A(1,3) - A(1,1) * A(2,3))
        A_INV(3,1)=Det_A_INV*(A(2,1) * A(3,2) - A(3,1) * A(2,2))
        A_INV(3,2)=Det_A_INV*(A(3,1) * A(1,2) - A(1,1) * A(3,2))
        A_INV(3,3)=Det_A_INV*(A(1,1) * A(2,2) - A(2,1) * A(1,2))
      END SUBROUTINE ComplexMATINV3D         
                
!==================== Euclidean norm ====================    
      SUBROUTINE EUNORM(VECTOR,NORM,N)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: N
        REAL(KIND=8),DIMENSION(N),INTENT(IN) :: VECTOR
        REAL(KIND=8),INTENT(OUT) :: NORM
        NORM=SQRT(SUM(VECTOR**2))
      END SUBROUTINE EUNORM           

!==================== MATRIX TO TENSOR ====================
      SUBROUTINE MATRIX2TENSOR(F_MATRIX,R)
        IMPLICIT NONE
        REAL(KIND=8),INTENT(IN) :: F_MATRIX(9,9)  ! Assuming F_MATRIX is a 9x9 matrix
        REAL(KIND=8),INTENT(OUT) :: R(3,3,3,3)  ! 4D Tensor output
        INTEGER :: INDEX(9,2)
        INTEGER :: I

        ! Define the index matrix
        INDEX=RESHAPE((/  
     &  1,2,3,2,3,1,3,1,2,
     &   1,2,3,3,1,2,2,3,1  
     &   /),(/9,2/))

        ! Initialize R to zeros
        R=0.0D0

        ! Loop to fill the tensor R
        DO I=1,9
            R(INDEX(I,1),INDEX(I,2),1,1)=F_MATRIX(I,1)
            R(INDEX(I,1),INDEX(I,2),2,2)=F_MATRIX(I,2)
            R(INDEX(I,1),INDEX(I,2),3,3)=F_MATRIX(I,3)
            R(INDEX(I,1),INDEX(I,2),2,3)=F_MATRIX(I,4)
            R(INDEX(I,1),INDEX(I,2),3,1)=F_MATRIX(I,5)
            R(INDEX(I,1),INDEX(I,2),1,2)=F_MATRIX(I,6)
            R(INDEX(I,1),INDEX(I,2),3,2)=F_MATRIX(I,7)
            R(INDEX(I,1),INDEX(I,2),1,3)=F_MATRIX(I,8)
            R(INDEX(I,1),INDEX(I,2),2,1)=F_MATRIX(I,9)
        ENDDO
      END SUBROUTINE MATRIX2TENSOR    

!==================== TENSOR TO MATRIX ==================== 
      SUBROUTINE TENSOR2MATRIX22(F,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3,3,3) :: F
        REAL(KIND=8),DIMENSION(9,9) :: RES
        INTEGER :: I,J
        INTEGER,DIMENSION(9,2) :: INDEX
        INDEX=RESHAPE((/ 
     &       1,2,3,2,3,1,3,1,2,
     &       1,2,3,3,1,2,2,3,1  
     &       /),(/9,2/))
        RES=0.0D0
        DO I=1,9
        DO J=1,9
            RES(I,J)=F(INDEX(I,1),INDEX(I,2),INDEX(J,1),INDEX(J,2))
        ENDDO
        ENDDO
      END SUBROUTINE TENSOR2MATRIX22 
        
      SUBROUTINE TENSOR2MATRIX12(F,RES)
        IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3,3) :: F
        REAL(KIND=8),DIMENSION(3,9) :: RES
        INTEGER :: I,J
        INTEGER,DIMENSION(9,2) :: INDEX
        INDEX=RESHAPE((/ 
     &       1,2,3,2,3,1,3,1,2,
     &       1,2,3,3,1,2,2,3,1  
     &       /),(/9,2/))
        RES=0.0D0
        DO I=1,3
        DO J=1,9
            RES(I,J)=F(I,INDEX(J,1),INDEX(J,2))
        ENDDO
        ENDDO
      END SUBROUTINE TENSOR2MATRIX12
    
      SUBROUTINE TENSOR2MATRIX21(F,RES)
        IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3,3) :: F
        REAL(KIND=8),DIMENSION(9,3) :: RES
        INTEGER :: I,J
        INTEGER,DIMENSION(9,2) :: INDEX
        INDEX=RESHAPE((/ 
     &        1,2,3,2,3,1,3,1,2,
     &       1,2,3,3,1,2,2,3,1  
     &       /),(/9,2/))
        RES=0.0D0
        DO I=1,9
        DO J=1,3
            RES(I,J)=F(INDEX(I,1),INDEX(I,2),J)
        ENDDO
        ENDDO
      END SUBROUTINE TENSOR2MATRIX21    

!==================== KGF ====================     
      SUBROUTINE KGF(PROPS,F,D,FV,DELTAT,D0,RES)  
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUNEQ,ETAN,LPA,LPER,IN0,INE,
     &    GCONST,DELTAT
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L,L0,
     &   GINE,INEF,INVL
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,K21,FEF,
     &   TEMPK22,K22
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(9,9) :: RES99
        INTEGER :: INFO
        INTEGER :: INDALPHA,INDBETA,INDKAPPA,INDA,INDI,
     &   INDJ,INDGAMMA
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        ETAN=PROPS(6)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL PARGPARINE(INE,L0,FE,INVL,PROPS,GINE)
        CALL PARINEPARF(L0,INVL,FE,FV,INEF)
        K21=0.0D0
        DO INDALPHA=1,3
          DO INDBETA=1,3
            DO INDKAPPA=1,3
             DO INDA=1,3
            K21(INDALPHA,INDBETA,INDKAPPA,INDA)=
     &       K21(INDALPHA,INDBETA,INDKAPPA,INDA) + 
     &        GINE(INDALPHA,INDBETA) * INEF(INDKAPPA,INDA)
             ENDDO
            ENDDO
          ENDDO
         ENDDO
        CALL PARFEPARF(FV,FEF)
        TEMPK22=0.0D0
        DO INDI=1,3
        DO INDJ=1,3
        DO INDGAMMA=1,3
        DO INDBETA=1,3
        DO INDKAPPA=1,3
        DO INDA=1,3
        DO INDALPHA=1,3
            TEMPK22(INDALPHA,INDBETA,INDKAPPA,INDA)=
     &       TEMPK22(INDALPHA,INDBETA,INDKAPPA,INDA) + 
     &           INVL(INDI,INDJ) * L0(INDGAMMA,INDBETA) * 
     &           (FEF(INDI,INDALPHA,INDKAPPA,INDA) * 
     &           FE(INDJ,INDGAMMA) + FE(INDI,INDALPHA) * 
     &           FEF(INDJ,INDGAMMA,INDKAPPA,INDA))
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        K22=IM / (IM - (INE - 3.0D0)) * TEMPK22
        GCONST=-MUNEQ * DELTAT / ETAN
        RES=GCONST * (K21 + K22)    
        !CALL TENSOR2MATRIX22(res,RES99)
        !DO i=1,9
        !    PRINT *,(RES99(i,j),j=1,9)
        !ENDDO
      END SUBROUTINE KGF 

!==================== KGd ==================== 
      SUBROUTINE KGD(PROPS,F,D,FV,DELTAT,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,
     &     GINE,INVL,KROND
        REAL(KIND=8),DIMENSION(3) :: D,D0,INED
        REAL(KIND=8),DIMENSION(3,3,3,3) :: GINVL
        REAL(KIND=8),DIMENSION(3,3,3) :: INVLD,K1,K2,RES
        REAL(KIND=8),DIMENSION(9,9) :: RES99
        REAL(KIND=8),DIMENSION(9,3) :: RES93
        REAL(KIND=8) :: DELTAT,LPA,LPER,IN0,INE,GCONST,ETAN,MUNEQ
        INTEGER :: INDA,INDB,INDL,INDI,INDALPHA,INDBETA,INDJ
        INTEGER :: INFO
        MUNEQ =PROPS(5)
        ETAN=PROPS(6)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL PARGPARINE(INE,L0,FE,INVL,PROPS,GINE)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        CALL PARINEPARD(D,LPA,LPER,L0,FE,KROND,INED)
        CALL PARGPARINVL(INE,L0,FE,PROPS,GINVL)
        CALL PARINVLPARD(D,LPA,LPER,KROND,INVLD)
        !CALL TENSOR2MATRIX21(INVLD,RES93)
        !THE FIRST TERM
        K1=0.0D0
        DO INDA=1,3
        DO INDB=1,3
        DO INDL=1,3
        K1(INDA,INDB,INDL)=K1(INDA,INDB,INDL) + 
     &    GINE(INDA,INDB) * INED(INDL)
        ENDDO
        ENDDO
        ENDDO
        !THE SECOND TERM
        K2=0.0D0
        DO INDALPHA=1,3
         DO INDBETA=1,3
          DO INDL=1,3
           DO INDI=1,3
            DO INDJ=1,3
        K2(INDALPHA,INDBETA,INDL)=K2(INDALPHA,INDBETA,INDL) + 
     &       GINVL(INDALPHA,INDBETA,INDI,INDJ) * 
     &        INVLD(INDI,INDJ,INDL)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        GCONST=-MUNEQ*DELTAT/ETAN
        RES=GCONST*(K1+K2)
        !CALL TENSOR2MATRIX21(RES,RES93)
        !DO i=1,9
        !    PRINT '(9F10.5)',(RES93(i,j),j=1,3)
        !ENDDO
        END SUBROUTINE KGD

        
!==================== partial invl partial d ====================
      SUBROUTINE PARINVLPARD(D,LPA,LPER,KROND,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: KROND 
        REAL(KIND=8),DIMENSION(3) :: D
        REAL(KIND=8) :: LPA,LPER,FACTOR
        REAL(KIND=8),DIMENSION(3,3,3) :: RES
        INTEGER :: INDI,INDK,INDJ
        Factor=(1.0D0/LPA - 1.0D0/LPER)
        RES=0.0D0
        DO INDI=1,3
         DO INDK=1,3
          DO INDJ=1,3
        Res(INDI,INDK,INDJ)=Res(INDI,INDK,INDJ) 
     &       + Factor * (KronD(INDI,INDJ)*D(INDK) + 
     &         KronD(INDK,INDJ)*D(INDI))
          ENDDO
         ENDDO
        ENDDO    
      END SUBROUTINE PARINVLPARD
    
!==================== partial g partial invl ====================
      SUBROUTINE PARGPARINVL(INE,L0,FE,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: INE,IM 
        REAL(KIND=8),DIMENSION(3,3) :: L0,FE
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,TEMP
        INTEGER ::INDI,INDALPHA,INDJ,INDGAMMA,INDBETA
        IM=PROPS(2)
        TEMP=0.0D0
        DO INDI=1,3
         DO INDALPHA=1,3
          DO INDJ=1,3
           DO INDGAMMA=1,3
            DO INDBETA=1,3
        TEMP(INDALPHA,INDBETA,INDI,INDJ)=
     &     TEMP(INDALPHA,INDBETA,INDI,INDJ) 
     &     + FE(INDI,INDALPHA)*FE(INDJ,INDGAMMA)*L0(INDGAMMA,INDBETA)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        RES=IM / (IM - (INe - 3.0D0))*TEMP;
      END SUBROUTINE PARGPARINVL    
      
    !==================== KFd ====================
      SUBROUTINE KFD(PROPS,FN,FVN,DN,F,FV,D,DELTAT,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,IN0,INE,DELTAT,
     &      ETAD,FCONST,LPA,LPER,MEQCONST,MNEQCONST,DCONST 
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV,
     &   CAUSTRESS,L0,L,INVF,INVL
        REAL(KIND=8),DIMENSION(3,3) :: G,GE,FE,KROND
        REAL(KIND=8),DIMENSION(3) :: DN,D0,D,MEQ,MNEQ
        REAL(KIND=8),DIMENSION(3) :: DPENALTY
        REAL(KIND=8) :: TEMP0
        REAL(KIND=8),DIMENSION(3) :: TEMP1,IND,MEQIN,INED,MNEQINE
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2 
        REAL(KIND=8),DIMENSION(3,3) :: RES
        REAL(KIND=8),DIMENSION(3,3) :: K1,K21,K22,K2,K31,
     &    K32,K3,KD
        INTEGER :: INFO
        INTEGER :: INDI,INDJ,INDA,INDM,INDL
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        DPEN=PROPS(11)
        CALL SIGMA(PROPS,FN,FVN,DN,D0,CAUSTRESS)
        CALL FUNCETAD(PROPS,CAUSTRESS,ETAD)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        CALL MATRIX_INVERSE(F,INVF,3,INFO)
        FCONST=DELTAT/ETAD
        !THE FIRST TERM
        K1=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDA=1,3
        K1(INDI,INDJ)=K1(INDI,INDJ) - 0.5D0 * (INVF(INDA,INDI) * 
     &     FN(INDJ,INDA) - FN(INDI,INDA) * INVF(INDA,INDJ))
          ENDDO
         ENDDO
        ENDDO
        K1=KROND + K1
        !THE SECOND TERM
        K21=0.0D0
        CALL PARINPARD(D,LPA,LPER,L0,F,KROND,IND)
        CALL PARMEQPARIN(D,IN0,L0,F,PROPS,MEQIN)
        DO INDI=1,3
         DO INDJ=1,3
        K21(INDI,INDJ)=K21(INDI,INDJ) + MEQIN(INDI) * IND(INDJ)
         ENDDO
        ENDDO
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        MEQCONST=-MUEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &  IM/(IM-(IN0-3.0D0))
        TEMP2=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDM=1,3
           DO INDL=1,3
        TEMP2(INDI,INDJ)=TEMP2(INDI,INDJ) + G(INDL,INDM)*( 
     &    KROND(INDL,INDJ) * D(INDM) * D(INDI) + 
     &    KROND(INDM,INDJ) * D(INDL) * D(INDI) + 
     &    KROND(INDI,INDJ) * D(INDL) * D(INDM) )
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        K22=MEQCONST*(G-TEMP2)
        K2=K21+K22
        ! THE THIRD TERM
        K31=0.0D0
        CALL PARINEPARD(D,LPA,LPER,L0,FE,KROND,INED)
        CALL PARMNEQPARINE(D,INE,L0,FE,PROPS,MNEQINE)
        DO INDI=1,3
         DO INDJ=1,3
        K31(INDI,INDJ)=K31(INDI,INDJ) + MNEQINE(INDI) * INED(INDJ)
         ENDDO
        ENDDO
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        MNEQCONST=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &   IM/(IM-(INE-3.0D0))
        TEMP2=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDM=1,3
           DO INDL=1,3
        TEMP2(INDI,INDJ)=TEMP2(INDI,INDJ) + GE(INDL,INDM) * ( 
     &  KROND(INDL,INDJ) * D(INDM) * D(INDI) + 
     &  KROND(INDM,INDJ) * D(INDL) * D(INDI) + 
     &  KROND(INDI,INDJ) * D(INDL) * D(INDM) )
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        K32=MNEQCONST*(G-TEMP2)
        K3=K31+K32
        ! THE FOURTH TERM 
        CALL INNER_PRODUCT(D,D,TEMP0)
        DCONST=DPEN*DELTAT/ETAD*(TEMP0 - 1.0D0) 
        KD=0.0D0  ! Initialize the KD matrix to zero
        DO INDI=1,3
         DO INDJ=1,3
        KD(INDI,INDJ)=KD(INDI,INDJ) + 4.0D0*D(INDI)*D(INDJ) 
     &       + (TEMP0 - 1.0D0) * KROND(INDI,INDJ)
         ENDDO
        ENDDO
        KD=DCONST*KD
        RES=K1 + FCONST*(K2+K3)-KD
      END SUBROUTINE KFD
      
!==================== Numerical NKFD ====================  
      SUBROUTINE NKFD(PROPS,FN,FVN,DN,F,FV,D,DELTAT,D0,RES)
      IMPLICIT NONE 
      REAL(KIND=8),DIMENSION(12) :: PROPS
      REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV
      REAL(KIND=8),DIMENSION(3) :: DN,D0, D
      COMPLEX,DIMENSION(3) :: pertF, pertD
      REAL(KIND=8),DIMENSION(3,3) :: RES
      REAL(KIND=8) :: h,DELTAT
      COMPLEX :: eps, Det_pertF
      INTEGER :: INDI,INDJ,ISTAT
      
      h = 1E-10
      eps = CMPLX(0.0, h, kind=kind(h))
      RES = 0.0
      
      DO INDI = 1,3
        DO INDJ = 1,3
        pertD = D
        pertD(INDJ) = D(INDJ) + eps
        CALL ComplexFUNCF(PROPS,FN,DN,FVN,F,FV,pertD,DELTAT,D0,pertF)
        RES(INDI,INDJ) = aimag(pertF(INDI))/h
        END DO
      END DO

      END SUBROUTINE NKFD
      
!==================== PARTIAL F PARTIAL F ====================    
      SUBROUTINE PARFPARF(RES)
      IMPLICIT NONE
            REAL(KIND=8),DIMENSION(3,3,3,3) :: RES
            REAL(KIND=8),DIMENSION(3,3) :: KROND
            INTEGER :: INDA,INDJ,INDI,INDB
            CALL CREATE_IDENTITY_MATRIX(3,KROND)
            RES=0.0D0
            DO INDI=1,3
                DO INDA=1,3
                    DO INDJ=1,3
                        DO INDB=1,3
                            RES(INDI,INDA,INDJ,INDB)=
     &                        RES(INDI,INDA,INDJ,INDB) + 
     &                           KROND(INDI,INDJ) * KROND(INDA,INDB)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
      END SUBROUTINE PARFPARF

!==================== PARTIAL INVF PARTIAL F ====================    
      SUBROUTINE PARINVFPARF(F,RES)
      IMPLICIT NONE
            REAL(KIND=8),DIMENSION(3,3) :: F,INVF
            REAL(KIND=8),DIMENSION(3,3,3,3) :: RES
            INTEGER :: INFO
            INTEGER :: INDA,INDJ,INDK,INDB
            CALL MATRIX_INVERSE(F,INVF,3,INFO)
            RES=0.0D0
            DO INDA=1,3
             DO INDJ=1,3
              DO INDK=1,3
               DO INDB=1,3
            RES(INDA,INDJ,INDK,INDB)=RES(INDA,INDJ,INDK,INDB) - 
     &         INVF(INDA,INDK) * INVF(INDB,INDJ)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
      END SUBROUTINE PARINVFPARF    
      
!==================== KGV ====================
      SUBROUTINE KGV(PROPS,FVN,F,D,FV,DELTAT,D0,RES) 
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(3,3) :: FVN,F,FV,FE ,L0,L,
     &    INVFV,INVL,GINE,INEFV
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,IN0,INE,DELTAT,
     &   ETAD,FCONST,LPA,LPER,MEQCONST,MNEQCONST,GCONST,ETAN 
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,INVFVFV,K1,K21,
     &   FEFV,TEMP4,K22
        REAL(KIND=8),DIMENSION(9,9) :: RES99
        INTEGER :: INFO
        INTEGER :: INDALPHA,INDBETA,INDKAPPA,INDA,INDB,
     &   INDI,INDGAMMA,INDJ
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        ETAN=PROPS(6)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL MATRIX_INVERSE(FV,INVFV,3,INFO)
        !THE FIRST TERM
        K1=0.0D0
        CALL PARINVFPARF(FV,INVFVFV)
        DO INDALPHA=1,3
         DO INDBETA=1,3
          DO INDKAPPA=1,3
           DO INDA=1,3
            DO INDB=1,3
        K1(INDALPHA,INDBETA,INDKAPPA,INDA)=
     &    K1(INDALPHA,INDBETA,INDKAPPA,INDA) - 
     &     FVN(INDALPHA,INDB) * INVFVFV(INDB,INDBETA,INDKAPPA,INDA)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        !THE SECOND TERM
        K21=0.0D0
        CALL PARGPARINE(INE,L0,FE,INVL,PROPS,GINE)
        CALL PARINEPARFV(L0,INVL,FE,FV,INEFV)
        DO INDKAPPA=1,3
         DO INDA=1,3
          DO INDALPHA=1,3
           DO INDBETA=1,3
        K21(INDALPHA,INDBETA,INDKAPPA,INDA)=
     &   K21(INDALPHA,INDBETA,INDKAPPA,INDA) + 
     &   INEFV(INDALPHA,INDBETA) * GINE(INDKAPPA,INDA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL PARFEPARFV(FE,FV,FEFV)
        TEMP4=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDGAMMA=1,3
           DO INDBETA=1,3
            DO INDKAPPA=1,3
             DO INDA=1,3
              DO INDALPHA=1,3
        TEMP4(INDALPHA,INDBETA,INDKAPPA,INDA)=
     &   TEMP4(INDALPHA,INDBETA,INDKAPPA,INDA) + 
     &    INVL(INDI,INDJ) * L0(INDGAMMA,INDBETA) * 
     &    (FEFV(INDI,INDALPHA,INDKAPPA,INDA) * FE(INDJ,INDGAMMA) + 
     &     FE(INDI,INDALPHA) * FEFV(INDJ,INDGAMMA,INDKAPPA,INDA))
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        K22=IM/(IM-INE+3.0D0)*TEMP4
        GCONST=-MUNEQ*DELTAT/ETAN
        RES=K1 + GCONST*(K21+K22)
      END SUBROUTINE KGV      
        
!==================== PARITAL G PARTIAL INE ====================  
      SUBROUTINE PARGPARINE(INE,L0,FE,INVL,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: IM,INE
        REAL(KIND=8),DIMENSION(3,3) :: FE,L0,INVL,TEMP,RES
        INTEGER :: INFO,INDALPHA,INDBETA,INDI,INDJ,INDGAMMA
        IM=PROPS(2)
        TEMP=0.0D0
        DO INDALPHA=1,3
         DO INDBETA=1,3
           DO INDI=1,3
            DO INDJ=1,3
             DO INDGAMMA=1,3
        TEMP(INDALPHA,INDBETA)=TEMP(INDALPHA,INDBETA) + 
     &         FE(INDI,INDALPHA) * INVL(INDI,INDJ) * 
     &   FE(INDJ,INDGAMMA) * L0(INDGAMMA,INDBETA)
             ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        RES=(IM / (IM - (INE - 3.0D0))**2.0D0) * TEMP
      END SUBROUTINE PARGPARINE        

      
!==================== KFF ====================
      SUBROUTINE KFF(PROPS,FN,FVN,DN,F,FV,D,DELTAT,D0,RES)    
      IMPLICIT NONE
      REAL(KIND=8),DIMENSION(12) :: PROPS
      REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,IN0,INE,DELTAT,ETAD,
     &   FCONST,LPA,LPER,MEQCONST,MNEQCONST,DCONST 
      REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV,CAUSTRESS,
     &   L0,L,INVF,INVL
      REAL(KIND=8),DIMENSION(3,3) :: G,GE,FE,KROND
      REAL(KIND=8),DIMENSION(3) :: DN,D0,D,MEQ,MNEQ
      REAL(KIND=8),DIMENSION(3) :: DPENALTY,FETAD
      REAL(KIND=8) :: TEMP0
      REAL(KIND=8),DIMENSION(3) :: TEMP1,IND,MEQIN,INED,MNEQINE
      REAL(KIND=8),DIMENSION(3,3) :: TEMP2,INEFV,INF,INEF,
     &  ETADSIGMA 
      REAL(KIND=8),DIMENSION(3,3,3) :: RES,MNEQGE,MEQG
      REAL(KIND=8),DIMENSION(3,3,3) :: K1,K21,K22,K2,K31,
     & K32,K3,K4
      REAL(KIND=8),DIMENSION(3,3,3,3) :: PARFF,GF,GEF,SIGMAF
      REAL(KIND=8),DIMENSION(3,9) :: RES39
      REAL(KIND=8),DIMENSION(9,9) :: RES99
      INTEGER :: INFO
      INTEGER :: INDI,INDKAPPA,INDA,INDM,INDL,INDN,INDB,INDJ
      Q=PROPS(1)
      CALL SIGMA(PROPS,FN,FVN,DN,D0,CAUSTRESS)
      CALL FUNCETAD(PROPS,CAUSTRESS,ETAD)
      CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
      FCONST=DELTAT/ETAD
      CALL MATRIX_INVERSE(F,INVF,3,INFO)
!THE FIRST TERM
      K1=0.0D0
      DO INDI=1,3
          DO INDKAPPA=1,3
              DO INDA=1,3
                  DO INDB=1,3
                      DO INDJ=1,3
                          K1(INDI,INDKAPPA,INDA)=
     & K1(INDI,INDKAPPA,INDA) + 0.5D0 * 
     & (INVF(INDB,INDKAPPA) * INVF(INDA,INDI) * FN(INDJ,INDB) - 
     & FN(INDI,INDB) * INVF(INDA,INDJ) * INVF(INDB,INDKAPPA)) * D(INDJ)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

!THE SECOND TERM
      CALL TENSOR2MATRIX12(K1,RES39)
      CALL PARFPARF(PARFF)
      CALL PARINPARF(L0,INVL,F,PARFF,INF)
      CALL PARMEQPARIN(D,IN0,L0,F,PROPS,MEQIN)
      K21=0.0D0
      DO INDI=1,3
          DO INDKAPPA=1,3
              DO INDA=1,3
                  K21(INDI,INDKAPPA,INDA)=K21(INDI,INDKAPPA,INDA)
     &            + MEQIN(INDI) * INF(INDKAPPA,INDA)
              ENDDO
          ENDDO
      ENDDO
      

      K22=0.0D0
      CALL PARMEQPARG(D,IN0,PROPS,MEQG)
      CALL PARGPARF(L0,F,GF)
      DO INDI=1,3
          DO INDKAPPA=1,3
              DO INDA=1,3
                  DO INDM=1,3
                      DO INDN=1,3
                          K22(INDI,INDKAPPA,INDA)=
     &                     K22(INDI,INDKAPPA,INDA) + 
     &                      MEQG(INDI,INDM,INDN) * 
     &                     GF(INDM,INDN,INDKAPPA,INDA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      K2=FCONST * (K21 + K22)

!THE THIRD TERM
      CALL PARINEPARF(L0,INVL,FE,FV,INEF)
      CALL PARMNEQPARINE(D,INE,L0,FE,PROPS,MNEQINE)
      K31=0.0D0
      DO INDI=1,3
          DO INDKAPPA=1,3
              DO INDA=1,3
                  K31(INDI,INDKAPPA,INDA)=K31(INDI,INDKAPPA,INDA)
     &            + MNEQINE(INDI) * INEF(INDKAPPA,INDA)
              ENDDO
          ENDDO
      ENDDO
            

      
      K32=0.0D0
      CALL PARMNEQPARGE(D,INE,PROPS ,MNEQGE)
      CALL PARGEPARF(L0,FE,FV,GEF)
      DO INDI=1,3
          DO INDKAPPA=1,3
              DO INDA=1,3
                  DO INDM=1,3
                      DO INDN=1,3
                          K32(INDI,INDKAPPA,INDA)=
     &                    K32(INDI,INDKAPPA,INDA) + 
     &                  MNEQGE(INDI,INDM,INDN) * 
     &                  GEF(INDM,INDN,INDKAPPA,INDA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      K3=FCONST*(K31 + K32)
      

!THE FOURTH TERM
      CALL PARFPARETAD(D,IN0,L0,F,INE,FE,PROPS,ETAD,
     &  FCONST,FETAD)
      CALL PARETADPARSIGMA(PROPS,CAUSTRESS,ETAD,ETADSIGMA)
      CALL PARSIGMAPARF(PROPS,F,FV,D,D0,SIGMAF)
      K4=0.0D0
      DO INDI=1,3
          DO INDKAPPA=1,3
              DO INDA=1,3
                  DO INDM=1,3
                      DO INDN=1,3
                          K4(INDI,INDKAPPA,INDA)=
     &                     K4(INDI,INDKAPPA,INDA) + 
     &                      FETAD(INDI) * ETADSIGMA(INDM,INDN) *
     &                      SIGMAF(INDM,INDN,INDKAPPA,INDA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      RES=K1 + K2 + K3
      END SUBROUTINE KFF    
        
!==================== Partial P partial d ====================        
      SUBROUTINE PARPPARD(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
          REAL(KIND=8),DIMENSION(12) :: PROPS
          REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,SEQIN,
     &    INVL,INF,KROND,PEQCONST,SNEQINE,INEF,G,GE
          REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,
     &     LPA,LPER,IN0,INE,Det,CONST
          REAL(KIND=8),DIMENSION(3) :: D,IND,D0,INED
          REAL(KIND=8),DIMENSION(3,3,3) :: K11,RES,TEMP,
     &     K12,SEQD,SNEQD,SIGMAD
          REAL(KIND=8),DIMENSION(9,3) :: RES93
          INTEGER :: INFO
          INTEGER :: INDI,INDJ,INDKAPPA,INDBETA
          CALL PARSIGMAPARD(PROPS,F,FV,D,D0,SIGMAD)
          CALL MATRIX_INVERSE(F,INVF,3,INFO)
          CALL MDet(F,Det)
          RES=0.0D0
          DO INDI=1,3
              DO INDJ=1,3
                  DO INDKAPPA=1,3
                      DO INDBETA=1,3
                          RES(INDI,INDJ,INDKAPPA)=
     &                    RES(INDI,INDJ,INDKAPPA) + 
     &                      SIGMAD(INDI,INDBETA,INDKAPPA) 
     &                       * INVF(INDJ,INDBETA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
          RES=Det * RES
    !CALL TENSOR2MATRIX21(RES,RES93)
    !DO INDI=1,9
    !    PRINT *,(RES93(INDI,INDJ),INDJ=1,3)
    !END 
      END SUBROUTINE PARPPARD

!==================== Partial sigma partial d ====================        
      SUBROUTINE PARSIGMAPARD(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
      REAL(KIND=8),DIMENSION(12) :: PROPS
      REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,SEQIN,
     &  INVL,INF,KROND,PEQCONST,SNEQINE,INEF,G,GE
      REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,INE,
     &   Det,CONST
      REAL(KIND=8),DIMENSION(3) :: D,IND,D0,INED
      REAL(KIND=8),DIMENSION(3,3,3) :: K11,RES,TEMP,K12,
     & SEQD,SNEQD
      REAL(KIND=8),DIMENSION(9,3) :: RES93
      INTEGER :: INFO
      INTEGER :: INDI,INDJ,INDL,INDK
      CALL CREATE_IDENTITY_MATRIX(3,KROND)
      Q=PROPS(1)
      IM=PROPS(2)
      MUEQ=PROPS(4)
      MUNEQ=PROPS(5)
      CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
!THE FIRST TERM
      CALL PARSEQPARIN(PROPS,F,FV,D,D0,SEQIN)
      CALL PARINPARD(D,LPA,LPER,L0,F,KROND,IND)
      K11=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDL=1,3
                  K11(INDI,INDJ,INDL)=K11(INDI,INDJ,INDL) 
     &             + SEQIN(INDI,INDJ) * IND(INDL)
              ENDDO
          ENDDO
      ENDDO
      G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
      CONST=MUEQ / (1.0D0 - Q) * IM / (IM - (IN0 - 3.0D0))
     &  * (-3.0D0 * Q / (2.0D0 * (1.0D0 + 2.0D0 * Q)))
      TEMP=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDL=1,3
                  DO INDK=1,3
                      TEMP(INDI,INDJ,INDL)=TEMP(INDI,INDJ,INDL) 
     &                 + KROND(INDI,INDL) * G(INDJ,INDK) * D(INDK) + 
     &                 KROND(INDK,INDL) * G(INDJ,INDK) * D(INDI) + 
     &                     KROND(INDJ,INDL) * G(INDI,INDK) * D(INDK) + 
     &                     KROND(INDK,INDL) * G(INDI,INDK) * D(INDJ)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      K12=CONST*TEMP
      SEQD=K11 + K12
!THE SECOND TERM
      CALL PARSNEQPARINE(PROPS,F,FV,D,D0,SNEQINE)
      CALL PARINEPARD(D,LPA,LPER,L0,FE,KROND,INED)
      K11=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDL=1,3
                  K11(INDI,INDJ,INDL)=K11(INDI,INDJ,INDL)
     &            + SNEQINE(INDI,INDJ) * INED(INDL)
              ENDDO
          ENDDO
      ENDDO
      GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
      CONST=MUNEQ / (1.0D0 - Q) * IM / (IM - (INE - 3.0D0)) 
     & * (-3.0D0 * Q / (2.0D0 * (1.0D0 + 2.0D0 * Q)))
      TEMP=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDL=1,3
                  DO INDK=1,3
                      TEMP(INDI,INDJ,INDL)=TEMP(INDI,INDJ,INDL) +
     &                KROND(INDI,INDL) * GE(INDJ,INDK) * D(INDK) + 
     &                KROND(INDK,INDL) * GE(INDJ,INDK) * D(INDI) + 
     &                KROND(INDJ,INDL) * GE(INDI,INDK) * D(INDK) + 
     &                KROND(INDK,INDL) * GE(INDI,INDK) * D(INDJ)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      K12=CONST*TEMP
      SNEQD=K11 + K12
      RES=SEQD + SNEQD
!CALL TENSOR2MATRIX21(RES,RES93)
!DO INDI=1,9
!    PRINT *,(RES93(INDI,INDJ),INDJ=1,3)
!END 
      END SUBROUTINE PARSIGMAPARD

!==================== Partial p partial FV ====================        
      SUBROUTINE PARPPARFV(PROPS,F,FV,D,D0,RES)    
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,SEQIN,
     &      INVL,INF,KROND,PEQCONST,SNEQINE,INEFV
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,PARFF,K1,GF,SEQG,
     &   K12,K2,SEQF,SNEQGE,GEFV,SNEQF,SIGMAFV
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(9,9) :: RES99
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,
     &   INE,Det
        INTEGER :: INFO   
        INTEGER :: INDI,INDJ,INDKAPPA,INDBETA,INDALPHA
        CALL PARSIGMAPARFV(PROPS,F,FV,D,D0,SIGMAFV)
        CALL MATRIX_INVERSE(F,INVF,3,INFO)
        CALL MDet(F,Det)
        RES=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDALPHA=1,3
           DO INDBETA=1,3
            DO INDKAPPA=1,3
            RES(INDI,INDJ,INDALPHA,INDBETA)=
     &       RES(INDI,INDJ,INDALPHA,INDBETA) + 
     &        SIGMAFV(INDI,INDKAPPA,INDALPHA,INDBETA) 
     &       * INVF(INDJ,INDKAPPA)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        RES=Det*RES
        !CALL TENSOR2MATRIX22(RES,RES99)
        !DO INDI=1,9
        !    PRINT *,(RES99(INDI,INDJ),INDJ=1,9)
        !ENDDO
      END SUBROUTINE PARPPARFV

!==================== Partial sigma partial FV ====================        
      SUBROUTINE PARSIGMAPARFV(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
      REAL(KIND=8),DIMENSION(12) :: PROPS
      REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,SEQIN,
     &    INVL,INF,KROND,PEQCONST,SNEQINE,INEFV
      REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,PARFF,K1,GF,SEQG,
     & K12,K2,SEQF,SNEQGE,GEFV,SNEQF
      REAL(KIND=8),DIMENSION(3) :: D,D0
      REAL(KIND=8),DIMENSION(9,9) :: RES99
      REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,INE,Det
      INTEGER :: INFO   
      INTEGER :: INDI,INDJ,INDL,INDA,INDM,INDN
      INTEGER :: INDALPHA
      Q=PROPS(1)
!THE FIRST TERM
      CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
      CALL PARSNEQPARINE(PROPS,F,FV,D,D0,SNEQINE)
      CALL PARINEPARFV(L0,INVL,FE,FV,INEFV)
      K1=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDALPHA=1,3
                  DO INDA=1,3
                      K1(INDI,INDJ,INDALPHA,INDA)=
     &                 K1(INDI,INDJ,INDALPHA,INDA) + 
     &                     SNEQINE(INDI,INDJ) * INEFV(INDALPHA,INDA)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
!THE SECOND TERM
      CALL PARGEPARFV(L0,FE,FV,GEFV)
      CALL PARSNEQPARGE(PROPS,F,FV,D,D0,SNEQGE)
      K2=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDL=1,3
                  DO INDA=1,3
                      DO INDM=1,3
                          DO INDN=1,3
                              K2(INDI,INDJ,INDL,INDA)=
     &                        K2(INDI,INDJ,INDL,INDA) + 
     &                         SNEQGE(INDI,INDJ,INDM,INDN) 
     &                         * GEFV(INDM,INDN,INDL,INDA)
                          ENDDO
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      RES=K1 + K2
!CALL TENSOR2MATRIX22(res,RES99)
!DO i=1,9
!    PRINT *,(RES99(i,j),j=1,9)
!END 
      END SUBROUTINE PARSIGMAPARFV    
      
!==================== Partial P partial F ====================        
      SUBROUTINE PARPPARF(PROPS,F,FV,D,D0,RES)   
      IMPLICIT NONE
      REAL(KIND=8),DIMENSION(12) :: PROPS
      REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,SEQIN,
     &    INVL,INF,KROND,PEQCONST,SNEQINE,INEF,CAUSTRESS,JF
      REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,PARFF,K1,GF,
     &  SEQG,K2,K3,SEQF,SNEQGE,GEF,SNEQF,SIGMAF,INVFF
      REAL(KIND=8),DIMENSION(3) :: D,D0
      REAL(KIND=8),DIMENSION(9,9) :: RES99
      REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,INE,Det
      INTEGER :: INFO
      INTEGER :: INDI,INDJ,INDKAPPA,INDBETA,INDALPHA
      CALL MATRIX_INVERSE(F,INVF,3,INFO)
      CALL PARSIGMAPARF(PROPS,F,FV,D,D0,SIGMAF)
      K1=0.0D0
      DO INDI=1,3
          DO INDJ=1,3
              DO INDALPHA=1,3
                  DO INDBETA=1,3
                      DO INDKAPPA=1,3
                          K1(INDI,INDJ,INDALPHA,INDBETA)=
     &                    K1(INDI,INDJ,INDALPHA,INDBETA) + 
     &                SIGMAF(INDI,INDKAPPA,INDALPHA,INDBETA)
     &                * INVF(INDJ,INDKAPPA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      K2=0.0D0
      CALL SIGMA(PROPS,F,FV,D,D0,CAUSTRESS)
      CALL PARINVFPARF(F,INVFF)
      DO INDI=1,3
          DO INDJ=1,3
              DO INDALPHA=1,3
                  DO INDBETA=1,3
                      DO INDKAPPA=1,3
                          K2(INDI,INDJ,INDALPHA,INDBETA)=
     &                    K2(INDI,INDJ,INDALPHA,INDBETA) + 
     &                         CAUSTRESS(INDI,INDKAPPA) * 
     &                         INVFF(INDJ,INDKAPPA,INDALPHA,INDBETA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      K3=0.0D0
      JF=0.0D0
      CALL MDet(F,Det)
      DO INDI=1,3
          DO INDJ=1,3
              JF(INDI,INDJ)=JF(INDI,INDJ) + Det * INVF(INDJ,INDI)
          ENDDO
      ENDDO
      DO INDI=1,3
          DO INDJ=1,3
              DO INDALPHA=1,3
                  DO INDBETA=1,3
                      DO INDKAPPA=1,3
                          K3(INDI,INDJ,INDALPHA,INDBETA)=
     &                     K3(INDI,INDJ,INDALPHA,INDBETA) + 
     &                         JF(INDALPHA,INDBETA) * 
     &               CAUSTRESS(INDI,INDKAPPA) * INVF(INDJ,INDKAPPA)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      RES=Det*(K1 + K2) + K3
!CALL TENSOR2MATRIX22(RES,RES99)
!DO i=1,9
!    PRINT *,(RES99(i,j),j=1,9)
!END 
      END SUBROUTINE PARPPARF          
      
      
!==================== Partial sigma partial F ====================        
      SUBROUTINE PARSIGMAPARF(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,SEQIN,
     &    INVL,INF,KROND,PEQCONST,SNEQINE,INEF
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,PARFF,K11,GF,
     &   SEQG,K12,K2,SEQF,SNEQGE,GEF,SNEQF
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(9,9) :: RES99
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,
     &   INE,Det
        INTEGER :: INFO
        INTEGER :: INDI,INDJ,INDL,INDA,INDM,INDN
        CALL MATRIX_INVERSE(F,INVF,3,INFO)
        Q=PROPS(1)
        MUEQ=PROPS(4)
        PPEN=PROPS(10)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL PARSEQPARIN(PROPS,F,FV,D,D0,SEQIN)
        CALL PARFPARF(PARFF)
        CALL PARINPARF(L0,INVL,F,PARFF,INF)
        !THE FIRST TERM
        K11=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDL=1,3
           DO INDA=1,3
            K11(INDI,INDJ,INDL,INDA)=K11(INDI,INDJ,INDL,INDA) +
     &       SEQIN(INDI,INDJ) * INF(INDL,INDA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL PARGPARF(L0,F,GF)
        CALL PARSEQPARG(PROPS,F,FV,D,D0,SEQG)
        K12=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDL=1,3
           DO INDA=1,3
            DO INDM=1,3
             DO INDN=1,3
            K12(INDI,INDJ,INDL,INDA)=K12(INDI,INDJ,INDL,INDA) +
     &       SEQG(INDI,INDJ,INDM,INDN) * GF(INDM,INDN,INDL,INDA)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        !THE SECOND TERM
        K2=0.0D0
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        PEQCONST=PPEN*MUEQ*KROND
        CALL MDet(F,Det)
        DO INDI=1,3
         DO INDJ=1,3
          DO INDL=1,3
           DO INDA=1,3
        K2(INDI,INDJ,INDL,INDA)=K2(INDI,INDJ,INDL,INDA) + 
     &                Det * INVF(INDL,INDA) * PEQCONST(INDI,INDJ)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        SEQF=K11 + K12 + K2
        !THE THIRD TERM
        CALL PARSNEQPARINE(PROPS,F,FV,D,D0,SNEQINE)
        CALL PARINEPARF(L0,INVL,FE,FV,INEF)
        K11=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDL=1,3
           DO INDA=1,3
            K11(INDI,INDJ,INDL,INDA)=K11(INDI,INDJ,INDL,INDA) + 
     &                         SNEQINE(INDI,INDJ) * INEF(INDL,INDA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL PARSNEQPARGE(PROPS,F,FV,D,D0,SNEQGE)
        CALL PARGEPARF(L0,FE,FV,GEF)
        K12=0.0D0
        DO INDI=1,3
        DO INDJ=1,3
        DO INDL=1,3
        DO INDA=1,3
        DO INDM=1,3
        DO INDN=1,3
            K12(INDI,INDJ,INDL,INDA)=K12(INDI,INDJ,INDL,INDA)
     &        + SNEQGE(INDI,INDJ,INDM,INDN) * GEF(INDM,INDN,INDL,INDA)
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        SNEQF=K11 + K12
        RES=SEQF + SNEQF
      END SUBROUTINE PARSIGMAPARF      
      
!==================== Partial Sneq partial Ge ====================         
      SUBROUTINE PARSNEQPARGE(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,L0,L,FE,KROND,INVL
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8) :: Q,IM,MUNEQ,LPA,LPER,IN0,INE,CONST
        REAL(KIND=8),DIMENSION(3,3,3,3) :: TEMP,RES
        INTEGER :: INDI,INDJ,INDM,INDN
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        TEMP=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDM=1,3
           DO INDN=1,3
        TEMP(INDI,INDJ,INDM,INDN)=KROND(INDI,INDM) * KROND(INDJ,INDN)-
     &       3.0D0 * Q / (2.0D0 * (1.0D0 + 2.0D0 * Q)) * ( 
     &       D(INDI) * KROND(INDJ,INDM) * D(INDN) + 
     &       D(INDJ) * KROND(INDI,INDM) * D(INDN) )
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CONST=MUNEQ / (1.0D0 - Q) * IM / (IM - (INE - 3.0D0))
        RES=CONST * TEMP
      END SUBROUTINE PARSNEQPARGE  
      
!==================== Partial Sneq partial INe ====================     
      SUBROUTINE PARSNEQPARINE(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,GE,
     &    RES,TEMP1,TEMP2,INVL
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,
     &   IN0,INE,CONST
        INTEGER :: INDI,INDJ,INDK
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        TEMP1=0.0D0
        DO INDI=1,3
        DO INDJ=1,3
        DO INDK=1,3
            TEMP1(INDI,INDJ)=TEMP1(INDI,INDJ) + 
     &       (3.0D0*Q) / (2.0D0*(1.0D0+2.0D0*Q)) 
     &       * (D(INDI)*GE(INDJ,INDK)*D(INDK)
     &       + D(INDJ)*GE(INDI,INDK)*D(INDK))
        ENDDO
        ENDDO
        ENDDO
        TEMP2=GE - TEMP1
        CONST=MUNEQ / (1.0D0 - Q) * IM / (IM - (INE - 3.0D0))**2.0D0
        RES=CONST * TEMP2
      END SUBROUTINE PARSNEQPARINE          
      
!==================== Partial Seq partial G ====================         
      SUBROUTINE PARSEQPARG(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,L0,L,FE,KROND,INVL
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8) :: Q,IM,MUEQ,LPA,LPER,IN0,INE,CONST
        REAL(KIND=8),DIMENSION(3,3,3,3) :: TEMP,RES
        INTEGER :: INDI,INDJ,INDM,INDN
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        TEMP=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDM=1,3
           DO INDN=1,3
        TEMP(INDI,INDJ,INDM,INDN)=KROND(INDI,INDM) * KROND(INDJ,INDN)-
     &      3.0D0*Q/(2.0D0 * (1.0D0 + 2.0D0 * Q))*( 
     &       D(INDI) * KROND(INDJ,INDM) * D(INDN) + 
     &       D(INDJ) * KROND(INDI,INDM) * D(INDN))
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CONST=MUEQ / (1.0D0 - Q) * IM / (IM - (IN0 - 3.0D0))
        RES=CONST * TEMP
      END SUBROUTINE PARSEQPARG      
      
!==================== Partial Seq partial IN ====================     
      SUBROUTINE PARSEQPARIN(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,FE,L0,L,INVF,
     &     G,RES,TEMP1,TEMP2,INVL
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,
     &   LPER,IN0,INE,CONST
        INTEGER :: INDI,INDJ,INDK
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        TEMP1=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDK=1,3
            TEMP1(INDI,INDJ)=TEMP1(INDI,INDJ) + 
     &       (3.0D0*Q) / (2.0D0*(1.0D0+2.0D0*Q))
     &         * (D(INDI)*G(INDJ,INDK)*D(INDK)
     &         + D(INDJ)*G(INDI,INDK)*D(INDK))
          ENDDO
         ENDDO
        ENDDO
        TEMP2=G - TEMP1
        CONST=MUEQ / (1.0D0 - Q) * IM / (IM - (IN0 - 3.0D0))**2.0D0
        RES=CONST * TEMP2
      END SUBROUTINE PARSEQPARIN      
      
!==================== Partial f partial etaD ====================    
      SUBROUTINE PARETADPARSIGMA(PROPS,CAUSTRESS,ETAD,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: ETAD,TEMP1,KS,FROBENIUS 
        REAL(KIND=8),DIMENSION(3,3) :: CAUSTRESS,TEMP2,RES
        KS=PROPS(8)
        FROBENIUS=sqrt(sum(CAUSTRESS**2.0D0))
        TEMP1=ETAD * (-1.0D0) / (2.0D0 * KS * FROBENIUS)
        TEMP2=2.0D0*CAUSTRESS
        RES=TEMP1*TEMP2
      END SUBROUTINE PARETADPARSIGMA      
      
!==================== Partial f partial etaD ====================    
      SUBROUTINE PARFPARETAD(D,IN0,L0,F,INE,FE,
     &     PROPS,ETAD,FCONST,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,IN0,INE,DELTAT,
     &    ETAD,FCONST,LPA,LPER
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV,
     &   CAUSTRESS,L0,L,INVF
        REAL(KIND=8),DIMENSION(3,3) :: G,GE,FE
        REAL(KIND=8),DIMENSION(3) :: DN,D0,D,MEQ,MNEQ
        REAL(KIND=8),DIMENSION(3) :: RES
        REAL(KIND=8) :: TEMP0
        REAL(KIND=8),DIMENSION(3) :: TEMP1 
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2 
        INTEGER :: INFO
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        ! CALCULATE G AND GE
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        ! CALCULATE MUEQ AND MUNEQ
        TEMP1=MATMUL(G,D)
        CALL INNER_PRODUCT(D,TEMP1,TEMP0)
        MEQ=-MUEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &  IM/(IM-(IN0-3.0D0))*(MATMUL(G,D)-TEMP0*D)
        TEMP1=MATMUL(GE,D)
        CALL INNER_PRODUCT(D,TEMP1,TEMP0)
        MNEQ=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &   IM/(IM-(INE-3.0D0))*(MATMUL(GE,D)-TEMP0*D)
        RES=FCONST*(-1.0D0/ETAD)*(MEQ + MNEQ)
      END SUBROUTINE PARFPARETAD      
      
!==================== KFV ====================
      SUBROUTINE KFV(PROPS,FN,FVN,DN,F,FV,D,DELTAT,D0,RES) 
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,IN0,INE,DELTAT,
     &   ETAD,FCONST,LPA,LPER,MEQCONST,MNEQCONST,DCONST 
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV,CAUSTRESS,
     &   L0,L,INVF,INVL
        REAL(KIND=8),DIMENSION(3,3) :: G,GE,FE,KROND
        REAL(KIND=8),DIMENSION(3) :: DN,D0,D,MEQ
        REAL(KIND=8),DIMENSION(3) :: DPENALTY
        REAL(KIND=8) :: TEMP0
        REAL(KIND=8),DIMENSION(3) :: TEMP1,IND,MEQIN,INED,MNEQINE
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2,INEFV 
        REAL(KIND=8),DIMENSION(3,3,3) :: RES,MNEQGE
        REAL(KIND=8),DIMENSION(3,3,3) :: K1,K2
        REAL(KIND=8),DIMENSION(3,3,3,3) :: GEFV
        REAL(KIND=8),DIMENSION(3,9) :: RES39
        INTEGER :: INFO
        INTEGER :: INDI,INDKAPPA,INDA,INDM,INDL,INDN
        Q=PROPS(1)
        CALL SIGMA(PROPS,FN,FVN,DN,D0,CAUSTRESS)
        CALL FUNCETAD(PROPS,CAUSTRESS,ETAD)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        FCONST=DELTAT/ETAD
        CALL PARMNEQPARINE(D,INE,L0,FE,PROPS,MNEQINE)
        CALL PARINEPARFV(L0,INVL,FE,FV,INEFV)
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        MNEQCONST=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)
     &    *IM/(IM-(INE-3.0D0))
        ! THE FIRST TERM
        K1=0.0D0
        DO INDI=1,3
         DO INDKAPPA=1,3
          DO INDA=1,3
        K1(INDI,INDKAPPA,INDA)=K1(INDI,INDKAPPA,INDA) + 
     &    MNEQINE(INDI) * INEFV(INDKAPPA,INDA)
          ENDDO
         ENDDO
        ENDDO

        ! THE SECOND TERM
        K2=0.0D0
        CALL PARMNEQPARGE(D,INE,PROPS,MNEQGE)
        CALL PARGEPARFV(L0,FE,FV,GEFV)
        DO INDI=1,3
         DO INDM=1,3
           DO INDN=1,3
            DO INDKAPPA=1,3
             DO INDA=1,3
        K2(INDI,INDKAPPA,INDA)=K2(INDI,INDKAPPA,INDA) 
     &    + MNEQGE(INDI,INDM,INDN) * GEFV(INDM,INDN,INDKAPPA,INDA)
             ENDDO
            ENDDO
           ENDDO
         ENDDO
        ENDDO
        RES=FCONST*(K1+K2)
      END SUBROUTINE KFV      
      
!==================== Partial G partial F ====================     
      SUBROUTINE PARGPARF(L0,F,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: L0,F
        REAL(KIND=8),DIMENSION(3,3,3,3) :: PARFF,RES
        INTEGER :: INDI,INDKAPPA,INDALPHA,INDJ,INDA,INDBETA
        CALL PARFPARF(PARFF)
        RES=0.0D0
        DO INDI=1,3
         DO INDKAPPA=1,3
          DO INDALPHA=1,3
           DO INDJ=1,3
            DO INDA=1,3
             DO INDBETA=1,3
        RES(INDI,INDJ,INDKAPPA,INDA)=
     &    RES(INDI,INDJ,INDKAPPA,INDA) + 
     &       L0(INDALPHA,INDBETA) * 
     &       (PARFF(INDI,INDALPHA,INDKAPPA,INDA) * F(INDJ,INDBETA) +
     &       F(INDI,INDALPHA) * PARFF(INDJ,INDBETA,INDKAPPA,INDA))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARGPARF      
      
!==================== Partial GE partial F ====================     
      SUBROUTINE PARGEPARF(L0,FE,FV,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: L0,FE,FV
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,FEF
        INTEGER :: INDI,INDJ,INDKAPPA,INDALPHA,INDA,INDB
        CALL PARFEPARF(FV,FEF)
        RES=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
          DO INDKAPPA=1,3
           DO INDALPHA=1,3
            DO INDA=1,3
             DO INDB=1,3
        RES(INDI,INDJ,INDKAPPA,INDALPHA)=
     &    RES(INDI,INDJ,INDKAPPA,INDALPHA) + 
     &     L0(INDA,INDB) * (FEF(INDI,INDA,INDKAPPA,INDALPHA)
     &   * FE(INDJ,INDB) + 
     &   FE(INDI,INDA) * FEF(INDJ,INDB,INDKAPPA,INDALPHA))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARGEPARF      
      
!==================== Partial GE partial FV ====================     
      SUBROUTINE PARGEPARFV(L0,FE,FV,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: L0,FE,FV
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES,FEFV
        INTEGER :: INDI,INDJ,INDKAPPA,INDALPHA,INDA,INDBETA
        CALL PARFEPARFV(FE,FV,FEFV)
        RES=0.0D0
        DO INDI=1,3
         DO INDKAPPA=1,3
          DO INDALPHA=1,3
           DO INDJ=1,3
            DO INDA=1,3
             DO INDBETA=1,3
        RES(INDI,INDJ,INDKAPPA,INDA)=
     &    RES(INDI,INDJ,INDKAPPA,INDA) 
     &   + L0(INDALPHA,INDBETA) * (FEFV(INDI,INDALPHA,INDKAPPA,INDA)
     &   * FE(INDJ,INDBETA) 
     &   + FE(INDI,INDALPHA) * FEFV(INDJ,INDBETA,INDKAPPA,INDA))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO    
      END SUBROUTINE PARGEPARFV       
      
!==================== Partial MNEQ partial GE ==================== 
      SUBROUTINE PARMNEQPARGE(D,INE,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3) :: D
        REAL(KIND=8),DIMENSION(3,3) :: KROND
        REAL(KIND=8),DIMENSION(3,3,3) :: RES,TEMP3
        REAL(KIND=8) :: INE,MNEQCONST,MUNEQ,Q,IM  
        INTEGER :: INDI,INDM,INDN
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        MNEQCONST=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)
     &    *IM/(IM-(INE-3.0D0))
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        TEMP3=0.0D0
        DO INDI=1,3
         DO INDM=1,3
          DO INDN=1,3
        TEMP3(INDI,INDM,INDN)=TEMP3(INDI,INDM,INDN) + 
     &    KROND(INDI,INDM)*D(INDN) - D(INDM)*D(INDN)*D(INDI)
          ENDDO
         ENDDO
        ENDDO
        RES=MNEQCONST * TEMP3
      END SUBROUTINE PARMNEQPARGE     
      
!==================== Partial INE partial F ====================  
      SUBROUTINE PARINEPARF(L0,INVL,FE,FV,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: FE,FV,L0,L,INVL,RES
        REAL(KIND=8),DIMENSION(3,3,3,3) :: FEF
        INTEGER :: INDI,INDJ,INDKAPPA,INDALPHA,INDBETA,INDA
        CALL PARFEPARF(FV,FEF)
        RES=0.0D0
        DO INDKAPPA=1,3
         DO INDA=1,3
          DO INDALPHA=1,3
           DO INDJ=1,3
            DO INDBETA=1,3
             DO INDI=1,3
        RES(INDKAPPA,INDA)=RES(INDKAPPA,INDA) + 
     &      L0(INDBETA,INDALPHA) * INVL(INDI,INDJ) * 
     &       (FEF(INDI,INDALPHA,INDKAPPA,INDA) * FE(INDJ,INDBETA) + 
     &       FE(INDI,INDALPHA) * FEF(INDJ,INDBETA,INDKAPPA,INDA))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARINEPARF     
      
!==================== Partial INE partial FV ====================  
      SUBROUTINE PARINEPARFV(L0,INVL,FE,FV,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: FE,FV,L0,L,INVL,RES
        REAL(KIND=8),DIMENSION(3,3,3,3) :: FEFV
        INTEGER :: INDI,INDJ,INDKAPPA,INDALPHA,INDBETA,INDA
        CALL PARFEPARFV(FE,FV,FEFV)
        RES=0.0D0
        DO INDKAPPA=1,3
         DO INDA=1,3
          DO INDI=1,3
           DO INDJ=1,3
            DO INDALPHA=1,3
             DO INDBETA=1,3
        RES(INDKAPPA,INDA)=RES(INDKAPPA,INDA) + 
     &   L0(INDBETA,INDALPHA) * INVL(INDI,INDJ) * 
     &  (FEFV(INDI,INDALPHA,INDKAPPA,INDA) * FE(INDJ,INDBETA) + 
     &  FE(INDI,INDALPHA) * FEFV(INDJ,INDBETA,INDKAPPA,INDA))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARINEPARFV      
      
!==================== Partial FE partial F ==================== 
      SUBROUTINE PARFEPARF(FV,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: FV,INVFV,KROND
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES
        INTEGER :: INDI,INDJ,INDB,INDALPHA
        INTEGER :: INFO
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        RES=0.0D0
        CALL MATRIX_INVERSE(FV,INVFV,3,INFO)
        DO INDI=1,3
         DO INDJ=1,3
          DO INDB=1,3
           DO INDALPHA=1,3
        RES(INDI,INDALPHA,INDJ,INDB)=
     &     RES(INDI,INDALPHA,INDJ,INDB)
     &     + KROND(INDI,INDJ)*INVFV(INDB,INDALPHA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARFEPARF      
      
!==================== Partial FE partial FV ==================== 
      SUBROUTINE PARFEPARFV(FE,FV,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: FE,FV,INVFV
        REAL(KIND=8),DIMENSION(3,3,3,3) :: RES
        INTEGER :: INDI,INDKAPPA,INDB,INDALPHA
        INTEGER :: INFO
        RES=0.0D0
        CALL MATRIX_INVERSE(FV,INVFV,3,INFO)
        DO INDI=1,3
         DO INDKAPPA=1,3
          DO INDB=1,3
           DO INDALPHA=1,3
        RES(INDI,INDALPHA,INDKAPPA,INDB)=
     &   RES(INDI,INDALPHA,INDKAPPA,INDB) - 
     &   FE(INDI,INDKAPPA)*INVFV(INDB,INDALPHA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARFEPARFV      
      
!==================== Partial INE partial d ====================  
      SUBROUTINE PARINEPARD(D,LPA,LPER,L0,FE,KROND,RES)
      IMPLICIT NONE
        REAL(KIND=8) :: LPA,LPER
        REAL(KIND=8),DIMENSION(3,3) :: L0,FE,KROND,INEINVL
        REAL(KIND=8),DIMENSION(3) :: D,RES,INED
        REAL(KIND=8),DIMENSION(3,3,3) :: INVLD
        INTEGER :: INDALPHA,INDBETA,INDI,INDK,INDJ
        INEINVL=0.0D0
        DO INDALPHA=1,3
         DO INDBETA=1,3
          DO INDI=1,3
           DO INDK=1,3
        INEINVL(INDI,INDK)=INEINVL(INDI,INDK) + 
     &   L0(INDALPHA,INDBETA) * FE(INDI,INDBETA) * FE(INDK,INDALPHA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        INVLD=0.0D0
        DO INDI=1,3
         DO INDK=1,3
          DO INDJ=1,3
        INVLD(INDI,INDK,INDJ)=INVLD(INDI,INDK,INDJ) + 
     &   (1.0D0/LPA - 1.0D0/LPER) * 
     &   (KROND(INDI,INDJ) * D(INDK) + D(INDI) * KROND(INDK,INDJ))
          ENDDO
         ENDDO
        ENDDO
        INED=0.0D0
        DO INDJ=1,3
         DO INDI=1,3
          DO INDK=1,3
        INED(INDJ)=INED(INDJ) + 
     &   INEINVL(INDI,INDK) * INVLD(INDI,INDK,INDJ)
          ENDDO
         ENDDO
        ENDDO
        RES=INED
      END SUBROUTINE PARINEPARD      
      
      
!==================== Partial IN partial d ====================     
      SUBROUTINE PARINPARD(D,LPA,LPER,L0,F,KROND,RES)
      IMPLICIT NONE
        REAL(KIND=8) :: LPA,LPER
        REAL(KIND=8),DIMENSION(3,3) :: L0,F,KROND,ININVL
        REAL(KIND=8),DIMENSION(3) :: D,RES,IND
        REAL(KIND=8),DIMENSION(3,3,3) :: INVLD
        INTEGER :: INDA,INDB,INDI,INDK,INDJ
        ININVL=0.0D0
        DO INDA=1,3
         DO INDB=1,3
          DO INDI=1,3
           DO INDK=1,3
            ININVL(INDI,INDK)=ININVL(INDI,INDK) + 
     &      L0(INDA,INDB) * F(INDI,INDB) * F(INDK,INDA)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        INVLD=0.0D0
        DO INDI=1,3
         DO INDK=1,3
          DO INDJ=1,3
        INVLD(INDI,INDK,INDJ)=INVLD(INDI,INDK,INDJ) + 
     &   (1.0D0/LPA - 1.0D0/LPER) * 
     &   (KROND(INDI,INDJ) * D(INDK) + D(INDI) * KROND(INDK,INDJ))
          ENDDO
         ENDDO
        ENDDO
        IND=0.0D0
        DO INDJ=1,3
         DO INDI=1,3
          DO INDK=1,3
        IND(INDJ)=IND(INDJ) + 
     &   ININVL(INDI,INDK) * INVLD(INDI,INDK,INDJ)
          ENDDO
         ENDDO
        ENDDO
        RES=IND
      END SUBROUTINE PARINPARD      
      
!==================== Partial IN partial F ====================  
      SUBROUTINE PARINPARF(L0,INVL,F,PARFF,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3) :: L0,F,INVL,RES
        REAL(KIND=8),DIMENSION(3,3,3,3) :: PARFF
        INTEGER :: INDI,INDKAPPA,INDA,INDB,INDJ,INDC
        RES=0.0D0
        DO INDKAPPA=1,3
         DO INDA=1,3
          DO INDB=1,3
           DO INDI=1,3
            DO INDJ=1,3
             DO INDC=1,3
        RES(INDKAPPA,INDA)=RES(INDKAPPA,INDA) + 
     &   L0(INDC,INDB) * INVL(INDI,INDJ) * ( 
     &   PARFF(INDI,INDB,INDKAPPA,INDA) * F(INDJ,INDC) + 
     &   F(INDI,INDB) * PARFF(INDJ,INDC,INDKAPPA,INDA))
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      END SUBROUTINE PARINPARF      
      
!==================== Partial MEQ partial G ====================  
      SUBROUTINE PARMEQPARG(D,IN0,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: IN0,IM,MUEQ,MEQCONST,Q
        REAL(KIND=8),DIMENSION(3) :: D
        REAL(KIND=8),DIMENSION(3,3,3) :: TEMP,RES
        REAL(KIND=8),DIMENSION(3,3) :: KROND
        INTEGER :: INDI,INDM,INDN
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MEQCONST=-MUEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)
     &    *IM/(IM-(IN0-3.0D0))
        TEMP=0.0D0
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        DO INDI=1,3
         DO INDM=1,3
          DO INDN=1,3
        TEMP(INDI,INDM,INDN)=TEMP(INDI,INDM,INDN) + 
     &   KROND(INDI,INDM) * D(INDN) - D(INDM) * D(INDN) * D(INDI)
          ENDDO
         ENDDO
        ENDDO
        RES=MEQCONST*TEMP
      END SUBROUTINE PARMEQPARG  
      
!==================== Partial MEQ partial IN ====================  
      SUBROUTINE PARMEQPARIN(D,IN0,L0,F,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: IN0,IM
        REAL(KIND=8),DIMENSION(3,3) :: L0,F
        REAL(KIND=8),DIMENSION(3) :: D,MEQ,RES
        IM=PROPS(2)
        CALL FUNCMEQ(D,IN0,L0,F,PROPS,MEQ) 
        RES=1.0D0 / (IM - IN0 + 3.0D0)
        RES=RES*MEQ  
      END SUBROUTINE PARMEQPARIN      
      
!==================== Partial MNEQ partial IN ====================  
      SUBROUTINE PARMNEQPARINE(D,INE,L0,FE,PROPS,RES)
      IMPLICIT NONE
            REAL(KIND=8),DIMENSION(12) :: PROPS
            REAL(KIND=8) :: INE,IM,Q,MUNEQ 
            REAL(KIND=8),DIMENSION(3,3) :: L0,FE
            REAL(KIND=8),DIMENSION(3) :: D,MNEQ,RES
            IM=PROPS(2)
            CALL FUNCMNEQ(D,INE,L0,FE,PROPS,MNEQ) 
            RES=1.0D0/(IM-INE+3.0D0)*MNEQ  
      END SUBROUTINE PARMNEQPARINE
      
      
!==================== Mneq ==================== 
      SUBROUTINE FUNCMNEQ(D,INE,L0,FE,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: INE,IM,Q,MUNEQ,TEMP0 
        REAL(KIND=8),DIMENSION(3,3) :: L0,FE,GE
        REAL(KIND=8),DIMENSION(3) :: D,TEMP1,MNEQ,RES
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        ! CALCULATE G AND GE
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        ! CALCULATE MUEQ AND MUNEQ
        TEMP1=MATMUL(GE,D)
        CALL INNER_PRODUCT(D,TEMP1,TEMP0)
        MNEQ=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &    IM/(IM-(INE-3.0D0))*(MATMUL(GE,D)-TEMP0*D)
        RES=MNEQ
      END SUBROUTINE FUNCMNEQ      
      
!==================== Meq ==================== 
      SUBROUTINE FUNCMEQ(D,IN0,L0,F,PROPS,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: IN0,IM,Q,MUEQ,TEMP0 
        REAL(KIND=8),DIMENSION(3,3) :: L0,F,G
        REAL(KIND=8),DIMENSION(3) :: D,TEMP1,MEQ,RES
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        ! CALCULATE G AND GE
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        ! CALCULATE MUEQ AND MUNEQ
        TEMP1=MATMUL(G,D)
        CALL INNER_PRODUCT(D,TEMP1,TEMP0)
        MEQ=-MUEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &    IM/(IM-(IN0-3.0D0))*(MATMUL(G,D)-TEMP0*D)
        RES=MEQ
      END SUBROUTINE FUNCMEQ  
      
!==================== Matrix to vector ====================     
      SUBROUTINE MATRIX2VECTOR(A,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3,3),INTENT(IN) :: A
        REAL(KIND=8),DIMENSION(9),INTENT(OUT) :: RES
        RES(1)=A(1,1)
        RES(2)=A(2,2)
        RES(3)=A(3,3)
        RES(4)=A(2,3)
        RES(5)=A(3,1)
        RES(6)=A(1,2)
        RES(7)=A(3,2)
        RES(8)=A(1,3)
        RES(9)=A(2,1)
      END SUBROUTINE MATRIX2VECTOR    
    
!==================== Vector to matrix ====================     
      SUBROUTINE VECTOR2MATRIX(A,RES)
        IMPLICIT NONE
        REAL(KIND=8),DIMENSION(9),INTENT(IN) :: A
        REAL(KIND=8),DIMENSION(3,3),INTENT(OUT) :: RES
        ! Assign the elements of the input vector 'A' to the output matrix 'RES'
        RES(1,1)=A(1)
        RES(1,2)=A(6)
        RES(1,3)=A(8)
        RES(2,1)=A(9)
        RES(2,2)=A(2)
        RES(2,3)=A(4)
        RES(3,1)=A(5)
        RES(3,2)=A(7)
        RES(3,3)=A(3)
      END SUBROUTINE VECTOR2MATRIX   
      
      
!==================== Function f ====================  
      SUBROUTINE FUNCF(PROPS,FN,DN,FVN,F,FV,D,DELTAT,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,IN0,INE,DELTAT,
     &    ETAD,FCONST,LPA,LPER
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV,
     &   CAUSTRESS,L0,L,INVF,INVL
        REAL(KIND=8),DIMENSION(3,3) :: G,GE,FE
        REAL(KIND=8),DIMENSION(3) :: DN,D0,D,MEQ,MNEQ
        REAL(KIND=8),DIMENSION(3) :: RES,DPENALTY
        REAL(KIND=8) :: TEMP0
        REAL(KIND=8),DIMENSION(3) :: TEMP1 
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2 
        INTEGER :: INFO
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        DPEN=PROPS(11)
        CALL SIGMA(PROPS,FN,FVN,DN,D0,CAUSTRESS)
        CALL FUNCETAD(PROPS,CAUSTRESS,ETAD)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        ! CALCULATE G AND GE
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        ! CALCULATE MUEQ AND MUNEQ
        TEMP1=MATMUL(G,D)
        CALL INNER_PRODUCT(D,TEMP1,TEMP0)
        MEQ=-MUEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &   IM/(IM-(IN0-3.0D0))*(MATMUL(G,D)-TEMP0*D)
        TEMP1=MATMUL(GE,D)
        CALL INNER_PRODUCT(D,TEMP1,TEMP0)
        MNEQ=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &   IM/(IM-(INE-3.0D0))*(MATMUL(GE,D)-TEMP0*D)
        FCONST=DELTAT/ETAD
        CALL INNER_PRODUCT(D,D,TEMP0)
        DPENALTY=DPEN*FCONST*(TEMP0- 1.0D0)**2.0D0*D
        CALL MATRIX_INVERSE(F,INVF,3,INFO)
        TEMP1=0.5D0*MATMUL((MATMUL(TRANSPOSE(INVF),TRANSPOSE(FN))
     &   -MATMUL(FN,INVF)),D)
        RES=D-DN-TEMP1+FCONST*(MEQ+MNEQ)-DPENALTY
      END SUBROUTINE FUNCF       

!==================== Function f ====================  
      SUBROUTINE ComplexFUNCF(PROPS,FN,DN,FVN,F,FV,D,DELTAT,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,DPEN,DELTAT,
     &    ETAD,FCONST,LPA,LPER
        COMPLEX :: IN0,INE
        REAL(KIND=8),DIMENSION(3,3) :: FN,FVN,F,FV,
     &   CAUSTRESS,INVF,INVL
        COMPLEX,DIMENSION(3,3) :: L0,L
        REAL(KIND=8),DIMENSION(3,3) :: G,GE,FE
        REAL(KIND=8),DIMENSION(3) :: DN,D0
        Complex,DIMENSION(3) :: D,MEQ,MNEQ,RES
        REAL(KIND=8),DIMENSION(3) :: DPENALTY
        COMPLEX :: TEMP0
        COMPLEX,DIMENSION(3) :: TEMP1 
        INTEGER :: INFO
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        DPEN=PROPS(11)
        CALL SIGMA(PROPS,FN,FVN,DN,D0,CAUSTRESS)
        CALL FUNCETAD(PROPS,CAUSTRESS,ETAD)
        CALL ComplexIN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        ! CALCULATE G AND GE
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        ! CALCULATE MUEQ AND MUNEQ
        TEMP1=MATMUL(G,D)
        CALL ComplexINNER_PRODUCT(D,TEMP1,TEMP0)
        MEQ=-MUEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &   IM/(IM-(IN0-3.0D0))*(MATMUL(G,D)-TEMP0*D)
        TEMP1=MATMUL(GE,D)
        CALL ComplexINNER_PRODUCT(D,TEMP1,TEMP0)
        MNEQ=-MUNEQ*(3.0D0*Q)/(1.0D0-Q)/(1.0D0+2.0D0*Q)*
     &   IM/(IM-(INE-3.0D0))*(MATMUL(GE,D)-TEMP0*D)
        FCONST=DELTAT/ETAD
        CALL ComplexINNER_PRODUCT(D,D,TEMP0)
        DPENALTY=DPEN*FCONST*(TEMP0- 1.0D0)**2.0D0*D
        CALL MATRIX_INVERSE(F,INVF,3,INFO)
        TEMP1=0.5D0*MATMUL((MATMUL(TRANSPOSE(INVF),TRANSPOSE(FN))
     &   -MATMUL(FN,INVF)),D)
        RES=D-DN-TEMP1+FCONST*(MEQ+MNEQ)-DPENALTY
      END SUBROUTINE ComplexFUNCF            
      
!==================== Function g ====================
      SUBROUTINE FUNCG(PROPS,FVN,F,FV,D,DELTAT,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8) :: Q,IM,MUNEQ,ETAN,DELTAT,IN0,INE,LPA,LPER
        REAL(KIND=8),DIMENSION(3,3) :: FVN,F,FV,L0,L,FE,
     &   RES,INVFV,KROND,INVL
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2
        INTEGER :: INFO
        Q=PROPS(1)
        IM=PROPS(2)
        MUNEQ=PROPS(5)
        ETAN=PROPS(6)
        CALL IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,L0,L,FE,INVL)
        CALL MATRIX_INVERSE(FV,INVFV,3,INFO)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        TEMP2=MATMUL(MATMUL(MATMUL(TRANSPOSE(FE),INVL),FE),L0)
        RES=KROND-MATMUL(FVN,INVFV) - MUNEQ*DELTAT/ETAN*
     &  (IM/(IM-INE+3.0D0)*TEMP2-KROND)
      END SUBROUTINE FUNCG  
      
!==================== Function etaD ====================   
      SUBROUTINE FUNCETAD(PROPS,CAUSTRESS,RES) 
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: CAUSTRESS
        REAL(KIND=8) :: ETAD0,KS,SUM_S,RES
        INTEGER :: INDI,INDJ
        ETAD0=PROPS(7)
        KS=PROPS(8)
        SUM_S=0.0D0
        DO INDI=1,3
         DO INDJ=1,3
            SUM_S=SUM_S + CAUSTRESS(INDI,INDJ)**2.0D0
         ENDDO
        ENDDO
        RES=ETAD0 * EXP(-SQRT(SUM_S)/ks)
      END SUBROUTINE FUNCETAD    
        
!==================== IN AND INE ==================== 
      SUBROUTINE IN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,
     &   L0,L,FE,INVL)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,RES,INVFV,FE
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(3,3) :: L,L0,G,GE,SEQ,SNEQ,
     &   KROND,INVL
        REAL(KIND=8),DIMENSION(3) :: GD,GED
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,
     &   INE,PPEQ,PNEQ,Det
        INTEGER :: INFO
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2,TEMP332
        REAL(KIND=8),DIMENSION(3) :: TEMP1
        REAL(KIND=8) :: TEMP0,TEMP112,Det_FV
        INTEGER :: ISTAT
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        PPEN=PROPS(10)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        ! CALCULATE Fe=F * Fv_inv
        CALL MATINV3D(FV,invFV,Det_FV,ISTAT)
        FE = MATMUL(F,invFV)
        ! CALCULATE LPA AND LPER
        LPA=1.0D0 + 2.0D0*Q
        LPER=1.0D0 - Q
        ! CALCULATE L AND L0
        L=(LPA - LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + LPER * KROND
        L0=(LPA - LPER) * MATMUL(RESHAPE(D0,[3,1]),
     &   TRANSPOSE(RESHAPE(D0,[3,1]))) + LPER * KROND
        ! CALCULATE IN AND INE
        INVL = (1.0D0/LPA - 1.0D0/LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + 1.0D0/LPER * KROND
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(F)),INVL),F)
        CALL CALCULATE_TRACE(TEMP2,3,IN0)
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(FE)),INVL),FE)
        CALL CALCULATE_TRACE(TEMP2,3,INE)
      END SUBROUTINE IN0INE      
     
!==================== IN AND INE ==================== 
      SUBROUTINE ComplexIN0INE(PROPS,F,FV,D,D0,LPA,LPER,IN0,INE,
     &   L0,L,FE)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,RES,INVFV,FE
        REAL(KIND=8),DIMENSION(3) :: D0
        COMPLEX,DIMENSION(3) :: D
        REAL(KIND=8),DIMENSION(3,3) :: G,GE,SEQ,SNEQ,
     &   KROND
        COMPLEX,DIMENSION(3,3) :: L,L0,INVL
        COMPLEX :: DET_L,IN0,INE
        REAL(KIND=8),DIMENSION(3) :: GD,GED
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,
     &   PPEQ,PNEQ,Det
        INTEGER :: INFO
        COMPLEX,DIMENSION(3,3) :: TEMP2
        REAL(KIND=8),DIMENSION(3) :: TEMP1
        REAL(KIND=8) :: TEMP0,TEMP112,Det_FV
        INTEGER :: ISTAT
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        PPEN=PROPS(10)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        ! CALCULATE Fe=F * Fv_inv
        CALL MATINV3D(FV,invFV,Det_FV,ISTAT)
        FE = MATMUL(F,invFV)
        ! CALCULATE LPA AND LPER
        LPA=1.0D0 + 2.0D0*Q
        LPER=1.0D0 - Q
        ! CALCULATE L AND L0
        L=(LPA - LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + LPER * KROND
        L0=(LPA - LPER) * MATMUL(RESHAPE(D0,[3,1]),
     &   TRANSPOSE(RESHAPE(D0,[3,1]))) + LPER * KROND
        ! CALCULATE IN AND INE
        CALL ComplexMATINV3D(L,INVL,Det_L,ISTAT)
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(F)),INVL),F)
        CALL COMPLEXCALCULATE_TRACE(TEMP2,3,IN0)
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(FE)),INVL),FE)
        CALL COMPLEXCALCULATE_TRACE(TEMP2,3,INE)
      END SUBROUTINE ComplexIN0INE      
        
!==================== Sigma ====================      
      SUBROUTINE SIGMA(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: F,FV,RES,INVFV,FE
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(3,3) :: L,L0,G,GE,SEQ,
     &    SNEQ,KROND,INVL
        REAL(KIND=8),DIMENSION(3) :: GD,GED
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,IN0,
     &   INE,PPEQ,PNEQ,Det,PEQ
        INTEGER :: INFO
        REAL(KIND=8),DIMENSION(3,3) :: TEMP2,TEMP332
        REAL(KIND=8),DIMENSION(3) :: TEMP1
        REAL(KIND=8) :: TEMP0,TEMP112
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        PPEN=PROPS(10)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        ! CALCULATE Fe=F * Fv_inv
        CALL MATRIX_DIVIDE(F,FV,FE,3)
        ! CALCULATE LPA AND LPER
        LPA=1.0D0 + 2.0D0*Q
        LPER=1.0D0 - Q
        ! CALCULATE L AND L0
        L=(LPA - LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + LPER * KROND
        L0=(LPA - LPER) * MATMUL(RESHAPE(D0,[3,1]),
     &   TRANSPOSE(RESHAPE(D0,[3,1]))) + LPER * KROND
        ! CALCULATE IN AND INE
        INVL=(1.0D0/LPA - 1.0D0/LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + 1.0D0/LPER * KROND
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(FE)),INVL),FE)
        CALL CALCULATE_TRACE(TEMP2,3,INE)
        ! CALCULATE G AND GE
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(F)),INVL),F)
        CALL CALCULATE_TRACE(TEMP2,3,IN0)
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        ! CALCULATE SEQ
        CALL MDet(F,Det)
        PEQ=MUEQ - PPEN*MUEQ*(Det-1.0D0)
        TEMP0=1.0D0 / (IM - IN0 + 3.0D0)
        TEMP0=1.0D0 / (IM - IN0 + 3.0D0)
        TEMP0=MUEQ/(1.0D0-Q)*IM*TEMP0
        TEMP112=3.0D0*Q/2.0D0/(1.0D0+2.0D0*Q)
        TEMP1=MATMUL(G,D)
        CALL OUTER_PRODUCT(D,TEMP1,TEMP2)
        CALL OUTER_PRODUCT(TEMP1,D,TEMP332)
        SEQ=TEMP0*(G - TEMP112*(TEMP2+TEMP332)) - PEQ*KROND
        ! CALCULATE SNEQ
        TEMP0=MUNEQ/(1.0D0-Q)*IM/(IM-INE+3.0D0)
        TEMP1=MATMUL(GE,D)
        CALL OUTER_PRODUCT(D,TEMP1,TEMP2)
        CALL OUTER_PRODUCT(TEMP1,D,TEMP332)
        SNEQ=TEMP0*(GE - TEMP112*(TEMP2+TEMP332)) - MUNEQ*KROND
        RES=SEQ + SNEQ
      END SUBROUTINE SIGMA     
      
!==================== ComplexSigma ==================== 
      SUBROUTINE ComplexSIGMA(PROPS,F,FV,D,D0,RES)
      IMPLICIT NONE
        REAL(KIND=8),DIMENSION(12) :: PROPS
        REAL(KIND=8),DIMENSION(3,3) :: FV,invFV
        COMPLEX,DIMENSION(3,3) :: F,RES,FE,TEMP2,G,GE,SEQ,SNEQ
        REAL(KIND=8),DIMENSION(3) :: D,D0
        REAL(KIND=8),DIMENSION(3,3) :: L,L0,
     &    KROND,INVL
        REAL(KIND=8),DIMENSION(3) :: GD,GED
        REAL(KIND=8) :: Q,IM,MUEQ,MUNEQ,PPEN,LPA,LPER,
     &   PPEQ,PNEQ
        COMPLEX :: PEQ,DET,INE,IN0
        INTEGER :: INFO,ISTAT
        COMPLEX,DIMENSION(3,3) :: TEMP332
        COMPLEX,DIMENSION(3) :: TEMP1
        REAL(KIND=8) :: TEMP0,TEMP112,Det_FV
        Q=PROPS(1)
        IM=PROPS(2)
        MUEQ=PROPS(4)
        MUNEQ=PROPS(5)
        PPEN=PROPS(10)
        CALL CREATE_IDENTITY_MATRIX(3,KROND)
        ! CALCULATE Fe=F * Fv_inv
        CALL MATINV3D(FV,invFV,Det_FV,ISTAT)
        FE = MATMUL(F,invFV)
        ! CALCULATE LPA AND LPER
        LPA=1.0D0 + 2.0D0*Q
        LPER=1.0D0 - Q
        ! CALCULATE L AND L0
        L=(LPA - LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + LPER * KROND
        L0=(LPA - LPER) * MATMUL(RESHAPE(D0,[3,1]),
     &   TRANSPOSE(RESHAPE(D0,[3,1]))) + LPER * KROND
        ! CALCULATE IN AND INE
        L=(1.0D0/LPA - 1.0D0/LPER) * MATMUL(RESHAPE(D,[3,1]),
     &   TRANSPOSE(RESHAPE(D,[3,1]))) + 1.0D0/LPER * KROND
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(FE)),INVL),FE)
        CALL ComplexCALCULATE_TRACE(TEMP2,3,INE)
        ! CALCULATE G AND GE
        G=MATMUL(MATMUL(F,L0),TRANSPOSE(F))
        TEMP2=MATMUL(MATMUL(MATMUL(L0,TRANSPOSE(F)),INVL),F)
        CALL ComplexCALCULATE_TRACE(TEMP2,3,IN0)
        GE=MATMUL(MATMUL(FE,L0),TRANSPOSE(FE))
        ! CALCULATE SEQ
        CALL ComplexMDet(F,Det)
        PEQ=MUEQ - PPEN*MUEQ*(Det-1.0D0)
        TEMP0=1.0D0 / (IM - IN0 + 3.0D0)
        TEMP0=1.0D0 / (IM - IN0 + 3.0D0)
        TEMP0=MUEQ/(1.0D0-Q)*IM*TEMP0
        TEMP112=3.0D0*Q/2.0D0/(1.0D0+2.0D0*Q)
        TEMP1=MATMUL(G,D)
        CALL ComplexROUTER_PRODUCT(D,TEMP1,TEMP2)
        CALL ComplexLOUTER_PRODUCT(TEMP1,D,TEMP332)
        SEQ=TEMP0*(G - TEMP112*(TEMP2+TEMP332)) - PEQ*KROND
        ! CALCULATE SNEQ
        TEMP0=MUNEQ/(1.0D0-Q)*IM/(IM-INE+3.0D0)
        TEMP1=MATMUL(GE,D)
        CALL ComplexROUTER_PRODUCT(D,TEMP1,TEMP2)
        CALL ComplexLOUTER_PRODUCT(TEMP1,D,TEMP332)
        SNEQ=TEMP0*(GE - TEMP112*(TEMP2+TEMP332)) - MUNEQ*KROND
        RES=SEQ + SNEQ
      END SUBROUTINE ComplexSIGMA          
      
!==================== Outer Product ====================     
      SUBROUTINE OUTER_PRODUCT(A,B,RES)
      IMPLICIT NONE
        INTEGER,PARAMETER :: n=3
        REAL(KIND=8),DIMENSION(n) :: A,B
        REAL(KIND=8),DIMENSION(n,n) :: RES
        INTEGER :: INDI,INDJ
        RES=0.0D0
        DO INDI=1,n
        DO INDJ=1,n
            RES(INDI,INDJ)=A(INDI) * B(INDJ)
        ENDDO
        ENDDO 
      END SUBROUTINE OUTER_PRODUCT
      
!==================== Outer Product ====================     
      SUBROUTINE ComplexLOUTER_PRODUCT(A,B,RES)
      IMPLICIT NONE
        INTEGER,PARAMETER :: n=3
        Complex,DIMENSION(n) :: A
        REAL(KIND=8),DIMENSION(n) :: B
        COMPLEX,DIMENSION(n,n) :: RES
        INTEGER :: INDI,INDJ
        RES=0.0D0
        DO INDI=1,n
        DO INDJ=1,n
            RES(INDI,INDJ)=A(INDI) * B(INDJ)
        ENDDO
        ENDDO 
      END SUBROUTINE ComplexLOUTER_PRODUCT  
      
!==================== Outer Product ====================     
      SUBROUTINE ComplexROUTER_PRODUCT(A,B,RES)
      IMPLICIT NONE
        INTEGER,PARAMETER :: n=3
        Complex,DIMENSION(n) :: B
        REAL(KIND=8),DIMENSION(n) :: A
        COMPLEX,DIMENSION(n,n) :: RES
        INTEGER :: INDI,INDJ
        RES=0.0D0
        DO INDI=1,n
        DO INDJ=1,n
            RES(INDI,INDJ)=A(INDI) * B(INDJ)
        ENDDO
        ENDDO 
      END SUBROUTINE ComplexROUTER_PRODUCT        
      
!==================== Inner Product ====================    
      SUBROUTINE COMPLEXINNER_PRODUCT(V1,V2,RESULT)
        IMPLICIT NONE
        COMPLEX,DIMENSION(3),INTENT(IN) :: V1,V2
        COMPLEX,INTENT(OUT) :: RESULT
        INTEGER :: INDI
        RESULT=0.0D0
        DO INDI=1,3
        RESULT=RESULT + V1(INDI) * V2(INDI)
        ENDDO
      END SUBROUTINE COMPLEXINNER_PRODUCT   
      
!==================== Inner Product ====================    
      SUBROUTINE INNER_PRODUCT(V1,V2,RESULT)
        IMPLICIT NONE
        REAL(KIND=8),DIMENSION(3),INTENT(IN) :: V1,V2
        REAL(KIND=8),INTENT(OUT) :: RESULT
        INTEGER :: INDI
        RESULT=0.0D0
        DO INDI=1,3
        RESULT=RESULT + V1(INDI) * V2(INDI)
        ENDDO
      END SUBROUTINE INNER_PRODUCT         
      
!==================== Matrix Inverse ==================== 
      SUBROUTINE MATRIX_INVERSE(A,A_INV,N,INFO)
      IMPLICIT NONE
!     Use LAPACK to compute the matrix inverse
      INTEGER :: N,INFO,IPIV(N),I,J
      DOUBLE PRECISION A(N,N),A_INV(N,N)
      DOUBLE PRECISION WORK(N)
!     Copy A to A_INV as DGETRF and DGETRI overWRITE their input
      DO I=1,N
         DO J=1,N
            A_INV(I,J)=A(I,J)
         ENDDO
      ENDDO
!     Perform LU decomposition of A_INV
      CALL DGETRF(N,N,A_INV,N,IPIV,INFO)
      IF (INFO .NE. 0) THEN
         PRINT *,'LU decomposition failed. Info=',INFO
         RETURN
      END IF
C     Compute the inverse of A from its LU decomposition
      CALL DGETRI(N,A_INV,N,IPIV,WORK,N,INFO)
      IF (INFO .NE. 0) THEN
         PRINT *,'Matrix inversion failed. Info=',INFO
      END IF
      END SUBROUTINE MATRIX_INVERSE

    
!==================== Matrix Divide ====================
       SUBROUTINE MATRIX_DIVIDE(A,B,C,N)
C     This SUBROUTINE divides matrix A by matrix B (A * B^-1) using a 
C     pre-defined MATRIX_INVERSE SUBROUTINE and efficient matrix multiplication.
      IMPLICIT NONE
      INTEGER N,INFO,I,J,K
      DOUBLE PRECISION A(N,N),B(N,N),C(N,N)
      DOUBLE PRECISION B_INV(N,N)
      DOUBLE PRECISION ALPHA,BETA,SUM
      
C     Compute the inverse of B into B_inv
      CALL MATRIX_INVERSE(B,B_INV,N,INFO)
      IF (INFO /= 0) THEN
          PRINT *,'Matrix inversion failed. Info=',INFO
          RETURN
      END IF
      
C     Initialize ALPHA and BETA for DGEMM usage
      ALPHA=1.0D0
      BETA=0.0D0
      
C     Efficiently multiply A and B_inv to get C using DGEMM (IF available)
C     Otherwise,use the manual loop for multiplication
C     Uncomment the following line IF DGEMM is available
C     CALL DGEMM('N','N',N,N,N,ALPHA,A,N,B_INV,N,BETA,C,N)
      
C     Manual loop for matrix multiplication IF DGEMM is not used
      DO I=1,N
          DO J=1,N
              SUM=0.0D0
              DO K=1,N
                  SUM=SUM + A(I,K) * B_INV(K,J)
              ENDDO
              C(I,J)=SUM
          ENDDO
      ENDDO
      END SUBROUTINE MATRIX_DIVIDE

    
!==================== Matrix Trace ====================
      SUBROUTINE CALCULATE_TRACE(MATRIX,N,TRACE_RESULT)
      IMPLICIT NONE
        INTEGER,INTENT(IN) :: N
        REAL(KIND=8),DIMENSION(N,N),INTENT(IN) :: MATRIX
        REAL(KIND=8),INTENT(OUT) :: TRACE_RESULT
        INTEGER :: I  
        TRACE_RESULT=0.0D0
        DO I=1,N
        TRACE_RESULT=TRACE_RESULT + MATRIX(I,I)
        ENDDO
      END SUBROUTINE CALCULATE_TRACE
      
      SUBROUTINE ComplexCALCULATE_TRACE(MATRIX,N,TRACE_RESULT)
      IMPLICIT NONE
        INTEGER,INTENT(IN) :: N
        COMPLEX,DIMENSION(N,N),INTENT(IN) :: MATRIX
        COMPLEX,INTENT(OUT) :: TRACE_RESULT
        INTEGER :: I  
        TRACE_RESULT=0.0D0
        DO I=1,N
        TRACE_RESULT=TRACE_RESULT + MATRIX(I,I)
        ENDDO
      END SUBROUTINE ComplexCALCULATE_TRACE      

!==================== Solve linear system ====================
      SUBROUTINE SOLVE_LINEAR_SYSTEM(KQ,KF,M,N,NRHS,INFO)
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: N,NRHS
          REAL(KIND=8),DIMENSION(N,N),INTENT(IN) :: KQ   ! The coefficient matrix
          REAL(KIND=8),DIMENSION(N,NRHS),INTENT(IN) :: KF ! The right-hand side matrix
          REAL(KIND=8),DIMENSION(N,NRHS),INTENT(OUT) :: M ! The solution matrix
          INTEGER,DIMENSION(N) :: IPIV               ! Pivot INDIces from LU decomposition
          INTEGER :: INFO
!         Copy KF to M,as DGESV overWRITEs the RHS matrix with the solution.
          M=KF
!         CALL DGESV to solve the system. KQ is overwritten by its LU decomposition.
          CALL DGESV(N,NRHS,KQ,N,IPIV,M,N,INFO)
      END SUBROUTINE SOLVE_LINEAR_SYSTEM


!==================== Create Identity Matrix ====================  
      SUBROUTINE CREATE_IDENTITY_MATRIX(N,MATRIX)
      IMPLICIT NONE
        INTEGER,INTENT(IN) :: N
        REAL(KIND=8),DIMENSION(N,N),INTENT(OUT) :: MATRIX
        INTEGER :: I,J
        MATRIX=0.0D0
        DO I=1,N
        MATRIX(I,I)=1.0D0
        ENDDO
      END SUBROUTINE CREATE_IDENTITY_MATRIX

!==================== xint2D1pt ====================                        
      SUBROUTINE xint2D1pt(xi,w,nInt)
      !
      ! This SUBROUTINE will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nINT,2): xi,eta coordinates for the integration pts
      !  w(nINT):    corresponding integration weights
      !
      IMPLICIT NONE
      !
      integer nInt
      !
      REAL(KIND=8) xi(1,2), w(1)
      
      ! Initialize
      !
      w=0.D0
      xi=0.D0


      ! Gauss weights
      !
      w=4.D0
      

      ! Gauss pt location in master element
      !
      xi(1,1)=0.D0
      xi(1,2)=0.D0


      return
      END SUBROUTINE xint2D1pt      
      
!==================== xint2D4pt ====================  

      SUBROUTINE xint2D4pt(xi,w,nINT)
      !
      ! This SUBROUTINE will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nINT,2): xi,eta coordinates for the integration pts
      !  w(nINT):    corresponding integration weights
      !
      IMPLICIT NONE
      !
      integer :: nINT
      !
      REAL(KIND=8) :: xi(4,2), w(4)

      ! Initialize
      !
      w=0.D0
      xi=0.D0

      ! Gauss weights
      !
      w(1)=1.D0
      w(2)=1.D0
      w(3)=1.D0
      w(4)=1.D0
      
      ! Gauss pt locations in master element
      !
      xi(1,1)=-dsqrt(1.D0/3.D0)
      xi(1,2)=-dsqrt(1.D0/3.D0)
      xi(2,1)=dsqrt(1.D0/3.D0)
      xi(2,2)=-dsqrt(1.D0/3.D0)
      xi(3,1)=-dsqrt(1.D0/3.D0)
      xi(3,2)=dsqrt(1.D0/3.D0)
      xi(4,1)=dsqrt(1.D0/3.D0)
      xi(4,2)=dsqrt(1.D0/3.D0)


      return
      END SUBROUTINE xint2D4pt
      
!==================== calcShape2DLinear ====================   

      SUBROUTINE calcShape2DLinear(nINT,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i)=shape function of node i at the intpt.
      ! dshxi(i,j)=derivative wrt j direction of shape fn of node i
      !
      IMPLICIT NONE
      !
      integer :: intpt,nDIM,nINT
      !
      REAL(KIND=8) :: xi_int(nINT,2),sh(4),dshxi(4,2),xi,eta
      !
      REAL(KIND=8),PARAMETER :: zero=0.D0,one=1.D0,fourth=1.D0/4.D0
      
      ! Location in the master element
      !
      xi=xi_int(intpt,1)
      eta=xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1)=fourth*(one - xi)*(one - eta)
      sh(2)=fourth*(one + xi)*(one - eta)
      sh(3)=fourth*(one + xi)*(one + eta)
      sh(4)=fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1)=-fourth*(one - eta)
      dshxi(1,2)=-fourth*(one - xi)
      dshxi(2,1)=fourth*(one - eta)
      dshxi(2,2)=-fourth*(one + xi)
      dshxi(3,1)=fourth*(one + eta)
      dshxi(3,2)=fourth*(one + xi)
      dshxi(4,1)=-fourth*(one + eta)
      dshxi(4,2)=fourth*(one - xi)

      return
      END SUBROUTINE calcShape2DLinear   
      
!==================== matInv2D ==================== 

      SUBROUTINE matInv2D(A,A_inv,Det_A,istat)
      !
      ! Returns A_inv, the inverse, and Det_A, the Determinant
      ! Note that the Det is of the original matrix, not the
      ! inverse
      !
      IMPLICIT NONE
      !
      integer istat
      !
      REAL(KIND=8) A(2,2),A_inv(2,2),Det_A,Det_A_inv
      istat=1
      
      Det_A=A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      IF (Det_A .le. 0.D0) then
        print *, 'WARNING: SUBROUTINE matInv2D:'
        print *, 'WARNING: Det of mat=',Det_A
        istat=0
        return
      end IF
            
      Det_A_inv=1.D0/Det_A
          
      A_inv(1,1)= Det_A_inv*A(2,2)
      A_inv(1,2)=-Det_A_inv*A(1,2)
      A_inv(2,1)=-Det_A_inv*A(2,1)
      A_inv(2,2)= Det_A_inv*A(1,1)


      return    
      END SUBROUTINE matInv2D      
      
!==================== MAPShape2D ====================  

      SUBROUTINE MAPShape2D(nNODE,dshxi,coords,dsh,DetMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      IMPLICIT NONE
      !
      integer :: i,j,k,nNODE,ieror,stat
      !
      REAL(KIND=8) :: dshxi(nNODE,2),dsh(nNODE,2),coords(2,nNODE),
     &  mapJ(2,2),mapJ_inv(2,2),DetmapJ
      !
      REAL(KIND=8),PARAMETER :: zero=0.D0,one=1.D0,two=2.D0,half=0.5D0,
     &     fourth=0.25D0,eighth=1.D0/8.D0
      !
      ! Calculate the mapping Jacobian matrix:
      !
      mapJ=zero
      DO i=1,2
        DO j=1,2
          DO k=1,nNODE
              mapJ(i,j)=mapJ(i,j) + dshxi(k,i)*coords(j,k)
          ENDDO
        ENDDO
      ENDDO

      ! Calculate the inverse and the Determinant of Jacobian
      !
      CALL matInv2D(mapJ,mapJ_inv,DetMapJ,stat)
      IF(stat.EQ.0) then
         print *, 'Problem: DetF.lt.zero in MAPShape2D'
         CALL XIT
      ENDIF


      ! Calculate first derivatives wrt x, y, z
      !
      dsh=transpose(matmul(mapJ_inv,transpose(dshxi)))
      
      return
      END SUBROUTINE MAPShape2D   
      
      
!==================== Det OF 3D MATRIX ====================
      SUBROUTINE MDet(A,Det)
            ! This SUBROUTINE calculates the Determinant
            ! of a 3 by 3 matrix [A]
            IMPLICIT NONE

            REAL(KIND=8) :: A(3,3),Det

            Det=A(1,1) * A(2,2) * A(3,3) + 
     &             A(1,2) * A(2,3) * A(3,1) + 
     &             A(1,3) * A(2,1) * A(3,2) - 
     &             A(3,1) * A(2,2) * A(1,3) - 
     &             A(3,2) * A(2,3) * A(1,1) - 
     &             A(3,3) * A(2,1) * A(1,2)
      RETURN      
      END SUBROUTINE MDet      

!==================== Det OF 3D MATRIX ====================
      SUBROUTINE ComplexMDet(A,Det)
            ! This SUBROUTINE calculates the Determinant
            ! of a 3 by 3 matrix [A]
            IMPLICIT NONE

            Complex :: A(3,3),Det

            Det=A(1,1) * A(2,2) * A(3,3) + 
     &             A(1,2) * A(2,3) * A(3,1) + 
     &             A(1,3) * A(2,1) * A(3,2) - 
     &             A(3,1) * A(2,2) * A(1,3) - 
     &             A(3,2) * A(2,3) * A(1,1) - 
     &             A(3,3) * A(2,1) * A(1,2)
      RETURN      
      END SUBROUTINE ComplexMDet                  
            
            
!==================== calcShape3DLinear ====================     
      subroutine calcShape3DLinear(nINT,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer :: intpt,nDim,nINT,i,j

      real(KIND=8) :: xi_int(nINT,3),sh(8),dshxi(8,3)
      real(KIND=8) :: d2shxi(8,3,3),xi,eta,zeta

      real(KIND=8),PARAMETER :: zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     &     fourth=0.25d0,eighth=1.d0/8.d0
  

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      END SUBROUTINE calcShape3DLinear           
      
!==================== mapShape3D ====================          
      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer :: i,j,k,nNode,ieror,stat

      real(kind=8) :: dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real(kind=8) :: mapJ(3,3),mapJ_inv(3,3),detmapJ

      real(kind=8),PARAMETER :: zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     &    fourth=0.25d0,eighth=1.d0/8.d0
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape3D'
         CALL XIT
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      END SUBROUTINE mapShape3D     
      
!==================== xint3D8pt ====================           
      subroutine xint3D8pt(xi,w,nInt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nINT,3): xi,eta,zeta coordinates for the integration pts
      !  w(nINT):    corresponding integration weights
      
      implicit none

      integer nInt

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      END SUBROUTINE xint3D8pt   
      
!==================== xint3D1pt ====================        
      subroutine xint3D1pt(xi,w,nINT)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nINT,3): xi,eta,zeta coordinates for the integration pts
      !  w(nINT):    corresponding integration weights
      
      implicit none

      integer nINT

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      END SUBROUTINE xint3D1pt         
            
            