! USER ELEMENT SUBROUTINE FOR QUADRATIC PLANE STRAIN COHESIVE ELEMENTS
! INCLUDES VISCOUS REGULARIZATION (GAO AND BOWER, 2004) 
! CODED BY TENG TONG (ASSISTANT PROF.SOUTHEAST UNIVERSITY)   

      
      
      
      
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPRO,PERIOD)
C-----------------------------------------------------------------------------------------C
      INCLUDE 'ABA_PARAM.INC' !IMPLICIT REAL(A-H O-Z)TM0
C-----------------------------------------------------------------------------------------C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),SVARS(*),
     1 ENERGY(*),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)     
C-----------------------------------------------------------------------------------------C
      IF (JTYPE .EQ. 1) THEN      ! CALL UEL FOR NONLOCAL MODEL
          CALL UEL1(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPRO,PERIOD)
      ELSE                        ! CALL UEL FOR COHESIVE MODEL
          CALL UEL6(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPRO,PERIOD)
      ENDIF
C-----------------------------------------------------------------------------------------C
      RETURN
      END   
C-----------------------------------------------------------------------------------------C
      
    
      
      
      
      
      
      
      
      
      
      
C-----------------------------------------------------------------------------------------C
      SUBROUTINE UEL6(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPRO,PERIOD)
C-----------------------------------------------------------------------------------------C
      INCLUDE 'ABA_PARAM.INC' !IMPLICIT REAL(A-H O-Z)
C-----------------------------------------------------------------------------------------C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),SVARS(*),
     1 ENERGY(*),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C-----------------------------------------------------------------------------------------C
       PARAMETER (ZERO = 0.D0, HALF=0.5D0, ONE= 1.0D0, TWO=2.0D0, 
     1     THREE= 3.0D0, TOL=-1E-5,NINTP = 3) 
       PARAMETER(N_ELEM=10000,NSDV=120,NGP=9)
C-----------------------------------------------------------------------------------------C
      DIMENSION SC(NDOFEL,NDOFEL),FC(NDOFEL,NRHS),T(MCRD,NRHS),
     1 T_D(MCRD,MCRD),R(MCRD,MCRD),BC(MCRD,NDOFEL),BCT(NDOFEL,MCRD),
     2 DEL(MCRD),TMP(NDOFEL,MCRD),RB(NDOFEL,NDOFEL),DELD(2)

C	GAUSS INTEGRATION VARIABLES (3 INTEG POINT)
      DIMENSION GAUSS3(3), WEIGHT3(3)

C	ARRAYS FOR QUADRATIC LINE ELEMENT
      DIMENSION DNDXI(3), DU_CONT(MCRD), DU_LOC(MCRD)
      DIMENSION DV_CONT(MCRD),DV_LOC(MCRD)
      DIMENSION H(MCRD,6), C_COOR(MCRD,NNODE), PSI(6,NDOFEL) 
      DIMENSION B(MCRD, NDOFEL), BT(NDOFEL, MCRD)
      DIMENSION A1(NDOFEL, MCRD), A2(NDOFEL, NDOFEL) 
      DIMENSION AV_COOR(MCRD, 3)
C
      COMMON/KUSER/USRVAR(N_ELEM,NSDV,NGP)
C-----------------------------------------------------------------------------------------C
C	INITILIZATION 
C     IMPORTANT!! FORTRAN DOES NOT PUT ZEROS IN THERE AUTOMATICALLY 
C
      CALL KASET2(AMATRX, NDOFEL, NDOFEL)
      CALL KASET1(RHS, MLVARX) 
C
      CALL KASET2(PSI, 6, NDOFEL)
	CALL KASET2(H, MCRD, 6)
      CALL KASET2(AV_COOR, MCRD, 3) 
      CALL KASET2(R, MCRD, MCRD)
      CALL KASET2(RB, NDOFEL, NDOFEL)
      CALL KASET2(T, MCRD, NRHS)
      CALL KASET2(T_D, MCRD, MCRD)
C     PAPARMETER INPUT	 
      WIDTH = 250.0d0	! WIDTH OF ELEMENTS (SAME AS SOLID SECTION WIDTH FOR SOLID ELEMENTS)
C-----------------------------------------------------------------------------------------C
C	RELATION MATRIX
      DO 10 K = 1, NDOFEL/2
          PSI(K, K) = -ONE
          PSI(K, K+NDOFEL/2) = ONE
10    END DO
C	COMPUTE NODAL COORDINATES IN DEFORMED STATE 
      DO 20 I=1,MCRD
          DO 30 J=1, NNODE
              NN=I+(J-1)*MCRD 
              C_COOR(I,J) = COORDS(I,J)  + U(NN)
30	    END DO
20    END DO
C	REFERENCE COORDINATE SYSTEM (MIDPOINT AVERAGES) 
      DO 31 I=1, MCRD
          DO 32 J=1, NNODE/2 
              AV_COOR(I,J)=ONE/TWO*(C_COOR(I,J)+C_COOR(I,J+NNODE/2))
32	    END DO
31    END DO
C-----------------------------------------------------------------------------------------C
C	GAUSSIAN POINT
      GAUSS3(1) = -SQRT(0.6)
      GAUSS3(2) = ZERO 
      GAUSS3(3) = SQRT(0.6)

      WEIGHT3(1) = 0.55555555555555
      WEIGHT3(2) = 0.88888888888888
      WEIGHT3(3) = 0.55555555555555
C-----------------------------------------------------------------------------------------C
C     TRANSFORMATION MATRIX
      CALL KCOORDTRANS (R,COORDS,U,NDOFEL,NNODE,MCRD)
C      
      DO I=1,MCRD
          DO J=1,MCRD
              DO K=0,NNODE-1
                  RB(I+2*K,J+2*K) = R(I,J)
              ENDDO
          ENDDO
      ENDDO
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
      DO IINTP=1,NINTP                                    !BEGIN LOOP
C-----------------------------------------------------------------------------------------C
C        
      POINT = GAUSS3(IINTP) 
      WEIGHT = WEIGHT3(IINTP)
      
C	SHAPE FUNCTION VALUE
      H1 = ONE/TWO*(-POINT + POINT**TWO) 
      H2 = ONE/TWO*( POINT + POINT**TWO) 
      H3 = ONE - POINT**TWO
      
C	DERIVATIVE OF SHAPE FUNCTION VALUE (3X1 MATRIX) 
      DNDXI(1) = -ONE/TWO + POINT
      DNDXI(2) = ONE/TWO + POINT 
      DNDXI(3) = -TWO*POINT

C	H MATRIX
      H(1,1) = H1
      H(2,2) = H1
      H(1,3) = H2
      H(2,4) = H2
      H(1,5) = H3
      H(2,6) = H3
      
C	B MATRIX: FROM NODAL DISPLACEMENT (U) TO DISPLCAMENE JUMP (DEALT_U)
      CALL KASET2(B, MCRD, NDOFEL) 
      DO 110 L=1, MCRD
        DO 120 J=1, NDOFEL
          DO 130 K=1, NDOFEL/2
              B(L,J) = B(L,J) + H(L,K)*PSI(K,J)
130	    END DO
120     END DO
110   END DO
      
C	TRANSPOSED B MATRIX 
      DO 140 L=1, MCRD
        DO 150 J=1, NDOFEL 
            BT(J,L) = B(L,J)
150	  END DO
140   END DO
C-----------------------------------------------------------------------------------------C
C 	CALCULATE GLOBAL DISPLACEMENT AT INTEGRATION POINT 
C     FROM CONTINUOUS DISPLACEMENT
	CALL KASET1(DU_CONT, MCRD)
      CALL KASET1(DV_CONT, MCRD)
	DO 160 L=1, MCRD
        DO 170 J=1, NDOFEL
            DU_CONT(L) = DU_CONT(L) + B(L,J)*U(J)
            DV_CONT(L) = DV_CONT(L) + B(L,J)*V(J)
170	  END DO
160   END DO
C-----------------------------------------------------------------------------------------C
C     LOCAL COORDINATE SYSTEM		
C     (USE AVERAGE OF DEFORMED X-POSITIONS OF TOP AND BOTTOM)
	X_XI = ZERO 
      Y_XI = ZERO
	DO 180 L=1,3
        X_XI = X_XI + DNDXI(L)*AV_COOR(1,L)
        Y_XI = Y_XI + DNDXI(L)*AV_COOR(2,L)
180	END DO

C	JACOBIAN (VECTOR LENGTH IN XI-DIRECTION) 
      DETJ = SQRT(X_XI**TWO + Y_XI**TWO)

C	RELATIVE DISPLACEMENT IN LOCAL COORDINATE SYSTEM
      CALL KASET1(DU_LOC, MCRD) 
      CALL KASET1(DV_LOC, MCRD)
      DO 181 I=1, MCRD
        DO 182 J=1, MCRD
            DU_LOC(I) = DU_LOC(I) + R(I,J)*DU_CONT(J)
            DV_LOC(I) = DV_LOC(I) + R(I,J)*DV_CONT(J)
182	  END DO
181   END DO
C      
      DEL  = DU_LOC
      DELD = DV_LOC
C-----------------------------------------------------------------------------------------C
C     DETERMINE LOAD CASE
      CALL KSEPLAW(PROPS,DEL,DELD,T,T_D,DTIME,damage,coeff,JELEM)
C     MODE II FRACTURE
      


C-----------------------------------------------------------------------------------------C

C-----------------------------------------------------------------------------------------C
C     ASSEMBLE AMATRX AND RHS
      BC  = MATMUL(B,RB)
      BCT = TRANSPOSE(BC)
      TMP = MATMUL(BCT,T_D)
      SC  = MATMUL(TMP,BC)
      FC  = MATMUL(BCT,T)
C       
      THICK = WIDTH*WEIGHT*DETJ
C
      CALL KMATRIXSCALAR(AMATRX,SC,THICK,NDOFEL,NDOFEL)
      CALL KMATRIXSCALAR(RHS,-FC,THICK,NDOFEL,NRHS)
C-----------------------------------------------------------------------------------------C
      
      USRVAR(JELEM,1,IINTP) = T(1,1)
      USRVAR(JELEM,1,IINTP+2) = T(1,1)
        
      USRVAR(JELEM,2,IINTP) = T(2,1)
      USRVAR(JELEM,2,IINTP+2) = T(2,1)
      
      USRVAR(JELEM,3,IINTP) = damage
      USRVAR(JELEM,3,IINTP+2) = damage
          
      USRVAR(JELEM,11,IINTP) = coeff
      USRVAR(JELEM,11,IINTP+2) =coeff


          
      ENDDO                                                       !END LOOP
C-----------------------------------------------------------------------------------------C
      RETURN
      END   

 
      
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C     SUBROUTINE KCOORDSTRANS TO OBTRAIN R, THE TWO DIMENSIONAL
C     TRANSFORMATION MATRIX AND EL_LENGTH, THE ELEMENT SIZE.
      SUBROUTINE KCOORDTRANS(R,COORDS,U,NDOFEL,NNODE,MCRD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(MCRD,MCRD),COORDS(MCRD,NNODE),U(NDOFEL)
      DIMENSION CO_DE(MCRD,NNODE),CO_DE_M(2,2)
C-----------------------------------------------------------------------------------------C
      DO I=1,MCRD
          DO J=1,NNODE
              CO_DE(I,J)=COORDS(I,J)+U(2*(J-1)+I)
          END DO
      END DO
C-----------------------------------------------------------------------------------------C
! CALCULATE OF THE DIRECTIONAL COSINE & THE TRANSFORMATION MATRIX
      D_X=CO_DE(1,3)-CO_DE(1,1)
      D_Y=CO_DE(2,3)-CO_DE(2,1)
      EL_LENGTH=(D_X**2+D_Y**2)**0.5D0
      COS_A=D_X/EL_LENGTH
      SIN_A=D_Y/EL_LENGTH
      R(1,1)=COS_A
      R(1,2)=SIN_A
      R(2,1)=-SIN_A
      R(2,2)=COS_A   
      RETURN
      END
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C     SUBROUTINE KMATRIX_PLUSSCALAR TO MULTIPLY A MATRIX
C     WITH A SCALAR NUMBER
      SUBROUTINE KMATRIXSCALAR(A,B,C,N,M)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(N,M),B(N,M)
      DO I=1,N
          DO J=1,M
              A(I,J)=A(I,J)+C*B(I,J)
          ENDDO
      ENDDO
      
      RETURN
      END
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
      SUBROUTINE KSEPLAW(PROPS,DEL,DELD,T,T_D,DTIME,damage,coeff,JELEM)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION PROPS(9),DEL(2),DELD(2),T(2,1),T_D(2,2)
C-----------------------------------------------------------------------------------------C      
      coeff = 1.0
        if (JELEM	.eq.	1833	)		coeff = 	0.771325974
        if (JELEM	.eq.	1834	)		coeff = 	0.771795337
        if (JELEM	.eq.	1835	)		coeff = 	0.930922658
        if (JELEM	.eq.	1836	)		coeff = 	0.689790467
        if (JELEM	.eq.	1837	)		coeff = 	0.911356732
        if (JELEM	.eq.	1838	)		coeff = 	0.920667932
        if (JELEM	.eq.	1839	)		coeff = 	0.975942318
        if (JELEM	.eq.	1840	)		coeff = 	1.005163124
        if (JELEM	.eq.	1841	)		coeff = 	0.785018643
        if (JELEM	.eq.	1842	)		coeff = 	0.848139616
        if (JELEM	.eq.	1843	)		coeff = 	0.856830769
        if (JELEM	.eq.	1844	)		coeff = 	1.174501132
        if (JELEM	.eq.	1845	)		coeff = 	1.069071236
        if (JELEM	.eq.	1846	)		coeff = 	1.071038429
        if (JELEM	.eq.	1847	)		coeff = 	0.997470041
        if (JELEM	.eq.	1848	)		coeff = 	0.79937783
        if (JELEM	.eq.	1849	)		coeff = 	0.777207256
        if (JELEM	.eq.	1850	)		coeff = 	1.057826993
        if (JELEM	.eq.	1851	)		coeff = 	1.042406501
        if (JELEM	.eq.	1852	)		coeff = 	0.660880621
        if (JELEM	.eq.	1853	)		coeff = 	0.879991899
        if (JELEM	.eq.	1854	)		coeff = 	0.821825677
        if (JELEM	.eq.	1855	)		coeff = 	0.91074392
        if (JELEM	.eq.	1856	)		coeff = 	0.958116319
        if (JELEM	.eq.	1857	)		coeff = 	1.058254722
        if (JELEM	.eq.	1858	)		coeff = 	0.973538449
        if (JELEM	.eq.	1859	)		coeff = 	1.052525714
        if (JELEM	.eq.	1860	)		coeff = 	1.075875218
        if (JELEM	.eq.	1861	)		coeff = 	1.089910981
        if (JELEM	.eq.	1862	)		coeff = 	1.05607143
        if (JELEM	.eq.	1863	)		coeff = 	1.127678254
        if (JELEM	.eq.	1864	)		coeff = 	0.909064987
        if (JELEM	.eq.	1865	)		coeff = 	0.986876414
        if (JELEM	.eq.	1866	)		coeff = 	0.799670381
        if (JELEM	.eq.	1867	)		coeff = 	0.989283276
        if (JELEM	.eq.	1868	)		coeff = 	0.978899423
        if (JELEM	.eq.	1869	)		coeff = 	0.977828417
        if (JELEM	.eq.	1870	)		coeff = 	0.942273045
        if (JELEM	.eq.	1871	)		coeff = 	0.930688913
        if (JELEM	.eq.	1872	)		coeff = 	0.683499757
        if (JELEM	.eq.	1873	)		coeff = 	0.996483261
        if (JELEM	.eq.	1874	)		coeff = 	1.087159602
        if (JELEM	.eq.	1875	)		coeff = 	1.003937015
        if (JELEM	.eq.	1876	)		coeff = 	1.07951303
        if (JELEM	.eq.	1877	)		coeff = 	1.038275997
        if (JELEM	.eq.	1878	)		coeff = 	0.870928604
        if (JELEM	.eq.	1879	)		coeff = 	0.950223242
        if (JELEM	.eq.	1880	)		coeff = 	0.929504977
        if (JELEM	.eq.	1881	)		coeff = 	0.908017237
        if (JELEM	.eq.	1882	)		coeff = 	0.741416071
        if (JELEM	.eq.	1883	)		coeff = 	1.023632617
        if (JELEM	.eq.	1884	)		coeff = 	0.785758548
        if (JELEM	.eq.	1885	)		coeff = 	1.05021177
        if (JELEM	.eq.	1886	)		coeff = 	0.954037312
        if (JELEM	.eq.	1887	)		coeff = 	0.967469173
        if (JELEM	.eq.	1888	)		coeff = 	1.047453841
        if (JELEM	.eq.	1889	)		coeff = 	0.845383515
        if (JELEM	.eq.	1890	)		coeff = 	0.834676204
        if (JELEM	.eq.	1891	)		coeff = 	1.057045876
        if (JELEM	.eq.	1892	)		coeff = 	1.12524065
        if (JELEM	.eq.	1893	)		coeff = 	0.892642783
        if (JELEM	.eq.	1894	)		coeff = 	1.053693267
        if (JELEM	.eq.	1895	)		coeff = 	0.960192019
        if (JELEM	.eq.	1896	)		coeff = 	1.08721226
        if (JELEM	.eq.	1897	)		coeff = 	1.067887088
        if (JELEM	.eq.	1898	)		coeff = 	0.982911363
        if (JELEM	.eq.	1899	)		coeff = 	1.078417674
        if (JELEM	.eq.	1900	)		coeff = 	0.869290018
        if (JELEM	.eq.	1901	)		coeff = 	0.884286095
        if (JELEM	.eq.	1902	)		coeff = 	1.010227969
        if (JELEM	.eq.	1903	)		coeff = 	0.940862707
        if (JELEM	.eq.	1904	)		coeff = 	1.041295545
        if (JELEM	.eq.	1905	)		coeff = 	0.931972772
        if (JELEM	.eq.	1906	)		coeff = 	1.087655219
        if (JELEM	.eq.	1907	)		coeff = 	1.000967448
        if (JELEM	.eq.	1908	)		coeff = 	1.112387295
        if (JELEM	.eq.	1909	)		coeff = 	0.990269538
        if (JELEM	.eq.	1910	)		coeff = 	0.957457267
        if (JELEM	.eq.	1911	)		coeff = 	1.068317139
        if (JELEM	.eq.	1912	)		coeff = 	0.980552727
        if (JELEM	.eq.	1913	)		coeff = 	1.038148502
        if (JELEM	.eq.	1914	)		coeff = 	0.960778787
        if (JELEM	.eq.	1915	)		coeff = 	0.976253176
        if (JELEM	.eq.	1916	)		coeff = 	0.980986538
        if (JELEM	.eq.	1917	)		coeff = 	0.879808883
        if (JELEM	.eq.	1918	)		coeff = 	0.722319835
        if (JELEM	.eq.	1919	)		coeff = 	0.883008638
        if (JELEM	.eq.	1920	)		coeff = 	0.788503392
        if (JELEM	.eq.	1921	)		coeff = 	0.959597946
        if (JELEM	.eq.	1922	)		coeff = 	0.984721649
        if (JELEM	.eq.	1923	)		coeff = 	0.894369568
        if (JELEM	.eq.	1924	)		coeff = 	1.119564661
        if (JELEM	.eq.	1925	)		coeff = 	0.92313854
        if (JELEM	.eq.	1926	)		coeff = 	1.048437714
        if (JELEM	.eq.	1927	)		coeff = 	0.949873064
        if (JELEM	.eq.	1928	)		coeff = 	0.902457291
        if (JELEM	.eq.	1929	)		coeff = 	0.850925154
        if (JELEM	.eq.	1930	)		coeff = 	0.975728269
        if (JELEM	.eq.	1931	)		coeff = 	1.117500237
        if (JELEM	.eq.	1932	)		coeff = 	1.123708937
        if (JELEM	.eq.	1933	)		coeff = 	0.931164945
        if (JELEM	.eq.	1934	)		coeff = 	0.989163466
        if (JELEM	.eq.	1935	)		coeff = 	0.938398054
        if (JELEM	.eq.	1936	)		coeff = 	1.093362599
        if (JELEM	.eq.	1937	)		coeff = 	0.950873295
        if (JELEM	.eq.	1938	)		coeff = 	0.983750358
        if (JELEM	.eq.	1939	)		coeff = 	0.797513276
        if (JELEM	.eq.	1940	)		coeff = 	1.099125328
        if (JELEM	.eq.	1941	)		coeff = 	1.061457242
        if (JELEM	.eq.	1942	)		coeff = 	1.0001273
        if (JELEM	.eq.	1943	)		coeff = 	1.175453703
        if (JELEM	.eq.	1944	)		coeff = 	0.811927122
        if (JELEM	.eq.	1945	)		coeff = 	0.975302624
        if (JELEM	.eq.	1946	)		coeff = 	0.942480873
        if (JELEM	.eq.	1947	)		coeff = 	1.057447687
        if (JELEM	.eq.	1948	)		coeff = 	0.935934373
        if (JELEM	.eq.	1949	)		coeff = 	1.041306308
        if (JELEM	.eq.	1950	)		coeff = 	0.891527664
        if (JELEM	.eq.	1951	)		coeff = 	1.0309402
        if (JELEM	.eq.	1952	)		coeff = 	1.050771748


C-----------------------------------------------------------------------------------------C      
      SIGMA_MAX = 1.0d0 
      TAU_MAX   = 1.55d0 * coeff
      DN = 0.01d0
      DT = 0.01d0
C      
      Q = 1.0 d0
      R = 0.0
      PHI_N = DEXP(1.D0)*SIGMA_MAX*DN
      
      DELT = ABS(DEL(1))
      DELN = DEL(2)
      
      if (del(1) .GE. 0) then
          sign_dt = 1.0d0
      else
          sign_dt = -1.0d0
      end if
      
      if ( del(2) .Ge. 0) then
          sign_dn = 1.0d0
      ELSE
          sign_dn = -1.0d0
      end if
C-----------------------------------------------------------------------------------------C      
      if (DELN .gt. 0.0d0) then
C         -----------------------------------------------------------------------------------------C      
          T22_TERM1 = SIGMA_MAX * DEXP(1.0d0)
          T22_TERM2 = DELN / DN
          T22_TERM3 = DEXP( -1.0D0 * DELN / DN )
          T22_TERM4 = DEXP( -1.0D0 * DELT * DELT / DT / DT )
	    T(2,1) = T22_TERM1 * T22_TERM2 * T22_TERM3 * T22_TERM4

          TD22_TERM1 = SIGMA_MAX * DEXP(1.0d0)
          TD22_TERM2 = 1.0D0/DN - DELN/DN/DN
          TD22_TERM3 = DEXP( -1.0D0 * DELN/ DN)
	    T_D(2,2) = TD22_TERM1 * TD22_TERM2 * TD22_TERM3 
          
          TD21_TERM1 = SIGMA_MAX * DEXP(1.0d0)
          TD21_TERM2 = DELN / DN
          TD21_TERM3 = DEXP( -1.0D0 * DELN / DN )
          TD21_TERM4 = DEXP( -1.0D0 * DELT * DELT / DT / DT )
          TD21_TERM5 = ( -2.0D0 * DELT * sign_dt / DT / DT )
	    T_D(2,1) = TD21_TERM1 * TD21_TERM2 * TD21_TERM3 * 
     +               TD21_TERM4 * TD21_TERM5

          T11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          T11_TERM2 = DELT / DT
          T11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
          T11_TERM4 = DEXP( -1.0D0 * DELN / DN)
	    T(1,1) = T11_TERM1 * T11_TERM2 * T11_TERM3 * sign_dt * T11_TERM4
      
          TD11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          TD11_TERM2 = 1.0D0/DT - 2.0D0*DELT*DELT/DT/DT/DT
          TD11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
          TD11_TERM4 = DEXP( -1.0D0 * DELN / DN)
	    T_D(1,1) = TD11_TERM1 * TD11_TERM2 * TD11_TERM3 * TD11_TERM4
          
          TD12_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          TD12_TERM2 = DELT / DT
          TD12_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
          TD12_TERM4 = DEXP( -1.0D0 * DELN / DN)
          TD12_TERM5 =( -1.0D0  / DN)
	    T_D(1,2) = TD12_TERM1 * TD12_TERM2 * TD12_TERM3 * 
     +               TD12_TERM4 * TD12_TERM5
      else
C         -----------------------------------------------------------------------------------------C      
	    T(2,1) = 1.0D5 * DELN

	    T_D(2,2) = 1.0D5
          
          T11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          T11_TERM2 = DELT / DT
          T11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
	    T(1,1) = T11_TERM1 * T11_TERM2 * T11_TERM3 * sign_dt
      
          TD11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          TD11_TERM2 = 1.0D0/DT - 2.0D0*DELT*DELT/DT/DT/DT
          TD11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
	    T_D(1,1) = TD11_TERM1 * TD11_TERM2 * TD11_TERM3
          
          FRIC_COEFF = 1.6D0
C
          IF (DELT .GE. 1.0d-15) THEN
              damage = 1.0d0 - abs(T(1,1)) / (TD11_TERM1/DT * DELT)
          ELSE
              damage = 0.0d0
          ENDIF
          
          T11_FRIC =  -1.0d0 * FRIC_COEFF * T(2,1) * sign_dt
     +                *  damage   
          T_D(1,2) =  -1.0d0 * FRIC_COEFF * T_D(2,2) * sign_dt 
     +                *  damage   
                
          T(1,1) = T(1,1) + T11_FRIC
          
          T_D(1,2) = 0.0d0

      end if 
      
      
      
C 	T_D(1,1) = 10000.0d0 
      
      
C     CROSS TERMS
      
C

      ZETA=0.1D0

	T(1,1)   = T(1,1)   + ZETA*TAU_MAX*DELD(1)/DT
	T_D(1,1) = T_D(1,1) + ZETA*TAU_MAX/DT/DTIME     
      
	T(2,1)   = T(2,1)   + ZETA*SIGMA_MAX*DELD(2)/DN
	T_D(2,2) = T_D(2,2) + ZETA*SIGMA_MAX/DN/DTIME   
      
      

C-----------------------------------------------------------------------------------------C
      RETURN
      END
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
      SUBROUTINE KASET2(DMATRIX,IDIMX,IDIMY)

      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO = 0.0D0) 
      DIMENSION DMATRIX(IDIMX,IDIMY)
      
      DO I = 1,IDIMX 
          DO J = 1,IDIMY
              DMATRIX(I,J) = ZERO
          END DO 
      END DO

      RETURN 
      END
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
C-----------------------------------------------------------------------------------------C
      SUBROUTINE KASET1(DMATRIX, IDIMX) 
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO = 0.0D0) 
      DIMENSION DMATRIX(IDIMX)

      DO I=1, IDIMX
          DMATRIX(I) = ZERO 
      END DO
      
      RETURN 
      END
C-----------------------------------------------------------------------------------------C
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
C ======================================================================
C User Subroutine UEL and UMAT for Localizing Gradient Damage
C Type: 8 noded SERENDIPITY element for displacement and e_eqv
C ======================================================================
C Please cite the related paper as:
C Authors- S Sarkar, IV Singh, BK Mishra, AS Shedbale, LH Poh
C Title- A comparative study and ABAQUS implementation of conventional and 
C        localizing gradient enhanced damage models
C Journal- Finite Elements in Analysis and Design 160 (2019) 1-31.
C DOI-10.1016/j.finel.2019.04.001
C ======================================================================
C Material properties to be given through the input file (*.inp), are
C
C For U1 element:
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = Averaging lenght parameter squared (c)
C PROPS(4) = Threshold fracture strain (ki)
C PROPS(5) = Second parameter in exp damage law (beta)
C
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C          by 2 - (N+N_UMAT)/2 (to be changed for each model)
C NSDV - solution dependent variables for the user element
C          strain(e)(x3)---> SDV = 1,2,3
C          loc_e_eqv(x1)---> SDV = 4
C          gp_e_eqv(x1) ---> SDV = 5
C          k_gp(x1)     ---> SDV = 6
C          k0(x1)       ---> SDV = 7
C          damage(D)(x1)---> SDV = 8
C          ....total=8
C NSDV - overall solution dependent variables (NSDV + 2), where
C        the additional 2 variables are the: time and iteration number
C ======================================================================
      SUBROUTINE UEL1(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     &     PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
C     ==================================================================
      PARAMETER(N_ELEM=10000,NSDV=120,NGP=9)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     &     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     &     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     &     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     &     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     &     JPROPS(*)
C
       DIMENSION GP(2,NGP),GW(2,NGP),XI(2),AJACOB(2,2),dNdxi(NNODE,2),
     & dHdxi(4,2),dNdx(NNODE,2),dHdx(4,2),AJABOBINV(2,2),AN(1,NNODE),
     & AH(1,4),BU(3,2*NNODE),BE(2,4),DMAT(3,3)
C
       DIMENSION E(3,NGP),ALOC_E_EQV(1,NGP),AK_GP(1,NGP),AK0(1,NGP)
       DIMENSION GP_E_EQV(1,NGP)
					  DIMENSION DISP(2*NNODE,1), E_EQV(4,1)
							DIMENSION DEEQV_DE(3), dJ2_de(3)
							DIMENSION del_oper(3), dI1_de(3), H_mat(3,3), del_E(2,1)
C
       DIMENSION AK_UU(2*NNODE,2*NNODE),AK_EE(4,4)
       DIMENSION AK_UE(2*NNODE,4),AK_EU(4,2*NNODE),F_EE(4,1)
C 
       COMMON/KUSER/USRVAR(N_ELEM,NSDV,NGP)
C
C     ******************************************************************
C     Constructing element U1
C     ******************************************************************
C     ==================================================================
C     Nodal variables saved in local array
C     ==================================================================
       coeff   = 1.0d0
       
        if (JELEM	.eq.	1	)		coeff = 	1.030988329
        if (JELEM	.eq.	2	)		coeff = 	0.809780715
        if (JELEM	.eq.	3	)		coeff = 	1.11617214
        if (JELEM	.eq.	4	)		coeff = 	1.009911783
        if (JELEM	.eq.	5	)		coeff = 	0.89099542
        if (JELEM	.eq.	6	)		coeff = 	0.964373205
        if (JELEM	.eq.	7	)		coeff = 	0.950820611
        if (JELEM	.eq.	8	)		coeff = 	1.040910799
        if (JELEM	.eq.	9	)		coeff = 	1.07679615
        if (JELEM	.eq.	10	)		coeff = 	0.968108434
        if (JELEM	.eq.	11	)		coeff = 	0.81830096
        if (JELEM	.eq.	12	)		coeff = 	0.922028663
        if (JELEM	.eq.	13	)		coeff = 	0.978229725
        if (JELEM	.eq.	14	)		coeff = 	1.071059972
        if (JELEM	.eq.	15	)		coeff = 	1.06949389
        if (JELEM	.eq.	16	)		coeff = 	0.667492107
        if (JELEM	.eq.	17	)		coeff = 	0.700668759
        if (JELEM	.eq.	18	)		coeff = 	0.983479254
        if (JELEM	.eq.	19	)		coeff = 	0.880366475
        if (JELEM	.eq.	20	)		coeff = 	0.947365446
        if (JELEM	.eq.	21	)		coeff = 	0.880109804
        if (JELEM	.eq.	22	)		coeff = 	0.68975315
        if (JELEM	.eq.	23	)		coeff = 	1.122401466
        if (JELEM	.eq.	24	)		coeff = 	1.004599531
        if (JELEM	.eq.	25	)		coeff = 	0.970466342
        if (JELEM	.eq.	26	)		coeff = 	1.057153817
        if (JELEM	.eq.	27	)		coeff = 	1.003261
        if (JELEM	.eq.	28	)		coeff = 	0.807759269
        if (JELEM	.eq.	29	)		coeff = 	1.103251256
        if (JELEM	.eq.	30	)		coeff = 	1.090483017
        if (JELEM	.eq.	31	)		coeff = 	1.046859795
        if (JELEM	.eq.	32	)		coeff = 	0.930350514
        if (JELEM	.eq.	33	)		coeff = 	1.026322908
        if (JELEM	.eq.	34	)		coeff = 	0.883762217
        if (JELEM	.eq.	35	)		coeff = 	1.017124731
        if (JELEM	.eq.	36	)		coeff = 	0.8763946
        if (JELEM	.eq.	37	)		coeff = 	1.032344163
        if (JELEM	.eq.	38	)		coeff = 	0.912206557
        if (JELEM	.eq.	39	)		coeff = 	1.07313772
        if (JELEM	.eq.	40	)		coeff = 	0.753852719
        if (JELEM	.eq.	41	)		coeff = 	0.979223055
        if (JELEM	.eq.	42	)		coeff = 	1.002690717
        if (JELEM	.eq.	43	)		coeff = 	1.058946496
        if (JELEM	.eq.	44	)		coeff = 	0.964269718
        if (JELEM	.eq.	45	)		coeff = 	0.977160577
        if (JELEM	.eq.	46	)		coeff = 	1.105030637
        if (JELEM	.eq.	47	)		coeff = 	0.808210185
        if (JELEM	.eq.	48	)		coeff = 	0.791707033
        if (JELEM	.eq.	49	)		coeff = 	0.868816337
        if (JELEM	.eq.	50	)		coeff = 	0.971757659
        if (JELEM	.eq.	51	)		coeff = 	0.862559708
        if (JELEM	.eq.	52	)		coeff = 	0.817178861
        if (JELEM	.eq.	53	)		coeff = 	0.91178505
        if (JELEM	.eq.	54	)		coeff = 	1.102521172
        if (JELEM	.eq.	55	)		coeff = 	1.076175259
        if (JELEM	.eq.	56	)		coeff = 	0.984434213
        if (JELEM	.eq.	57	)		coeff = 	0.868811154
        if (JELEM	.eq.	58	)		coeff = 	0.938322566
        if (JELEM	.eq.	59	)		coeff = 	0.971704912
        if (JELEM	.eq.	60	)		coeff = 	1.12976904
        if (JELEM	.eq.	61	)		coeff = 	1.055385108
        if (JELEM	.eq.	62	)		coeff = 	0.767167684
        if (JELEM	.eq.	63	)		coeff = 	0.937835193
        if (JELEM	.eq.	64	)		coeff = 	0.962246571
        if (JELEM	.eq.	65	)		coeff = 	1.042845651
        if (JELEM	.eq.	66	)		coeff = 	0.739844658
        if (JELEM	.eq.	67	)		coeff = 	1.001528882
        if (JELEM	.eq.	68	)		coeff = 	0.975421591
        if (JELEM	.eq.	69	)		coeff = 	0.844890044
        if (JELEM	.eq.	70	)		coeff = 	0.865136761
        if (JELEM	.eq.	71	)		coeff = 	0.980190065
        if (JELEM	.eq.	72	)		coeff = 	0.927583136
        if (JELEM	.eq.	73	)		coeff = 	0.734196775
        if (JELEM	.eq.	74	)		coeff = 	0.914757057
        if (JELEM	.eq.	75	)		coeff = 	1.028340114
        if (JELEM	.eq.	76	)		coeff = 	0.893453963
        if (JELEM	.eq.	77	)		coeff = 	0.946640182
        if (JELEM	.eq.	78	)		coeff = 	0.697848089
        if (JELEM	.eq.	79	)		coeff = 	1.018489995
        if (JELEM	.eq.	80	)		coeff = 	0.96484027
        if (JELEM	.eq.	81	)		coeff = 	1.023294393
        if (JELEM	.eq.	82	)		coeff = 	1.100302465
        if (JELEM	.eq.	83	)		coeff = 	1.081774127
        if (JELEM	.eq.	84	)		coeff = 	1.099678985
        if (JELEM	.eq.	85	)		coeff = 	0.572996804
        if (JELEM	.eq.	86	)		coeff = 	0.935156413
        if (JELEM	.eq.	87	)		coeff = 	0.865924924
        if (JELEM	.eq.	88	)		coeff = 	0.826837524
        if (JELEM	.eq.	89	)		coeff = 	0.846925456
        if (JELEM	.eq.	90	)		coeff = 	0.761015856
        if (JELEM	.eq.	91	)		coeff = 	1.127774952
        if (JELEM	.eq.	92	)		coeff = 	1.070055982
        if (JELEM	.eq.	93	)		coeff = 	0.892579671
        if (JELEM	.eq.	94	)		coeff = 	0.819647233
        if (JELEM	.eq.	95	)		coeff = 	0.998309538
        if (JELEM	.eq.	96	)		coeff = 	1.043553681
        if (JELEM	.eq.	97	)		coeff = 	0.942247067
        if (JELEM	.eq.	98	)		coeff = 	1.138643318
        if (JELEM	.eq.	99	)		coeff = 	1.121042578
        if (JELEM	.eq.	100	)		coeff = 	0.765247914
        if (JELEM	.eq.	101	)		coeff = 	0.990849172
        if (JELEM	.eq.	102	)		coeff = 	0.867785686
        if (JELEM	.eq.	103	)		coeff = 	1.035411968
        if (JELEM	.eq.	104	)		coeff = 	0.93241313
        if (JELEM	.eq.	105	)		coeff = 	0.8040694
        if (JELEM	.eq.	106	)		coeff = 	0.925878589
        if (JELEM	.eq.	107	)		coeff = 	0.775709375
        if (JELEM	.eq.	108	)		coeff = 	0.578587295
        if (JELEM	.eq.	109	)		coeff = 	1.066094655
        if (JELEM	.eq.	110	)		coeff = 	1.139686425
        if (JELEM	.eq.	111	)		coeff = 	0.788518887
        if (JELEM	.eq.	112	)		coeff = 	1.066881224
        if (JELEM	.eq.	113	)		coeff = 	1.09188294
        if (JELEM	.eq.	114	)		coeff = 	0.896272393
        if (JELEM	.eq.	115	)		coeff = 	1.057126409
        if (JELEM	.eq.	116	)		coeff = 	1.062545612
        if (JELEM	.eq.	117	)		coeff = 	0.969838379
        if (JELEM	.eq.	118	)		coeff = 	0.744942498
        if (JELEM	.eq.	119	)		coeff = 	0.886467622
        if (JELEM	.eq.	120	)		coeff = 	0.911757076
        if (JELEM	.eq.	121	)		coeff = 	1.102543304
        if (JELEM	.eq.	122	)		coeff = 	1.017018305
        if (JELEM	.eq.	123	)		coeff = 	0.76112883
        if (JELEM	.eq.	124	)		coeff = 	0.829587428
        if (JELEM	.eq.	125	)		coeff = 	0.805914703
        if (JELEM	.eq.	126	)		coeff = 	0.938490973
        if (JELEM	.eq.	127	)		coeff = 	1.061511355
        if (JELEM	.eq.	128	)		coeff = 	1.056769708
        if (JELEM	.eq.	129	)		coeff = 	0.883852363
        if (JELEM	.eq.	130	)		coeff = 	0.749939289
        if (JELEM	.eq.	131	)		coeff = 	0.916596
        if (JELEM	.eq.	132	)		coeff = 	0.951633946
        if (JELEM	.eq.	133	)		coeff = 	0.935748626
        if (JELEM	.eq.	134	)		coeff = 	1.037845334
        if (JELEM	.eq.	135	)		coeff = 	0.850799583
        if (JELEM	.eq.	136	)		coeff = 	0.866143518
        if (JELEM	.eq.	137	)		coeff = 	0.938348734
        if (JELEM	.eq.	138	)		coeff = 	0.981482475
        if (JELEM	.eq.	139	)		coeff = 	0.949932377
        if (JELEM	.eq.	140	)		coeff = 	1.00406251
        if (JELEM	.eq.	141	)		coeff = 	0.578771916
        if (JELEM	.eq.	142	)		coeff = 	1.00079293
        if (JELEM	.eq.	143	)		coeff = 	0.713952818
        if (JELEM	.eq.	144	)		coeff = 	1.139756136
        if (JELEM	.eq.	145	)		coeff = 	0.87097986
        if (JELEM	.eq.	146	)		coeff = 	0.980418258
        if (JELEM	.eq.	147	)		coeff = 	0.925535129
        if (JELEM	.eq.	148	)		coeff = 	1.014802416
        if (JELEM	.eq.	149	)		coeff = 	0.778158379
        if (JELEM	.eq.	150	)		coeff = 	1.04131925
        if (JELEM	.eq.	151	)		coeff = 	0.998205593
        if (JELEM	.eq.	152	)		coeff = 	1.016446535
        if (JELEM	.eq.	153	)		coeff = 	1.125950155
        if (JELEM	.eq.	154	)		coeff = 	0.933451341
        if (JELEM	.eq.	155	)		coeff = 	0.842083862
        if (JELEM	.eq.	156	)		coeff = 	0.938591636
        if (JELEM	.eq.	157	)		coeff = 	1.001203496
        if (JELEM	.eq.	158	)		coeff = 	1.069988139
        if (JELEM	.eq.	159	)		coeff = 	1.079531333
        if (JELEM	.eq.	160	)		coeff = 	0.672646023
        if (JELEM	.eq.	161	)		coeff = 	0.952248045
        if (JELEM	.eq.	162	)		coeff = 	1.093634507
        if (JELEM	.eq.	163	)		coeff = 	1.031392934
        if (JELEM	.eq.	164	)		coeff = 	1.045587322
        if (JELEM	.eq.	165	)		coeff = 	1.014639205
        if (JELEM	.eq.	166	)		coeff = 	0.795276901
        if (JELEM	.eq.	167	)		coeff = 	0.958217127
        if (JELEM	.eq.	168	)		coeff = 	1.061637332
        if (JELEM	.eq.	169	)		coeff = 	1.068549061
        if (JELEM	.eq.	170	)		coeff = 	1.106375006
        if (JELEM	.eq.	171	)		coeff = 	0.831015123
        if (JELEM	.eq.	172	)		coeff = 	0.693886142
        if (JELEM	.eq.	173	)		coeff = 	0.770369885
        if (JELEM	.eq.	174	)		coeff = 	0.988721345
        if (JELEM	.eq.	175	)		coeff = 	0.797120875
        if (JELEM	.eq.	176	)		coeff = 	0.968305601
        if (JELEM	.eq.	177	)		coeff = 	0.972379751
        if (JELEM	.eq.	178	)		coeff = 	1.05964235
        if (JELEM	.eq.	179	)		coeff = 	1.095546383
        if (JELEM	.eq.	180	)		coeff = 	0.92468643
        if (JELEM	.eq.	181	)		coeff = 	1.015653117
        if (JELEM	.eq.	182	)		coeff = 	0.973288628
        if (JELEM	.eq.	183	)		coeff = 	1.04517474
        if (JELEM	.eq.	184	)		coeff = 	1.206568873
        if (JELEM	.eq.	185	)		coeff = 	0.861859701
        if (JELEM	.eq.	186	)		coeff = 	0.999954338
        if (JELEM	.eq.	187	)		coeff = 	1.154057283
        if (JELEM	.eq.	188	)		coeff = 	0.896539899
        if (JELEM	.eq.	189	)		coeff = 	1.042801751
        if (JELEM	.eq.	190	)		coeff = 	0.906070696
        if (JELEM	.eq.	191	)		coeff = 	1.034396438
        if (JELEM	.eq.	192	)		coeff = 	1.02845864
        if (JELEM	.eq.	193	)		coeff = 	1.073074801
        if (JELEM	.eq.	194	)		coeff = 	1.006459219
        if (JELEM	.eq.	195	)		coeff = 	0.881051227
        if (JELEM	.eq.	196	)		coeff = 	0.807166215
        if (JELEM	.eq.	197	)		coeff = 	0.898228533
        if (JELEM	.eq.	198	)		coeff = 	1.034256226
        if (JELEM	.eq.	199	)		coeff = 	0.97737207
        if (JELEM	.eq.	200	)		coeff = 	0.767859243
        if (JELEM	.eq.	201	)		coeff = 	0.930204596
        if (JELEM	.eq.	202	)		coeff = 	1.030745284
        if (JELEM	.eq.	203	)		coeff = 	0.925343054
        if (JELEM	.eq.	204	)		coeff = 	0.903649455
        if (JELEM	.eq.	205	)		coeff = 	0.990306143
        if (JELEM	.eq.	206	)		coeff = 	1.050878312
        if (JELEM	.eq.	207	)		coeff = 	1.016746527
        if (JELEM	.eq.	208	)		coeff = 	0.904569328
        if (JELEM	.eq.	209	)		coeff = 	0.925579795
        if (JELEM	.eq.	210	)		coeff = 	1.03325521
        if (JELEM	.eq.	211	)		coeff = 	1.121110773
        if (JELEM	.eq.	212	)		coeff = 	0.955171955
        if (JELEM	.eq.	213	)		coeff = 	1.056715695
        if (JELEM	.eq.	214	)		coeff = 	0.957876968
        if (JELEM	.eq.	215	)		coeff = 	0.820203449
        if (JELEM	.eq.	216	)		coeff = 	0.930843727
        if (JELEM	.eq.	217	)		coeff = 	1.026902467
        if (JELEM	.eq.	218	)		coeff = 	0.831080631
        if (JELEM	.eq.	219	)		coeff = 	1.063368752
        if (JELEM	.eq.	220	)		coeff = 	1.030412796
        if (JELEM	.eq.	221	)		coeff = 	0.912354889
        if (JELEM	.eq.	222	)		coeff = 	1.040502278
        if (JELEM	.eq.	223	)		coeff = 	0.847518178
        if (JELEM	.eq.	224	)		coeff = 	1.012002053
        if (JELEM	.eq.	225	)		coeff = 	0.967231854
        if (JELEM	.eq.	226	)		coeff = 	0.877129655
        if (JELEM	.eq.	227	)		coeff = 	0.887816551
        if (JELEM	.eq.	228	)		coeff = 	0.923769005
        if (JELEM	.eq.	229	)		coeff = 	1.033458989
        if (JELEM	.eq.	230	)		coeff = 	0.889139279
        if (JELEM	.eq.	231	)		coeff = 	1.007750631
        if (JELEM	.eq.	232	)		coeff = 	1.017956268
        if (JELEM	.eq.	233	)		coeff = 	1.075074002
        if (JELEM	.eq.	234	)		coeff = 	1.164643424
        if (JELEM	.eq.	235	)		coeff = 	0.774873689
        if (JELEM	.eq.	236	)		coeff = 	0.789683297
        if (JELEM	.eq.	237	)		coeff = 	0.99612798
        if (JELEM	.eq.	238	)		coeff = 	0.901539336
        if (JELEM	.eq.	239	)		coeff = 	0.911754183
        if (JELEM	.eq.	240	)		coeff = 	1.07485609
        if (JELEM	.eq.	241	)		coeff = 	1.048423624
        if (JELEM	.eq.	242	)		coeff = 	0.796078035
        if (JELEM	.eq.	243	)		coeff = 	0.805137737
        if (JELEM	.eq.	244	)		coeff = 	0.997338695
        if (JELEM	.eq.	245	)		coeff = 	1.033145232
        if (JELEM	.eq.	246	)		coeff = 	0.934208137
        if (JELEM	.eq.	247	)		coeff = 	0.927136083
        if (JELEM	.eq.	248	)		coeff = 	1.146613919
        if (JELEM	.eq.	249	)		coeff = 	0.850817927
        if (JELEM	.eq.	250	)		coeff = 	0.95899233
        if (JELEM	.eq.	251	)		coeff = 	0.790208195
        if (JELEM	.eq.	252	)		coeff = 	1.078634749
        if (JELEM	.eq.	253	)		coeff = 	0.748713772
        if (JELEM	.eq.	254	)		coeff = 	0.772574188
        if (JELEM	.eq.	255	)		coeff = 	1.089472553
        if (JELEM	.eq.	256	)		coeff = 	0.987169681
        if (JELEM	.eq.	257	)		coeff = 	1.012344162
        if (JELEM	.eq.	258	)		coeff = 	0.836093919
        if (JELEM	.eq.	259	)		coeff = 	0.746100546
        if (JELEM	.eq.	260	)		coeff = 	0.891124383
        if (JELEM	.eq.	261	)		coeff = 	0.907779377
        if (JELEM	.eq.	262	)		coeff = 	0.922032983
        if (JELEM	.eq.	263	)		coeff = 	1.014284401
        if (JELEM	.eq.	264	)		coeff = 	0.889987748
        if (JELEM	.eq.	265	)		coeff = 	1.063299994
        if (JELEM	.eq.	266	)		coeff = 	0.87867753
        if (JELEM	.eq.	267	)		coeff = 	1.090898309
        if (JELEM	.eq.	268	)		coeff = 	0.710712752
        if (JELEM	.eq.	269	)		coeff = 	0.932282411
        if (JELEM	.eq.	270	)		coeff = 	0.802650188
        if (JELEM	.eq.	271	)		coeff = 	1.050658698
        if (JELEM	.eq.	272	)		coeff = 	0.642641157
        if (JELEM	.eq.	273	)		coeff = 	1.035808217
        if (JELEM	.eq.	274	)		coeff = 	1.041289289
        if (JELEM	.eq.	275	)		coeff = 	0.990755594
        if (JELEM	.eq.	276	)		coeff = 	0.869605777
        if (JELEM	.eq.	277	)		coeff = 	0.840187443
        if (JELEM	.eq.	278	)		coeff = 	0.988726874
        if (JELEM	.eq.	279	)		coeff = 	0.990646716
        if (JELEM	.eq.	280	)		coeff = 	0.873666798
        if (JELEM	.eq.	281	)		coeff = 	1.03346384
        if (JELEM	.eq.	282	)		coeff = 	0.823924305
        if (JELEM	.eq.	283	)		coeff = 	0.908359554
        if (JELEM	.eq.	284	)		coeff = 	0.685361847
        if (JELEM	.eq.	285	)		coeff = 	1.053534946
        if (JELEM	.eq.	286	)		coeff = 	1.04355518
        if (JELEM	.eq.	287	)		coeff = 	1.015110228
        if (JELEM	.eq.	288	)		coeff = 	0.882033206
        if (JELEM	.eq.	289	)		coeff = 	1.000038827
        if (JELEM	.eq.	290	)		coeff = 	0.933670416
        if (JELEM	.eq.	291	)		coeff = 	0.995779275
        if (JELEM	.eq.	292	)		coeff = 	0.949716158
        if (JELEM	.eq.	293	)		coeff = 	0.99523548
        if (JELEM	.eq.	294	)		coeff = 	0.914098428
        if (JELEM	.eq.	295	)		coeff = 	1.154836186
        if (JELEM	.eq.	296	)		coeff = 	0.854220479
        if (JELEM	.eq.	297	)		coeff = 	1.004913303
        if (JELEM	.eq.	298	)		coeff = 	0.896339047
        if (JELEM	.eq.	299	)		coeff = 	0.728281059
        if (JELEM	.eq.	300	)		coeff = 	0.961944432
        if (JELEM	.eq.	301	)		coeff = 	0.913850993
        if (JELEM	.eq.	302	)		coeff = 	0.983324068
        if (JELEM	.eq.	303	)		coeff = 	0.834157015
        if (JELEM	.eq.	304	)		coeff = 	0.917838766
        if (JELEM	.eq.	305	)		coeff = 	0.962826392
        if (JELEM	.eq.	306	)		coeff = 	0.999826964
        if (JELEM	.eq.	307	)		coeff = 	0.907289754
        if (JELEM	.eq.	308	)		coeff = 	1.087355054
        if (JELEM	.eq.	309	)		coeff = 	1.117747211
        if (JELEM	.eq.	310	)		coeff = 	0.821524969
        if (JELEM	.eq.	311	)		coeff = 	0.918619513
        if (JELEM	.eq.	312	)		coeff = 	0.931821322
        if (JELEM	.eq.	313	)		coeff = 	0.985360098
        if (JELEM	.eq.	314	)		coeff = 	0.962979853
        if (JELEM	.eq.	315	)		coeff = 	1.060953951
        if (JELEM	.eq.	316	)		coeff = 	1.090273322
        if (JELEM	.eq.	317	)		coeff = 	0.907263327
        if (JELEM	.eq.	318	)		coeff = 	1.068910447
        if (JELEM	.eq.	319	)		coeff = 	0.955950082
        if (JELEM	.eq.	320	)		coeff = 	1.188734505
        if (JELEM	.eq.	321	)		coeff = 	0.984517707
        if (JELEM	.eq.	322	)		coeff = 	1.041003437
        if (JELEM	.eq.	323	)		coeff = 	0.915941803
        if (JELEM	.eq.	324	)		coeff = 	0.963148619
        if (JELEM	.eq.	325	)		coeff = 	0.935016245
        if (JELEM	.eq.	326	)		coeff = 	1.048682441
        if (JELEM	.eq.	327	)		coeff = 	1.070958893
        if (JELEM	.eq.	328	)		coeff = 	1.140959152
        if (JELEM	.eq.	329	)		coeff = 	0.757690798
        if (JELEM	.eq.	330	)		coeff = 	0.903660501
        if (JELEM	.eq.	331	)		coeff = 	0.948631838
        if (JELEM	.eq.	332	)		coeff = 	1.020031969
        if (JELEM	.eq.	333	)		coeff = 	1.121765514
        if (JELEM	.eq.	334	)		coeff = 	0.83156181
        if (JELEM	.eq.	335	)		coeff = 	0.936043793
        if (JELEM	.eq.	336	)		coeff = 	1.148368772
        if (JELEM	.eq.	337	)		coeff = 	0.765554907
        if (JELEM	.eq.	338	)		coeff = 	1.165536469
        if (JELEM	.eq.	339	)		coeff = 	0.922180655
        if (JELEM	.eq.	340	)		coeff = 	0.866055727
        if (JELEM	.eq.	341	)		coeff = 	0.956576164
        if (JELEM	.eq.	342	)		coeff = 	1.021213827
        if (JELEM	.eq.	343	)		coeff = 	0.926225063
        if (JELEM	.eq.	344	)		coeff = 	0.984406465
        if (JELEM	.eq.	345	)		coeff = 	0.950903899
        if (JELEM	.eq.	346	)		coeff = 	1.051768885
        if (JELEM	.eq.	347	)		coeff = 	0.885096563
        if (JELEM	.eq.	348	)		coeff = 	1.044606961
        if (JELEM	.eq.	349	)		coeff = 	1.086922253
        if (JELEM	.eq.	350	)		coeff = 	0.801485908
        if (JELEM	.eq.	351	)		coeff = 	0.849450351
        if (JELEM	.eq.	352	)		coeff = 	0.812040971
        if (JELEM	.eq.	353	)		coeff = 	1.092805512
        if (JELEM	.eq.	354	)		coeff = 	0.901975231
        if (JELEM	.eq.	355	)		coeff = 	1.058797717
        if (JELEM	.eq.	356	)		coeff = 	0.843657874
        if (JELEM	.eq.	357	)		coeff = 	0.891573977
        if (JELEM	.eq.	358	)		coeff = 	0.87184444
        if (JELEM	.eq.	359	)		coeff = 	0.87722698
        if (JELEM	.eq.	360	)		coeff = 	1.110916804
        if (JELEM	.eq.	361	)		coeff = 	0.932010207
        if (JELEM	.eq.	362	)		coeff = 	0.799765576
        if (JELEM	.eq.	363	)		coeff = 	0.959502486
        if (JELEM	.eq.	364	)		coeff = 	0.962292356
        if (JELEM	.eq.	365	)		coeff = 	1.031534795
        if (JELEM	.eq.	366	)		coeff = 	0.70137631
        if (JELEM	.eq.	367	)		coeff = 	0.940726296
        if (JELEM	.eq.	368	)		coeff = 	0.709060143
        if (JELEM	.eq.	369	)		coeff = 	1.122770878
        if (JELEM	.eq.	370	)		coeff = 	1.052492506
        if (JELEM	.eq.	371	)		coeff = 	0.958481054
        if (JELEM	.eq.	372	)		coeff = 	0.922213523
        if (JELEM	.eq.	373	)		coeff = 	1.007413808
        if (JELEM	.eq.	374	)		coeff = 	0.919781043
        if (JELEM	.eq.	375	)		coeff = 	1.003388322
        if (JELEM	.eq.	376	)		coeff = 	0.833594856
        if (JELEM	.eq.	377	)		coeff = 	0.951928935
        if (JELEM	.eq.	378	)		coeff = 	0.911501858
        if (JELEM	.eq.	379	)		coeff = 	1.031252145
        if (JELEM	.eq.	380	)		coeff = 	1.077484967
        if (JELEM	.eq.	381	)		coeff = 	1.123484087
        if (JELEM	.eq.	382	)		coeff = 	0.80816313
        if (JELEM	.eq.	383	)		coeff = 	0.881548715
        if (JELEM	.eq.	384	)		coeff = 	0.859706679
        if (JELEM	.eq.	385	)		coeff = 	0.980086111
        if (JELEM	.eq.	386	)		coeff = 	1.099754889
        if (JELEM	.eq.	387	)		coeff = 	1.064555199
        if (JELEM	.eq.	388	)		coeff = 	1.034025908
        if (JELEM	.eq.	389	)		coeff = 	1.097502506
        if (JELEM	.eq.	390	)		coeff = 	0.831931226
        if (JELEM	.eq.	391	)		coeff = 	1.063472289
        if (JELEM	.eq.	392	)		coeff = 	0.929894414
        if (JELEM	.eq.	393	)		coeff = 	0.956024322
        if (JELEM	.eq.	394	)		coeff = 	1.05162828
        if (JELEM	.eq.	395	)		coeff = 	0.87881228
        if (JELEM	.eq.	396	)		coeff = 	0.930412936
        if (JELEM	.eq.	397	)		coeff = 	1.065357964
        if (JELEM	.eq.	398	)		coeff = 	0.976487466
        if (JELEM	.eq.	399	)		coeff = 	0.908390352
        if (JELEM	.eq.	400	)		coeff = 	1.024066889
        if (JELEM	.eq.	401	)		coeff = 	0.927321556
        if (JELEM	.eq.	402	)		coeff = 	0.754816813
        if (JELEM	.eq.	403	)		coeff = 	0.944963796
        if (JELEM	.eq.	404	)		coeff = 	0.904511989
        if (JELEM	.eq.	405	)		coeff = 	0.681613431
        if (JELEM	.eq.	406	)		coeff = 	0.799227382
        if (JELEM	.eq.	407	)		coeff = 	0.919731906
        if (JELEM	.eq.	408	)		coeff = 	0.666223965
        if (JELEM	.eq.	409	)		coeff = 	1.061070603
        if (JELEM	.eq.	410	)		coeff = 	0.961281729
        if (JELEM	.eq.	411	)		coeff = 	0.92653575
        if (JELEM	.eq.	412	)		coeff = 	0.860712705
        if (JELEM	.eq.	413	)		coeff = 	0.98649023
        if (JELEM	.eq.	414	)		coeff = 	1.11247751
        if (JELEM	.eq.	415	)		coeff = 	0.958654568
        if (JELEM	.eq.	416	)		coeff = 	1.072092006
        if (JELEM	.eq.	417	)		coeff = 	1.000246037
        if (JELEM	.eq.	418	)		coeff = 	0.90325943
        if (JELEM	.eq.	419	)		coeff = 	1.074078365
        if (JELEM	.eq.	420	)		coeff = 	0.867361002
        if (JELEM	.eq.	421	)		coeff = 	1.115734595
        if (JELEM	.eq.	422	)		coeff = 	0.852948913
        if (JELEM	.eq.	423	)		coeff = 	0.823130724
        if (JELEM	.eq.	424	)		coeff = 	1.049753672
        if (JELEM	.eq.	425	)		coeff = 	0.932767284
        if (JELEM	.eq.	426	)		coeff = 	0.56518321
        if (JELEM	.eq.	427	)		coeff = 	1.121980658
        if (JELEM	.eq.	428	)		coeff = 	1.070402057
        if (JELEM	.eq.	429	)		coeff = 	1.041434572
        if (JELEM	.eq.	430	)		coeff = 	0.940163051
        if (JELEM	.eq.	431	)		coeff = 	0.912175713
        if (JELEM	.eq.	432	)		coeff = 	0.893438279
        if (JELEM	.eq.	433	)		coeff = 	0.962775978
        if (JELEM	.eq.	434	)		coeff = 	0.920988492
        if (JELEM	.eq.	435	)		coeff = 	1.127289943
        if (JELEM	.eq.	436	)		coeff = 	1.050756029
        if (JELEM	.eq.	437	)		coeff = 	0.883044627
        if (JELEM	.eq.	438	)		coeff = 	1.016437009
        if (JELEM	.eq.	439	)		coeff = 	0.883246183
        if (JELEM	.eq.	440	)		coeff = 	0.832618201
        if (JELEM	.eq.	441	)		coeff = 	1.087284899
        if (JELEM	.eq.	442	)		coeff = 	0.935710568
        if (JELEM	.eq.	443	)		coeff = 	0.894423902
        if (JELEM	.eq.	444	)		coeff = 	0.944402222
        if (JELEM	.eq.	445	)		coeff = 	1.03769064
        if (JELEM	.eq.	446	)		coeff = 	0.789101824
        if (JELEM	.eq.	447	)		coeff = 	0.742079339
        if (JELEM	.eq.	448	)		coeff = 	0.959450604
        if (JELEM	.eq.	449	)		coeff = 	0.899462065
        if (JELEM	.eq.	450	)		coeff = 	0.861292721
        if (JELEM	.eq.	451	)		coeff = 	1.048840252
        if (JELEM	.eq.	452	)		coeff = 	1.083815354
        if (JELEM	.eq.	453	)		coeff = 	0.929297164
        if (JELEM	.eq.	454	)		coeff = 	1.143019277
        if (JELEM	.eq.	455	)		coeff = 	1.126187487
        if (JELEM	.eq.	456	)		coeff = 	1.080036857
        if (JELEM	.eq.	457	)		coeff = 	0.851357345
        if (JELEM	.eq.	458	)		coeff = 	0.997550928
        if (JELEM	.eq.	459	)		coeff = 	1.039514664
        if (JELEM	.eq.	460	)		coeff = 	0.999929199
        if (JELEM	.eq.	461	)		coeff = 	1.057290783
        if (JELEM	.eq.	462	)		coeff = 	0.87233659
        if (JELEM	.eq.	463	)		coeff = 	1.010360378
        if (JELEM	.eq.	464	)		coeff = 	1.014312952
        if (JELEM	.eq.	465	)		coeff = 	0.892907141
        if (JELEM	.eq.	466	)		coeff = 	0.93958123
        if (JELEM	.eq.	467	)		coeff = 	0.87101711
        if (JELEM	.eq.	468	)		coeff = 	1.002259746
        if (JELEM	.eq.	469	)		coeff = 	1.053801825
        if (JELEM	.eq.	470	)		coeff = 	0.942934469
        if (JELEM	.eq.	471	)		coeff = 	0.968692318
        if (JELEM	.eq.	472	)		coeff = 	0.985446366
        if (JELEM	.eq.	473	)		coeff = 	1.013882644
        if (JELEM	.eq.	474	)		coeff = 	0.698912355
        if (JELEM	.eq.	475	)		coeff = 	1.030765137
        if (JELEM	.eq.	476	)		coeff = 	0.944876697
        if (JELEM	.eq.	477	)		coeff = 	1.004308888
        if (JELEM	.eq.	478	)		coeff = 	1.015359441
        if (JELEM	.eq.	479	)		coeff = 	1.166558648
        if (JELEM	.eq.	480	)		coeff = 	0.913188148
        if (JELEM	.eq.	481	)		coeff = 	0.809876154
        if (JELEM	.eq.	482	)		coeff = 	0.985530047
        if (JELEM	.eq.	483	)		coeff = 	0.828552373
        if (JELEM	.eq.	484	)		coeff = 	0.746038428
        if (JELEM	.eq.	485	)		coeff = 	0.997765893
        if (JELEM	.eq.	486	)		coeff = 	1.047669886
        if (JELEM	.eq.	487	)		coeff = 	1.055694536
        if (JELEM	.eq.	488	)		coeff = 	1.016965378
        if (JELEM	.eq.	489	)		coeff = 	0.833666278
        if (JELEM	.eq.	490	)		coeff = 	0.956831091
        if (JELEM	.eq.	491	)		coeff = 	0.846182115
        if (JELEM	.eq.	492	)		coeff = 	1.083799052
        if (JELEM	.eq.	493	)		coeff = 	0.929744844
        if (JELEM	.eq.	494	)		coeff = 	0.944956989
        if (JELEM	.eq.	495	)		coeff = 	1.029966513
        if (JELEM	.eq.	496	)		coeff = 	1.19764569
        if (JELEM	.eq.	497	)		coeff = 	1.022351039
        if (JELEM	.eq.	498	)		coeff = 	0.850012456
        if (JELEM	.eq.	499	)		coeff = 	0.802518412
        if (JELEM	.eq.	500	)		coeff = 	0.578541148
        if (JELEM	.eq.	501	)		coeff = 	1.169222045
        if (JELEM	.eq.	502	)		coeff = 	0.889715957
        if (JELEM	.eq.	503	)		coeff = 	1.094618942
        if (JELEM	.eq.	504	)		coeff = 	0.924626566
        if (JELEM	.eq.	505	)		coeff = 	0.858236654
        if (JELEM	.eq.	506	)		coeff = 	0.845352886
        if (JELEM	.eq.	507	)		coeff = 	0.94484099
        if (JELEM	.eq.	508	)		coeff = 	1.03022896
        if (JELEM	.eq.	509	)		coeff = 	0.867963026
        if (JELEM	.eq.	510	)		coeff = 	1.059156169
        if (JELEM	.eq.	511	)		coeff = 	1.026083835
        if (JELEM	.eq.	512	)		coeff = 	1.181881021
        if (JELEM	.eq.	513	)		coeff = 	0.991955582
        if (JELEM	.eq.	514	)		coeff = 	1.137561046
        if (JELEM	.eq.	515	)		coeff = 	0.830981942
        if (JELEM	.eq.	516	)		coeff = 	1.047153691
        if (JELEM	.eq.	517	)		coeff = 	0.827426083
        if (JELEM	.eq.	518	)		coeff = 	1.111583137
        if (JELEM	.eq.	519	)		coeff = 	1.041261049
        if (JELEM	.eq.	520	)		coeff = 	1.119603788
        if (JELEM	.eq.	521	)		coeff = 	0.882452471
        if (JELEM	.eq.	522	)		coeff = 	0.902743774
        if (JELEM	.eq.	523	)		coeff = 	0.781246951
        if (JELEM	.eq.	524	)		coeff = 	0.790675377
        if (JELEM	.eq.	525	)		coeff = 	0.858021388
        if (JELEM	.eq.	526	)		coeff = 	0.447576558
        if (JELEM	.eq.	527	)		coeff = 	0.860343142
        if (JELEM	.eq.	528	)		coeff = 	1.012329458
        if (JELEM	.eq.	529	)		coeff = 	1.139256464
        if (JELEM	.eq.	530	)		coeff = 	0.889723677
        if (JELEM	.eq.	531	)		coeff = 	0.971526416
        if (JELEM	.eq.	532	)		coeff = 	0.993300973
        if (JELEM	.eq.	533	)		coeff = 	0.827641061
        if (JELEM	.eq.	534	)		coeff = 	1.110671208
        if (JELEM	.eq.	535	)		coeff = 	0.995478513
        if (JELEM	.eq.	536	)		coeff = 	1.069230639
        if (JELEM	.eq.	537	)		coeff = 	1.119284512
        if (JELEM	.eq.	538	)		coeff = 	0.902765952
        if (JELEM	.eq.	539	)		coeff = 	0.786403301
        if (JELEM	.eq.	540	)		coeff = 	0.995749273
        if (JELEM	.eq.	541	)		coeff = 	1.035266706
        if (JELEM	.eq.	542	)		coeff = 	0.792323842
        if (JELEM	.eq.	543	)		coeff = 	0.901312379
        if (JELEM	.eq.	544	)		coeff = 	0.802094767
        if (JELEM	.eq.	545	)		coeff = 	0.74606632
        if (JELEM	.eq.	546	)		coeff = 	0.761304345
        if (JELEM	.eq.	547	)		coeff = 	0.563333957
        if (JELEM	.eq.	548	)		coeff = 	0.992122962
        if (JELEM	.eq.	549	)		coeff = 	0.915943701
        if (JELEM	.eq.	550	)		coeff = 	0.849126025
        if (JELEM	.eq.	551	)		coeff = 	0.928107996
        if (JELEM	.eq.	552	)		coeff = 	0.901860012
        if (JELEM	.eq.	553	)		coeff = 	0.966919891
        if (JELEM	.eq.	554	)		coeff = 	0.961553132
        if (JELEM	.eq.	555	)		coeff = 	0.870141934
        if (JELEM	.eq.	556	)		coeff = 	1.123588407
        if (JELEM	.eq.	557	)		coeff = 	0.797634982
        if (JELEM	.eq.	558	)		coeff = 	0.890027142
        if (JELEM	.eq.	559	)		coeff = 	0.93016862
        if (JELEM	.eq.	560	)		coeff = 	0.880720146
        if (JELEM	.eq.	561	)		coeff = 	1.023147562
        if (JELEM	.eq.	562	)		coeff = 	1.046716694
        if (JELEM	.eq.	563	)		coeff = 	0.918514481
        if (JELEM	.eq.	564	)		coeff = 	0.771383901
        if (JELEM	.eq.	565	)		coeff = 	0.858488596
        if (JELEM	.eq.	566	)		coeff = 	1.004754669
        if (JELEM	.eq.	567	)		coeff = 	1.02726355
        if (JELEM	.eq.	568	)		coeff = 	0.885044626
        if (JELEM	.eq.	569	)		coeff = 	0.893515979
        if (JELEM	.eq.	570	)		coeff = 	0.829723191
        if (JELEM	.eq.	571	)		coeff = 	0.781215389
        if (JELEM	.eq.	572	)		coeff = 	1.027606625
        if (JELEM	.eq.	573	)		coeff = 	0.906270072
        if (JELEM	.eq.	574	)		coeff = 	0.929477803
        if (JELEM	.eq.	575	)		coeff = 	1.101833879
        if (JELEM	.eq.	576	)		coeff = 	1.052569465
        if (JELEM	.eq.	577	)		coeff = 	1.008833121
        if (JELEM	.eq.	578	)		coeff = 	1.01563441
        if (JELEM	.eq.	579	)		coeff = 	0.985831424
        if (JELEM	.eq.	580	)		coeff = 	1.09300259
        if (JELEM	.eq.	581	)		coeff = 	0.946993059
        if (JELEM	.eq.	582	)		coeff = 	0.91335829
        if (JELEM	.eq.	583	)		coeff = 	1.080063914
        if (JELEM	.eq.	584	)		coeff = 	1.014346322
        if (JELEM	.eq.	585	)		coeff = 	0.674153468
        if (JELEM	.eq.	586	)		coeff = 	0.809739075
        if (JELEM	.eq.	587	)		coeff = 	1.142303781
        if (JELEM	.eq.	588	)		coeff = 	0.825533369
        if (JELEM	.eq.	589	)		coeff = 	1.034700567
        if (JELEM	.eq.	590	)		coeff = 	0.982203194
        if (JELEM	.eq.	591	)		coeff = 	1.071757261
        if (JELEM	.eq.	592	)		coeff = 	1.043738949
        if (JELEM	.eq.	593	)		coeff = 	0.872358969
        if (JELEM	.eq.	594	)		coeff = 	1.093286019
        if (JELEM	.eq.	595	)		coeff = 	0.991193757
        if (JELEM	.eq.	596	)		coeff = 	0.938032081
        if (JELEM	.eq.	597	)		coeff = 	1.049876059
        if (JELEM	.eq.	598	)		coeff = 	0.903190704
        if (JELEM	.eq.	599	)		coeff = 	0.937169163
        if (JELEM	.eq.	600	)		coeff = 	0.973880462
        if (JELEM	.eq.	601	)		coeff = 	1.035856505
        if (JELEM	.eq.	602	)		coeff = 	0.987842337
        if (JELEM	.eq.	603	)		coeff = 	0.87421429
        if (JELEM	.eq.	604	)		coeff = 	1.128263828
        if (JELEM	.eq.	605	)		coeff = 	0.932719076
        if (JELEM	.eq.	606	)		coeff = 	1.097367329
        if (JELEM	.eq.	607	)		coeff = 	0.987137475
        if (JELEM	.eq.	608	)		coeff = 	0.655422258
        if (JELEM	.eq.	609	)		coeff = 	0.838207916
        if (JELEM	.eq.	610	)		coeff = 	1.114203498
        if (JELEM	.eq.	611	)		coeff = 	0.891781875
        if (JELEM	.eq.	612	)		coeff = 	1.080213076
        if (JELEM	.eq.	613	)		coeff = 	0.968061813
        if (JELEM	.eq.	614	)		coeff = 	0.822024268
        if (JELEM	.eq.	615	)		coeff = 	1.151316624
        if (JELEM	.eq.	616	)		coeff = 	0.851247538
        if (JELEM	.eq.	617	)		coeff = 	0.868262845
        if (JELEM	.eq.	618	)		coeff = 	1.050803643
        if (JELEM	.eq.	619	)		coeff = 	0.787036189
        if (JELEM	.eq.	620	)		coeff = 	0.907835183
        if (JELEM	.eq.	621	)		coeff = 	0.80781765
        if (JELEM	.eq.	622	)		coeff = 	0.897326102
        if (JELEM	.eq.	623	)		coeff = 	0.968196706
        if (JELEM	.eq.	624	)		coeff = 	0.843561493
        if (JELEM	.eq.	625	)		coeff = 	1.091479156
        if (JELEM	.eq.	626	)		coeff = 	1.138831595
        if (JELEM	.eq.	627	)		coeff = 	1.008789511
        if (JELEM	.eq.	628	)		coeff = 	0.967636264
        if (JELEM	.eq.	629	)		coeff = 	1.040686445
        if (JELEM	.eq.	630	)		coeff = 	0.948493868
        if (JELEM	.eq.	631	)		coeff = 	0.908720458
        if (JELEM	.eq.	632	)		coeff = 	0.842600328
        if (JELEM	.eq.	633	)		coeff = 	0.801903704
        if (JELEM	.eq.	634	)		coeff = 	0.981119505
        if (JELEM	.eq.	635	)		coeff = 	0.991389177
        if (JELEM	.eq.	636	)		coeff = 	0.980435132
        if (JELEM	.eq.	637	)		coeff = 	1.016031503
        if (JELEM	.eq.	638	)		coeff = 	0.978578847
        if (JELEM	.eq.	639	)		coeff = 	0.717698486
        if (JELEM	.eq.	640	)		coeff = 	0.974383771
        if (JELEM	.eq.	641	)		coeff = 	0.965125764
        if (JELEM	.eq.	642	)		coeff = 	1.143954586
        if (JELEM	.eq.	643	)		coeff = 	1.051274583
        if (JELEM	.eq.	644	)		coeff = 	0.800045831
        if (JELEM	.eq.	645	)		coeff = 	0.994193085
        if (JELEM	.eq.	646	)		coeff = 	1.086020055
        if (JELEM	.eq.	647	)		coeff = 	1.004956351
        if (JELEM	.eq.	648	)		coeff = 	1.055908343
        if (JELEM	.eq.	649	)		coeff = 	1.033473877
        if (JELEM	.eq.	650	)		coeff = 	0.922171304
        if (JELEM	.eq.	651	)		coeff = 	0.690232831
        if (JELEM	.eq.	652	)		coeff = 	0.779344892
        if (JELEM	.eq.	653	)		coeff = 	1.11320018
        if (JELEM	.eq.	654	)		coeff = 	0.795855126
        if (JELEM	.eq.	655	)		coeff = 	1.052034592
        if (JELEM	.eq.	656	)		coeff = 	0.612342978
        if (JELEM	.eq.	657	)		coeff = 	1.011336499
        if (JELEM	.eq.	658	)		coeff = 	1.040011914
        if (JELEM	.eq.	659	)		coeff = 	0.752132827
        if (JELEM	.eq.	660	)		coeff = 	0.855371712
        if (JELEM	.eq.	661	)		coeff = 	0.77539251
        if (JELEM	.eq.	662	)		coeff = 	0.895217471
        if (JELEM	.eq.	663	)		coeff = 	0.814611645
        if (JELEM	.eq.	664	)		coeff = 	0.824371806
        if (JELEM	.eq.	665	)		coeff = 	1.002846216
        if (JELEM	.eq.	666	)		coeff = 	1.071401491
        if (JELEM	.eq.	667	)		coeff = 	0.857894963
        if (JELEM	.eq.	668	)		coeff = 	0.908728639
        if (JELEM	.eq.	669	)		coeff = 	0.982803416
        if (JELEM	.eq.	670	)		coeff = 	0.747321401
        if (JELEM	.eq.	671	)		coeff = 	0.741222739
        if (JELEM	.eq.	672	)		coeff = 	1.13148702
        if (JELEM	.eq.	673	)		coeff = 	1.025532507
        if (JELEM	.eq.	674	)		coeff = 	1.104093039
        if (JELEM	.eq.	675	)		coeff = 	0.8469896
        if (JELEM	.eq.	676	)		coeff = 	0.992290418
        if (JELEM	.eq.	677	)		coeff = 	0.940710552
        if (JELEM	.eq.	678	)		coeff = 	1.006644735
        if (JELEM	.eq.	679	)		coeff = 	0.702356098
        if (JELEM	.eq.	680	)		coeff = 	0.933992425
        if (JELEM	.eq.	681	)		coeff = 	0.913352906
        if (JELEM	.eq.	682	)		coeff = 	0.825507571
        if (JELEM	.eq.	683	)		coeff = 	0.985863075
        if (JELEM	.eq.	684	)		coeff = 	1.100758431
        if (JELEM	.eq.	685	)		coeff = 	0.820124906
        if (JELEM	.eq.	686	)		coeff = 	0.903963324
        if (JELEM	.eq.	687	)		coeff = 	0.838916745
        if (JELEM	.eq.	688	)		coeff = 	0.910415683
        if (JELEM	.eq.	689	)		coeff = 	0.901535841
        if (JELEM	.eq.	690	)		coeff = 	0.802273908
        if (JELEM	.eq.	691	)		coeff = 	0.938843597
        if (JELEM	.eq.	692	)		coeff = 	1.011697562
        if (JELEM	.eq.	693	)		coeff = 	0.590591694
        if (JELEM	.eq.	694	)		coeff = 	1.046302768
        if (JELEM	.eq.	695	)		coeff = 	1.022637196
        if (JELEM	.eq.	696	)		coeff = 	0.843622349
        if (JELEM	.eq.	697	)		coeff = 	0.874321617
        if (JELEM	.eq.	698	)		coeff = 	1.048550452
        if (JELEM	.eq.	699	)		coeff = 	0.876684738
        if (JELEM	.eq.	700	)		coeff = 	0.930124898
        if (JELEM	.eq.	701	)		coeff = 	0.717949954
        if (JELEM	.eq.	702	)		coeff = 	0.962984346
        if (JELEM	.eq.	703	)		coeff = 	0.969901728
        if (JELEM	.eq.	704	)		coeff = 	0.939358907
        if (JELEM	.eq.	705	)		coeff = 	0.678309036
        if (JELEM	.eq.	706	)		coeff = 	1.005120342
        if (JELEM	.eq.	707	)		coeff = 	1.029946946
        if (JELEM	.eq.	708	)		coeff = 	0.938082873
        if (JELEM	.eq.	709	)		coeff = 	0.884735091
        if (JELEM	.eq.	710	)		coeff = 	1.096012238
        if (JELEM	.eq.	711	)		coeff = 	0.998078167
        if (JELEM	.eq.	712	)		coeff = 	1.093377472
        if (JELEM	.eq.	713	)		coeff = 	1.013652724
        if (JELEM	.eq.	714	)		coeff = 	0.976870899
        if (JELEM	.eq.	715	)		coeff = 	1.062053696
        if (JELEM	.eq.	716	)		coeff = 	1.109065664
        if (JELEM	.eq.	717	)		coeff = 	0.979146419
        if (JELEM	.eq.	718	)		coeff = 	1.022916889
        if (JELEM	.eq.	719	)		coeff = 	1.049006995
        if (JELEM	.eq.	720	)		coeff = 	0.961274722
        if (JELEM	.eq.	721	)		coeff = 	1.146918361
        if (JELEM	.eq.	722	)		coeff = 	0.845277603
        if (JELEM	.eq.	723	)		coeff = 	1.048613443
        if (JELEM	.eq.	724	)		coeff = 	0.984110098
        if (JELEM	.eq.	725	)		coeff = 	0.889620131
        if (JELEM	.eq.	726	)		coeff = 	0.956492936
        if (JELEM	.eq.	727	)		coeff = 	1.027522228
        if (JELEM	.eq.	728	)		coeff = 	0.875366682
        if (JELEM	.eq.	729	)		coeff = 	0.736027189
        if (JELEM	.eq.	730	)		coeff = 	1.079035584
        if (JELEM	.eq.	731	)		coeff = 	1.139573056
        if (JELEM	.eq.	732	)		coeff = 	0.985334547
        if (JELEM	.eq.	733	)		coeff = 	1.106167389
        if (JELEM	.eq.	734	)		coeff = 	1.075707272
        if (JELEM	.eq.	735	)		coeff = 	1.09064621
        if (JELEM	.eq.	736	)		coeff = 	1.03843013
        if (JELEM	.eq.	737	)		coeff = 	1.045135387
        if (JELEM	.eq.	738	)		coeff = 	0.93202156
        if (JELEM	.eq.	739	)		coeff = 	0.992776412
        if (JELEM	.eq.	740	)		coeff = 	0.892097333
        if (JELEM	.eq.	741	)		coeff = 	0.723814774
        if (JELEM	.eq.	742	)		coeff = 	0.882769957
        if (JELEM	.eq.	743	)		coeff = 	1.097768684
        if (JELEM	.eq.	744	)		coeff = 	0.822454232
        if (JELEM	.eq.	745	)		coeff = 	1.001389665
        if (JELEM	.eq.	746	)		coeff = 	0.862148055
        if (JELEM	.eq.	747	)		coeff = 	0.867114526
        if (JELEM	.eq.	748	)		coeff = 	1.009994844
        if (JELEM	.eq.	749	)		coeff = 	0.743042533
        if (JELEM	.eq.	750	)		coeff = 	0.846156393
        if (JELEM	.eq.	751	)		coeff = 	0.92597075
        if (JELEM	.eq.	752	)		coeff = 	1.019972696
        if (JELEM	.eq.	753	)		coeff = 	0.92978776
        if (JELEM	.eq.	754	)		coeff = 	0.640188215
        if (JELEM	.eq.	755	)		coeff = 	1.00651445
        if (JELEM	.eq.	756	)		coeff = 	1.036313138
        if (JELEM	.eq.	757	)		coeff = 	1.228523386
        if (JELEM	.eq.	758	)		coeff = 	0.790340266
        if (JELEM	.eq.	759	)		coeff = 	0.926474914
        if (JELEM	.eq.	760	)		coeff = 	0.867186955
        if (JELEM	.eq.	761	)		coeff = 	0.88338316
        if (JELEM	.eq.	762	)		coeff = 	0.929821115
        if (JELEM	.eq.	763	)		coeff = 	1.054252304
        if (JELEM	.eq.	764	)		coeff = 	0.837657651
        if (JELEM	.eq.	765	)		coeff = 	0.896160287
        if (JELEM	.eq.	766	)		coeff = 	0.918641257
        if (JELEM	.eq.	767	)		coeff = 	0.931093864
        if (JELEM	.eq.	768	)		coeff = 	0.818478948
        if (JELEM	.eq.	769	)		coeff = 	0.864195961
        if (JELEM	.eq.	770	)		coeff = 	0.892969067
        if (JELEM	.eq.	771	)		coeff = 	0.988820267
        if (JELEM	.eq.	772	)		coeff = 	0.545877093
        if (JELEM	.eq.	773	)		coeff = 	0.927654141
        if (JELEM	.eq.	774	)		coeff = 	1.004008823
        if (JELEM	.eq.	775	)		coeff = 	1.054482441
        if (JELEM	.eq.	776	)		coeff = 	1.012887544
        if (JELEM	.eq.	777	)		coeff = 	1.031537691
        if (JELEM	.eq.	778	)		coeff = 	0.755674323
        if (JELEM	.eq.	779	)		coeff = 	0.804290096
        if (JELEM	.eq.	780	)		coeff = 	0.642672734
        if (JELEM	.eq.	781	)		coeff = 	0.739176371
        if (JELEM	.eq.	782	)		coeff = 	0.818186414
        if (JELEM	.eq.	783	)		coeff = 	1.061634637
        if (JELEM	.eq.	784	)		coeff = 	1.067376305
        if (JELEM	.eq.	785	)		coeff = 	0.919298138
        if (JELEM	.eq.	786	)		coeff = 	0.865406249
        if (JELEM	.eq.	787	)		coeff = 	1.101398552
        if (JELEM	.eq.	788	)		coeff = 	1.034842999
        if (JELEM	.eq.	789	)		coeff = 	1.163256206
        if (JELEM	.eq.	790	)		coeff = 	1.038475393
        if (JELEM	.eq.	791	)		coeff = 	1.031760243
        if (JELEM	.eq.	792	)		coeff = 	0.787654764
        if (JELEM	.eq.	793	)		coeff = 	1.000090487
        if (JELEM	.eq.	794	)		coeff = 	1.125758239
        if (JELEM	.eq.	795	)		coeff = 	0.918960079
        if (JELEM	.eq.	796	)		coeff = 	0.852562786
        if (JELEM	.eq.	797	)		coeff = 	1.0950745
        if (JELEM	.eq.	798	)		coeff = 	1.161454379
        if (JELEM	.eq.	799	)		coeff = 	1.000615645
        if (JELEM	.eq.	800	)		coeff = 	0.987843692
        if (JELEM	.eq.	801	)		coeff = 	0.986963128
        if (JELEM	.eq.	802	)		coeff = 	0.974102942
        if (JELEM	.eq.	803	)		coeff = 	0.973383703
        if (JELEM	.eq.	804	)		coeff = 	0.950388622
        if (JELEM	.eq.	805	)		coeff = 	0.954819378
        if (JELEM	.eq.	806	)		coeff = 	0.961867459
        if (JELEM	.eq.	807	)		coeff = 	0.869522463
        if (JELEM	.eq.	808	)		coeff = 	1.021578803
        if (JELEM	.eq.	809	)		coeff = 	0.949712097
        if (JELEM	.eq.	810	)		coeff = 	0.942645555
        if (JELEM	.eq.	811	)		coeff = 	1.047476481
        if (JELEM	.eq.	812	)		coeff = 	0.94834058
        if (JELEM	.eq.	813	)		coeff = 	1.04420784
        if (JELEM	.eq.	814	)		coeff = 	0.950100622
        if (JELEM	.eq.	815	)		coeff = 	1.102883999
        if (JELEM	.eq.	816	)		coeff = 	0.827364088
        if (JELEM	.eq.	817	)		coeff = 	1.046032313
        if (JELEM	.eq.	818	)		coeff = 	1.036619884
        if (JELEM	.eq.	819	)		coeff = 	0.996551491
        if (JELEM	.eq.	820	)		coeff = 	1.022319276
        if (JELEM	.eq.	821	)		coeff = 	0.882644646
        if (JELEM	.eq.	822	)		coeff = 	1.013433678
        if (JELEM	.eq.	823	)		coeff = 	1.083575685
        if (JELEM	.eq.	824	)		coeff = 	0.857719587
        if (JELEM	.eq.	825	)		coeff = 	1.032589946
        if (JELEM	.eq.	826	)		coeff = 	0.962924797
        if (JELEM	.eq.	827	)		coeff = 	1.124141031
        if (JELEM	.eq.	828	)		coeff = 	0.927532707
        if (JELEM	.eq.	829	)		coeff = 	0.994694103
        if (JELEM	.eq.	830	)		coeff = 	1.016410258
        if (JELEM	.eq.	831	)		coeff = 	1.009803035
        if (JELEM	.eq.	832	)		coeff = 	0.767431279
        if (JELEM	.eq.	833	)		coeff = 	0.972841013
        if (JELEM	.eq.	834	)		coeff = 	1.253180544
        if (JELEM	.eq.	835	)		coeff = 	0.892449888
        if (JELEM	.eq.	836	)		coeff = 	0.991435798
        if (JELEM	.eq.	837	)		coeff = 	1.14208505
        if (JELEM	.eq.	838	)		coeff = 	0.923217353
        if (JELEM	.eq.	839	)		coeff = 	0.980643557
        if (JELEM	.eq.	840	)		coeff = 	0.937531852
        if (JELEM	.eq.	841	)		coeff = 	0.929481691
        if (JELEM	.eq.	842	)		coeff = 	0.976634804
        if (JELEM	.eq.	843	)		coeff = 	1.116650209
        if (JELEM	.eq.	844	)		coeff = 	0.894767509
        if (JELEM	.eq.	845	)		coeff = 	0.886857223
        if (JELEM	.eq.	846	)		coeff = 	0.989998146
        if (JELEM	.eq.	847	)		coeff = 	1.114922135
        if (JELEM	.eq.	848	)		coeff = 	0.746020506
        if (JELEM	.eq.	849	)		coeff = 	0.815949331
        if (JELEM	.eq.	850	)		coeff = 	0.846611327
        if (JELEM	.eq.	851	)		coeff = 	0.978121113
        if (JELEM	.eq.	852	)		coeff = 	0.810785861
        if (JELEM	.eq.	853	)		coeff = 	0.831587602
        if (JELEM	.eq.	854	)		coeff = 	0.919975105
        if (JELEM	.eq.	855	)		coeff = 	0.962121519
        if (JELEM	.eq.	856	)		coeff = 	1.099275964
        if (JELEM	.eq.	857	)		coeff = 	1.001488757
        if (JELEM	.eq.	858	)		coeff = 	1.038919921
        if (JELEM	.eq.	859	)		coeff = 	0.827534775
        if (JELEM	.eq.	860	)		coeff = 	0.9242934
        if (JELEM	.eq.	861	)		coeff = 	0.890256402
        if (JELEM	.eq.	862	)		coeff = 	1.063773109
        if (JELEM	.eq.	863	)		coeff = 	0.743566217
        if (JELEM	.eq.	864	)		coeff = 	0.991181551
        if (JELEM	.eq.	865	)		coeff = 	0.931534139
        if (JELEM	.eq.	866	)		coeff = 	0.861439486
        if (JELEM	.eq.	867	)		coeff = 	0.906896839
        if (JELEM	.eq.	868	)		coeff = 	1.258132018
        if (JELEM	.eq.	869	)		coeff = 	0.963272403
        if (JELEM	.eq.	870	)		coeff = 	0.65477192
        if (JELEM	.eq.	871	)		coeff = 	1.17370693
        if (JELEM	.eq.	872	)		coeff = 	1.123206049
        if (JELEM	.eq.	873	)		coeff = 	1.085380678
        if (JELEM	.eq.	874	)		coeff = 	0.933992194
        if (JELEM	.eq.	875	)		coeff = 	1.118882061
        if (JELEM	.eq.	876	)		coeff = 	0.98538006
        if (JELEM	.eq.	877	)		coeff = 	1.003272307
        if (JELEM	.eq.	878	)		coeff = 	0.902882561
        if (JELEM	.eq.	879	)		coeff = 	1.001333352
        if (JELEM	.eq.	880	)		coeff = 	1.011874931
        if (JELEM	.eq.	881	)		coeff = 	0.909991917
        if (JELEM	.eq.	882	)		coeff = 	0.940090745
        if (JELEM	.eq.	883	)		coeff = 	1.036068064
        if (JELEM	.eq.	884	)		coeff = 	0.998102392
        if (JELEM	.eq.	885	)		coeff = 	1.047204339
        if (JELEM	.eq.	886	)		coeff = 	1.030641507
        if (JELEM	.eq.	887	)		coeff = 	0.835464281
        if (JELEM	.eq.	888	)		coeff = 	0.983980736
        if (JELEM	.eq.	889	)		coeff = 	0.916927457
        if (JELEM	.eq.	890	)		coeff = 	0.955223044
        if (JELEM	.eq.	891	)		coeff = 	0.76089696
        if (JELEM	.eq.	892	)		coeff = 	1.00942459
        if (JELEM	.eq.	893	)		coeff = 	1.001895901
        if (JELEM	.eq.	894	)		coeff = 	0.884286573
        if (JELEM	.eq.	895	)		coeff = 	0.825113681
        if (JELEM	.eq.	896	)		coeff = 	1.143249876
        if (JELEM	.eq.	897	)		coeff = 	0.834141733
        if (JELEM	.eq.	898	)		coeff = 	0.856981988
        if (JELEM	.eq.	899	)		coeff = 	0.951935174
        if (JELEM	.eq.	900	)		coeff = 	0.982154034
        if (JELEM	.eq.	901	)		coeff = 	0.934289113
        if (JELEM	.eq.	902	)		coeff = 	1.00896096
        if (JELEM	.eq.	903	)		coeff = 	0.87086032
        if (JELEM	.eq.	904	)		coeff = 	1.08854434
        if (JELEM	.eq.	905	)		coeff = 	0.977218087
        if (JELEM	.eq.	906	)		coeff = 	0.877883767
        if (JELEM	.eq.	907	)		coeff = 	0.927015168
        if (JELEM	.eq.	908	)		coeff = 	1.205992119
        if (JELEM	.eq.	909	)		coeff = 	1.150535954
        if (JELEM	.eq.	910	)		coeff = 	0.975883268
        if (JELEM	.eq.	911	)		coeff = 	1.041867714
        if (JELEM	.eq.	912	)		coeff = 	0.999037838
        if (JELEM	.eq.	913	)		coeff = 	0.980133206
        if (JELEM	.eq.	914	)		coeff = 	1.108527731
        if (JELEM	.eq.	915	)		coeff = 	1.001043586
        if (JELEM	.eq.	916	)		coeff = 	0.997559107
        if (JELEM	.eq.	917	)		coeff = 	1.018386408
        if (JELEM	.eq.	918	)		coeff = 	0.781646381
        if (JELEM	.eq.	919	)		coeff = 	0.970273223
        if (JELEM	.eq.	920	)		coeff = 	1.055276134
        if (JELEM	.eq.	921	)		coeff = 	1.019572617
        if (JELEM	.eq.	922	)		coeff = 	0.845715984
        if (JELEM	.eq.	923	)		coeff = 	0.929438129
        if (JELEM	.eq.	924	)		coeff = 	0.914716428
        if (JELEM	.eq.	925	)		coeff = 	0.930037627
        if (JELEM	.eq.	926	)		coeff = 	1.029869411
        if (JELEM	.eq.	927	)		coeff = 	0.819369053
        if (JELEM	.eq.	928	)		coeff = 	1.179852031
        if (JELEM	.eq.	929	)		coeff = 	1.000713578
        if (JELEM	.eq.	930	)		coeff = 	0.937470634
        if (JELEM	.eq.	931	)		coeff = 	0.774969955
        if (JELEM	.eq.	932	)		coeff = 	0.97556941
        if (JELEM	.eq.	933	)		coeff = 	0.9847477
        if (JELEM	.eq.	934	)		coeff = 	0.917597824
        if (JELEM	.eq.	935	)		coeff = 	0.876911094
        if (JELEM	.eq.	936	)		coeff = 	0.996806748
        if (JELEM	.eq.	937	)		coeff = 	0.976642352
        if (JELEM	.eq.	938	)		coeff = 	0.921257813
        if (JELEM	.eq.	939	)		coeff = 	0.888810464
        if (JELEM	.eq.	940	)		coeff = 	0.924475324
        if (JELEM	.eq.	941	)		coeff = 	0.868684933
        if (JELEM	.eq.	942	)		coeff = 	1.129084729
        if (JELEM	.eq.	943	)		coeff = 	1.086574018
        if (JELEM	.eq.	944	)		coeff = 	1.122514173
        if (JELEM	.eq.	945	)		coeff = 	1.062256993
        if (JELEM	.eq.	946	)		coeff = 	0.842566904
        if (JELEM	.eq.	947	)		coeff = 	0.947119705
        if (JELEM	.eq.	948	)		coeff = 	1.058167915
        if (JELEM	.eq.	949	)		coeff = 	0.758647441
        if (JELEM	.eq.	950	)		coeff = 	0.952114114
        if (JELEM	.eq.	951	)		coeff = 	0.845192667
        if (JELEM	.eq.	952	)		coeff = 	1.013849638
        if (JELEM	.eq.	953	)		coeff = 	1.038414801
        if (JELEM	.eq.	954	)		coeff = 	0.901321124
        if (JELEM	.eq.	955	)		coeff = 	0.958474406
        if (JELEM	.eq.	956	)		coeff = 	1.022667564
        if (JELEM	.eq.	957	)		coeff = 	1.136257541
        if (JELEM	.eq.	958	)		coeff = 	1.065083312
        if (JELEM	.eq.	959	)		coeff = 	0.868240583
        if (JELEM	.eq.	960	)		coeff = 	0.958765056
        if (JELEM	.eq.	961	)		coeff = 	0.986746211
        if (JELEM	.eq.	962	)		coeff = 	1.083044711
        if (JELEM	.eq.	963	)		coeff = 	0.925239659
        if (JELEM	.eq.	964	)		coeff = 	1.075241599
        if (JELEM	.eq.	965	)		coeff = 	0.940869488
        if (JELEM	.eq.	966	)		coeff = 	1.053880424
        if (JELEM	.eq.	967	)		coeff = 	0.965425633
        if (JELEM	.eq.	968	)		coeff = 	1.008868597
        if (JELEM	.eq.	969	)		coeff = 	1.046524041
        if (JELEM	.eq.	970	)		coeff = 	1.084578714
        if (JELEM	.eq.	971	)		coeff = 	1.064775928
        if (JELEM	.eq.	972	)		coeff = 	0.818795473
        if (JELEM	.eq.	973	)		coeff = 	0.892989514
        if (JELEM	.eq.	974	)		coeff = 	1.18866726
        if (JELEM	.eq.	975	)		coeff = 	0.855588684
        if (JELEM	.eq.	976	)		coeff = 	1.034511848
        if (JELEM	.eq.	977	)		coeff = 	0.809688519
        if (JELEM	.eq.	978	)		coeff = 	0.840148236
        if (JELEM	.eq.	979	)		coeff = 	1.003353003
        if (JELEM	.eq.	980	)		coeff = 	1.014461509
        if (JELEM	.eq.	981	)		coeff = 	1.064571221
        if (JELEM	.eq.	982	)		coeff = 	0.994428333
        if (JELEM	.eq.	983	)		coeff = 	1.017205167
        if (JELEM	.eq.	984	)		coeff = 	0.850061997
        if (JELEM	.eq.	985	)		coeff = 	1.024462358
        if (JELEM	.eq.	986	)		coeff = 	0.83297366
        if (JELEM	.eq.	987	)		coeff = 	1.050283272
        if (JELEM	.eq.	988	)		coeff = 	1.003507182
        if (JELEM	.eq.	989	)		coeff = 	1.005975444
        if (JELEM	.eq.	990	)		coeff = 	0.939378544
        if (JELEM	.eq.	991	)		coeff = 	0.976047749
        if (JELEM	.eq.	992	)		coeff = 	1.125384271
        if (JELEM	.eq.	993	)		coeff = 	0.996705895
        if (JELEM	.eq.	994	)		coeff = 	0.950870302
        if (JELEM	.eq.	995	)		coeff = 	1.103788051
        if (JELEM	.eq.	996	)		coeff = 	1.044489561
        if (JELEM	.eq.	997	)		coeff = 	0.873602321
        if (JELEM	.eq.	998	)		coeff = 	0.925573494
        if (JELEM	.eq.	999	)		coeff = 	0.856436509
        if (JELEM	.eq.	1000	)		coeff = 	1.144023678
        if (JELEM	.eq.	1001	)		coeff = 	1.160835605
        if (JELEM	.eq.	1002	)		coeff = 	0.991621413
        if (JELEM	.eq.	1003	)		coeff = 	1.045612645
        if (JELEM	.eq.	1004	)		coeff = 	0.938693103
        if (JELEM	.eq.	1005	)		coeff = 	0.849246907
        if (JELEM	.eq.	1006	)		coeff = 	0.883236458
        if (JELEM	.eq.	1007	)		coeff = 	0.770234836
        if (JELEM	.eq.	1008	)		coeff = 	0.957130128
        if (JELEM	.eq.	1009	)		coeff = 	0.953120599
        if (JELEM	.eq.	1010	)		coeff = 	0.913470251
        if (JELEM	.eq.	1011	)		coeff = 	0.820217754
        if (JELEM	.eq.	1012	)		coeff = 	1.054958007
        if (JELEM	.eq.	1013	)		coeff = 	0.902365496
        if (JELEM	.eq.	1014	)		coeff = 	1.052398212
        if (JELEM	.eq.	1015	)		coeff = 	1.048798902
        if (JELEM	.eq.	1016	)		coeff = 	0.874008958
        if (JELEM	.eq.	1017	)		coeff = 	0.95561387
        if (JELEM	.eq.	1018	)		coeff = 	1.182749682
        if (JELEM	.eq.	1019	)		coeff = 	0.734444946
        if (JELEM	.eq.	1020	)		coeff = 	0.816025979
        if (JELEM	.eq.	1021	)		coeff = 	0.81949967
        if (JELEM	.eq.	1022	)		coeff = 	0.965671249
        if (JELEM	.eq.	1023	)		coeff = 	0.955271384
        if (JELEM	.eq.	1024	)		coeff = 	1.131302709
        if (JELEM	.eq.	1025	)		coeff = 	0.882794724
        if (JELEM	.eq.	1026	)		coeff = 	0.879399899
        if (JELEM	.eq.	1027	)		coeff = 	0.93999794
        if (JELEM	.eq.	1028	)		coeff = 	0.999391947
        if (JELEM	.eq.	1029	)		coeff = 	0.966232339
        if (JELEM	.eq.	1030	)		coeff = 	0.956383662
        if (JELEM	.eq.	1031	)		coeff = 	0.978097211
        if (JELEM	.eq.	1032	)		coeff = 	1.069014537
        if (JELEM	.eq.	1033	)		coeff = 	0.918330692
        if (JELEM	.eq.	1034	)		coeff = 	0.910957076
        if (JELEM	.eq.	1035	)		coeff = 	0.882955896
        if (JELEM	.eq.	1036	)		coeff = 	1.017068722
        if (JELEM	.eq.	1037	)		coeff = 	0.85859941
        if (JELEM	.eq.	1038	)		coeff = 	1.018745398
        if (JELEM	.eq.	1039	)		coeff = 	1.026133487
        if (JELEM	.eq.	1040	)		coeff = 	1.027177399
        if (JELEM	.eq.	1041	)		coeff = 	0.877228821
        if (JELEM	.eq.	1042	)		coeff = 	1.013534835
        if (JELEM	.eq.	1043	)		coeff = 	1.002960938
        if (JELEM	.eq.	1044	)		coeff = 	0.9106784
        if (JELEM	.eq.	1045	)		coeff = 	1.130016712
        if (JELEM	.eq.	1046	)		coeff = 	1.066119874
        if (JELEM	.eq.	1047	)		coeff = 	0.985467155
        if (JELEM	.eq.	1048	)		coeff = 	1.044594932
        if (JELEM	.eq.	1049	)		coeff = 	0.878676926
        if (JELEM	.eq.	1050	)		coeff = 	0.921418884
        if (JELEM	.eq.	1051	)		coeff = 	0.839938902
        if (JELEM	.eq.	1052	)		coeff = 	1.046220283
        if (JELEM	.eq.	1053	)		coeff = 	0.980916113
        if (JELEM	.eq.	1054	)		coeff = 	0.889998283
        if (JELEM	.eq.	1055	)		coeff = 	1.084768422
        if (JELEM	.eq.	1056	)		coeff = 	0.979228073
        if (JELEM	.eq.	1057	)		coeff = 	0.966632899
        if (JELEM	.eq.	1058	)		coeff = 	0.9040447
        if (JELEM	.eq.	1059	)		coeff = 	0.768339728
        if (JELEM	.eq.	1060	)		coeff = 	1.064215617
        if (JELEM	.eq.	1061	)		coeff = 	1.114002969
        if (JELEM	.eq.	1062	)		coeff = 	0.932371873
        if (JELEM	.eq.	1063	)		coeff = 	1.035317911
        if (JELEM	.eq.	1064	)		coeff = 	0.947449443
        if (JELEM	.eq.	1065	)		coeff = 	0.909638005
        if (JELEM	.eq.	1066	)		coeff = 	0.973247387
        if (JELEM	.eq.	1067	)		coeff = 	0.916454415
        if (JELEM	.eq.	1068	)		coeff = 	0.856621943
        if (JELEM	.eq.	1069	)		coeff = 	0.829899109
        if (JELEM	.eq.	1070	)		coeff = 	0.898065159
        if (JELEM	.eq.	1071	)		coeff = 	1.029263189
        if (JELEM	.eq.	1072	)		coeff = 	1.019835808
        if (JELEM	.eq.	1073	)		coeff = 	1.053426664
        if (JELEM	.eq.	1074	)		coeff = 	0.935223782
        if (JELEM	.eq.	1075	)		coeff = 	0.87077543
        if (JELEM	.eq.	1076	)		coeff = 	1.089005153
        if (JELEM	.eq.	1077	)		coeff = 	0.988597965
        if (JELEM	.eq.	1078	)		coeff = 	0.925297644
        if (JELEM	.eq.	1079	)		coeff = 	0.876450699
        if (JELEM	.eq.	1080	)		coeff = 	0.95964975
        if (JELEM	.eq.	1081	)		coeff = 	1.118412291
        if (JELEM	.eq.	1082	)		coeff = 	0.98029407
        if (JELEM	.eq.	1083	)		coeff = 	0.869442346
        if (JELEM	.eq.	1084	)		coeff = 	0.623601363
        if (JELEM	.eq.	1085	)		coeff = 	0.834289105
        if (JELEM	.eq.	1086	)		coeff = 	1.108987976
        if (JELEM	.eq.	1087	)		coeff = 	0.875844875
        if (JELEM	.eq.	1088	)		coeff = 	0.892784815
        if (JELEM	.eq.	1089	)		coeff = 	0.880929795
        if (JELEM	.eq.	1090	)		coeff = 	1.04827336
        if (JELEM	.eq.	1091	)		coeff = 	1.01207348
        if (JELEM	.eq.	1092	)		coeff = 	0.935923232
        if (JELEM	.eq.	1093	)		coeff = 	0.686214036
        if (JELEM	.eq.	1094	)		coeff = 	0.977361456
        if (JELEM	.eq.	1095	)		coeff = 	0.965038789
        if (JELEM	.eq.	1096	)		coeff = 	0.701243382
        if (JELEM	.eq.	1097	)		coeff = 	0.877435737
        if (JELEM	.eq.	1098	)		coeff = 	0.677723813
        if (JELEM	.eq.	1099	)		coeff = 	1.047819382
        if (JELEM	.eq.	1100	)		coeff = 	1.131911515
        if (JELEM	.eq.	1101	)		coeff = 	1.063117732
        if (JELEM	.eq.	1102	)		coeff = 	0.894764021
        if (JELEM	.eq.	1103	)		coeff = 	0.871227772
        if (JELEM	.eq.	1104	)		coeff = 	0.920881894
        if (JELEM	.eq.	1105	)		coeff = 	1.008146579
        if (JELEM	.eq.	1106	)		coeff = 	0.91770019
        if (JELEM	.eq.	1107	)		coeff = 	1.119304404
        if (JELEM	.eq.	1108	)		coeff = 	0.822767106
        if (JELEM	.eq.	1109	)		coeff = 	1.085594269
        if (JELEM	.eq.	1110	)		coeff = 	0.747533291
        if (JELEM	.eq.	1111	)		coeff = 	0.996846279
        if (JELEM	.eq.	1112	)		coeff = 	1.072096471
        if (JELEM	.eq.	1113	)		coeff = 	1.015121962
        if (JELEM	.eq.	1114	)		coeff = 	0.799588914
        if (JELEM	.eq.	1115	)		coeff = 	0.869657987
        if (JELEM	.eq.	1116	)		coeff = 	1.051024056
        if (JELEM	.eq.	1117	)		coeff = 	1.057499028
        if (JELEM	.eq.	1118	)		coeff = 	0.803963904
        if (JELEM	.eq.	1119	)		coeff = 	0.874096725
        if (JELEM	.eq.	1120	)		coeff = 	0.864458623
        if (JELEM	.eq.	1121	)		coeff = 	1.044384554
        if (JELEM	.eq.	1122	)		coeff = 	1.050914419
        if (JELEM	.eq.	1123	)		coeff = 	0.920301123
        if (JELEM	.eq.	1124	)		coeff = 	1.081822925
        if (JELEM	.eq.	1125	)		coeff = 	0.962364672
        if (JELEM	.eq.	1126	)		coeff = 	0.907404563
        if (JELEM	.eq.	1127	)		coeff = 	1.063526891
        if (JELEM	.eq.	1128	)		coeff = 	0.998507311
        if (JELEM	.eq.	1129	)		coeff = 	0.856610557
        if (JELEM	.eq.	1130	)		coeff = 	0.79031094
        if (JELEM	.eq.	1131	)		coeff = 	1.062195109
        if (JELEM	.eq.	1132	)		coeff = 	0.887708766
        if (JELEM	.eq.	1133	)		coeff = 	0.948207403
        if (JELEM	.eq.	1134	)		coeff = 	1.052750876
        if (JELEM	.eq.	1135	)		coeff = 	0.905612504
        if (JELEM	.eq.	1136	)		coeff = 	1.018787846
        if (JELEM	.eq.	1137	)		coeff = 	0.924971648
        if (JELEM	.eq.	1138	)		coeff = 	1.037485447
        if (JELEM	.eq.	1139	)		coeff = 	0.966967451
        if (JELEM	.eq.	1140	)		coeff = 	0.879436701
        if (JELEM	.eq.	1141	)		coeff = 	0.928785963
        if (JELEM	.eq.	1142	)		coeff = 	1.033867007
        if (JELEM	.eq.	1143	)		coeff = 	0.918177538
        if (JELEM	.eq.	1144	)		coeff = 	1.010151848
        if (JELEM	.eq.	1145	)		coeff = 	1.078388476
        if (JELEM	.eq.	1146	)		coeff = 	0.916969566
        if (JELEM	.eq.	1147	)		coeff = 	0.747023968
        if (JELEM	.eq.	1148	)		coeff = 	0.839072754
        if (JELEM	.eq.	1149	)		coeff = 	1.067490225
        if (JELEM	.eq.	1150	)		coeff = 	0.906502125
        if (JELEM	.eq.	1151	)		coeff = 	0.91451008
        if (JELEM	.eq.	1152	)		coeff = 	0.961085124
        if (JELEM	.eq.	1153	)		coeff = 	0.98189951
        if (JELEM	.eq.	1154	)		coeff = 	0.995737724
        if (JELEM	.eq.	1155	)		coeff = 	0.838413528
        if (JELEM	.eq.	1156	)		coeff = 	1.032953942
        if (JELEM	.eq.	1157	)		coeff = 	1.021148066
        if (JELEM	.eq.	1158	)		coeff = 	0.946672014
        if (JELEM	.eq.	1159	)		coeff = 	0.926752539
        if (JELEM	.eq.	1160	)		coeff = 	0.977826528
        if (JELEM	.eq.	1161	)		coeff = 	0.82478942
        if (JELEM	.eq.	1162	)		coeff = 	0.736776534
        if (JELEM	.eq.	1163	)		coeff = 	0.852480155
        if (JELEM	.eq.	1164	)		coeff = 	0.894627998
        if (JELEM	.eq.	1165	)		coeff = 	1.077189258
        if (JELEM	.eq.	1166	)		coeff = 	1.095495568
        if (JELEM	.eq.	1167	)		coeff = 	0.840238079
        if (JELEM	.eq.	1168	)		coeff = 	1.004723187
        if (JELEM	.eq.	1169	)		coeff = 	0.79362385
        if (JELEM	.eq.	1170	)		coeff = 	0.941688604
        if (JELEM	.eq.	1171	)		coeff = 	1.115160862
        if (JELEM	.eq.	1172	)		coeff = 	1.043721888
        if (JELEM	.eq.	1173	)		coeff = 	1.067377764
        if (JELEM	.eq.	1174	)		coeff = 	0.986887291
        if (JELEM	.eq.	1175	)		coeff = 	0.815827374
        if (JELEM	.eq.	1176	)		coeff = 	1.090378627
        if (JELEM	.eq.	1177	)		coeff = 	0.76944184
        if (JELEM	.eq.	1178	)		coeff = 	0.762359815
        if (JELEM	.eq.	1179	)		coeff = 	1.028670952
        if (JELEM	.eq.	1180	)		coeff = 	1.127940248
        if (JELEM	.eq.	1181	)		coeff = 	0.973257661
        if (JELEM	.eq.	1182	)		coeff = 	0.939753409
        if (JELEM	.eq.	1183	)		coeff = 	0.892080258
        if (JELEM	.eq.	1184	)		coeff = 	1.009404759
        if (JELEM	.eq.	1185	)		coeff = 	1.106608533
        if (JELEM	.eq.	1186	)		coeff = 	0.805294409
        if (JELEM	.eq.	1187	)		coeff = 	0.845899393
        if (JELEM	.eq.	1188	)		coeff = 	0.998622699
        if (JELEM	.eq.	1189	)		coeff = 	0.876386544
        if (JELEM	.eq.	1190	)		coeff = 	0.981172936
        if (JELEM	.eq.	1191	)		coeff = 	0.958900249
        if (JELEM	.eq.	1192	)		coeff = 	0.946899602
        if (JELEM	.eq.	1193	)		coeff = 	1.034546808
        if (JELEM	.eq.	1194	)		coeff = 	1.074598291
        if (JELEM	.eq.	1195	)		coeff = 	0.931813717
        if (JELEM	.eq.	1196	)		coeff = 	0.930857741
        if (JELEM	.eq.	1197	)		coeff = 	1.037745155
        if (JELEM	.eq.	1198	)		coeff = 	0.829164664
        if (JELEM	.eq.	1199	)		coeff = 	1.080062631
        if (JELEM	.eq.	1200	)		coeff = 	0.812239058
        if (JELEM	.eq.	1201	)		coeff = 	0.994296796
        if (JELEM	.eq.	1202	)		coeff = 	1.018786118
        if (JELEM	.eq.	1203	)		coeff = 	0.881661472
        if (JELEM	.eq.	1204	)		coeff = 	1.028726019
        if (JELEM	.eq.	1205	)		coeff = 	0.956846432
        if (JELEM	.eq.	1206	)		coeff = 	1.039490525
        if (JELEM	.eq.	1207	)		coeff = 	1.111404054
        if (JELEM	.eq.	1208	)		coeff = 	0.846540411
        if (JELEM	.eq.	1209	)		coeff = 	1.07984201
        if (JELEM	.eq.	1210	)		coeff = 	0.9521819
        if (JELEM	.eq.	1211	)		coeff = 	0.927378091
        if (JELEM	.eq.	1212	)		coeff = 	1.229384829
        if (JELEM	.eq.	1213	)		coeff = 	1.061142181
        if (JELEM	.eq.	1214	)		coeff = 	1.166243856
        if (JELEM	.eq.	1215	)		coeff = 	0.881625442
        if (JELEM	.eq.	1216	)		coeff = 	0.900697281
        if (JELEM	.eq.	1217	)		coeff = 	1.003141901
        if (JELEM	.eq.	1218	)		coeff = 	0.912412203
        if (JELEM	.eq.	1219	)		coeff = 	1.091354641
        if (JELEM	.eq.	1220	)		coeff = 	0.849282293
        if (JELEM	.eq.	1221	)		coeff = 	0.985015209
        if (JELEM	.eq.	1222	)		coeff = 	1.086233894
        if (JELEM	.eq.	1223	)		coeff = 	1.088126514
        if (JELEM	.eq.	1224	)		coeff = 	1.084034821
        if (JELEM	.eq.	1225	)		coeff = 	0.929258368
        if (JELEM	.eq.	1226	)		coeff = 	1.04841767
        if (JELEM	.eq.	1227	)		coeff = 	0.985021147
        if (JELEM	.eq.	1228	)		coeff = 	0.856286377
        if (JELEM	.eq.	1229	)		coeff = 	1.010668156
        if (JELEM	.eq.	1230	)		coeff = 	1.02634671
        if (JELEM	.eq.	1231	)		coeff = 	0.85128893
        if (JELEM	.eq.	1232	)		coeff = 	1.04664686
        if (JELEM	.eq.	1233	)		coeff = 	0.917728473
        if (JELEM	.eq.	1234	)		coeff = 	0.942246184
        if (JELEM	.eq.	1235	)		coeff = 	0.945869908
        if (JELEM	.eq.	1236	)		coeff = 	1.01462033
        if (JELEM	.eq.	1237	)		coeff = 	1.067765445
        if (JELEM	.eq.	1238	)		coeff = 	1.073821822
        if (JELEM	.eq.	1239	)		coeff = 	0.838197299
        if (JELEM	.eq.	1240	)		coeff = 	0.971025606
        if (JELEM	.eq.	1241	)		coeff = 	1.038733982
        if (JELEM	.eq.	1242	)		coeff = 	0.869544664
        if (JELEM	.eq.	1243	)		coeff = 	0.908276569
        if (JELEM	.eq.	1244	)		coeff = 	1.062003996
        if (JELEM	.eq.	1245	)		coeff = 	1.069480558
        if (JELEM	.eq.	1246	)		coeff = 	0.912979479
        if (JELEM	.eq.	1247	)		coeff = 	0.916345203
        if (JELEM	.eq.	1248	)		coeff = 	0.921256538
        if (JELEM	.eq.	1249	)		coeff = 	1.053109504
        if (JELEM	.eq.	1250	)		coeff = 	0.84957624
        if (JELEM	.eq.	1251	)		coeff = 	0.930802335
        if (JELEM	.eq.	1252	)		coeff = 	0.766697419
        if (JELEM	.eq.	1253	)		coeff = 	0.915448907
        if (JELEM	.eq.	1254	)		coeff = 	1.033440417
        if (JELEM	.eq.	1255	)		coeff = 	1.118488185
        if (JELEM	.eq.	1256	)		coeff = 	1.126574154
        if (JELEM	.eq.	1257	)		coeff = 	0.895305171
        if (JELEM	.eq.	1258	)		coeff = 	1.035174201
        if (JELEM	.eq.	1259	)		coeff = 	1.067760691
        if (JELEM	.eq.	1260	)		coeff = 	1.120673856
        if (JELEM	.eq.	1261	)		coeff = 	0.856382079
        if (JELEM	.eq.	1262	)		coeff = 	0.839320796
        if (JELEM	.eq.	1263	)		coeff = 	1.039403637
        if (JELEM	.eq.	1264	)		coeff = 	0.88374039
        if (JELEM	.eq.	1265	)		coeff = 	1.02386953
        if (JELEM	.eq.	1266	)		coeff = 	1.008095201
        if (JELEM	.eq.	1267	)		coeff = 	1.010588107
        if (JELEM	.eq.	1268	)		coeff = 	1.10477397
        if (JELEM	.eq.	1269	)		coeff = 	1.098644577
        if (JELEM	.eq.	1270	)		coeff = 	0.965447669
        if (JELEM	.eq.	1271	)		coeff = 	0.997768634
        if (JELEM	.eq.	1272	)		coeff = 	0.986810111
        if (JELEM	.eq.	1273	)		coeff = 	0.928949619
        if (JELEM	.eq.	1274	)		coeff = 	0.91473996
        if (JELEM	.eq.	1275	)		coeff = 	0.923934822
        if (JELEM	.eq.	1276	)		coeff = 	1.007559283
        if (JELEM	.eq.	1277	)		coeff = 	0.884989181
        if (JELEM	.eq.	1278	)		coeff = 	0.903658971
        if (JELEM	.eq.	1279	)		coeff = 	0.939094883
        if (JELEM	.eq.	1280	)		coeff = 	0.812071571
        if (JELEM	.eq.	1281	)		coeff = 	0.731275873
        if (JELEM	.eq.	1282	)		coeff = 	0.827599485
        if (JELEM	.eq.	1283	)		coeff = 	1.013038085
        if (JELEM	.eq.	1284	)		coeff = 	0.914720313
        if (JELEM	.eq.	1285	)		coeff = 	0.793627661
        if (JELEM	.eq.	1286	)		coeff = 	0.907233074
        if (JELEM	.eq.	1287	)		coeff = 	0.494574597
        if (JELEM	.eq.	1288	)		coeff = 	0.964180758
        if (JELEM	.eq.	1289	)		coeff = 	1.026477045
        if (JELEM	.eq.	1290	)		coeff = 	0.904407027
        if (JELEM	.eq.	1291	)		coeff = 	0.588703435
        if (JELEM	.eq.	1292	)		coeff = 	1.045395117
        if (JELEM	.eq.	1293	)		coeff = 	0.987041643
        if (JELEM	.eq.	1294	)		coeff = 	1.061719578
        if (JELEM	.eq.	1295	)		coeff = 	1.033285632
        if (JELEM	.eq.	1296	)		coeff = 	0.992998107
        if (JELEM	.eq.	1297	)		coeff = 	1.007568891
        if (JELEM	.eq.	1298	)		coeff = 	0.77093162
        if (JELEM	.eq.	1299	)		coeff = 	0.977008729
        if (JELEM	.eq.	1300	)		coeff = 	0.956558141
        if (JELEM	.eq.	1301	)		coeff = 	0.963691099
        if (JELEM	.eq.	1302	)		coeff = 	0.855633303
        if (JELEM	.eq.	1303	)		coeff = 	0.858358646
        if (JELEM	.eq.	1304	)		coeff = 	1.078978415
        if (JELEM	.eq.	1305	)		coeff = 	0.929367404
        if (JELEM	.eq.	1306	)		coeff = 	0.873239034
        if (JELEM	.eq.	1307	)		coeff = 	1.02304136
        if (JELEM	.eq.	1308	)		coeff = 	1.027574081
        if (JELEM	.eq.	1309	)		coeff = 	0.973353682
        if (JELEM	.eq.	1310	)		coeff = 	0.999084106
        if (JELEM	.eq.	1311	)		coeff = 	0.96975351
        if (JELEM	.eq.	1312	)		coeff = 	0.980954073
        if (JELEM	.eq.	1313	)		coeff = 	1.031650309
        if (JELEM	.eq.	1314	)		coeff = 	1.069502256
        if (JELEM	.eq.	1315	)		coeff = 	1.028433212
        if (JELEM	.eq.	1316	)		coeff = 	1.005288717
        if (JELEM	.eq.	1317	)		coeff = 	1.032452203
        if (JELEM	.eq.	1318	)		coeff = 	1.149647451
        if (JELEM	.eq.	1319	)		coeff = 	0.977546777
        if (JELEM	.eq.	1320	)		coeff = 	0.794552472
        if (JELEM	.eq.	1321	)		coeff = 	0.916900206
        if (JELEM	.eq.	1322	)		coeff = 	1.06861803
        if (JELEM	.eq.	1323	)		coeff = 	1.030684311
        if (JELEM	.eq.	1324	)		coeff = 	0.864857097
        if (JELEM	.eq.	1325	)		coeff = 	0.891835049
        if (JELEM	.eq.	1326	)		coeff = 	0.849036914
        if (JELEM	.eq.	1327	)		coeff = 	0.761786723
        if (JELEM	.eq.	1328	)		coeff = 	1.051644315
        if (JELEM	.eq.	1329	)		coeff = 	0.956451323
        if (JELEM	.eq.	1330	)		coeff = 	1.024253285
        if (JELEM	.eq.	1331	)		coeff = 	0.989027687
        if (JELEM	.eq.	1332	)		coeff = 	1.036250114
        if (JELEM	.eq.	1333	)		coeff = 	1.002810629
        if (JELEM	.eq.	1334	)		coeff = 	1.06088956
        if (JELEM	.eq.	1335	)		coeff = 	1.195824991
        if (JELEM	.eq.	1336	)		coeff = 	0.989495361
        if (JELEM	.eq.	1337	)		coeff = 	1.02453845
        if (JELEM	.eq.	1338	)		coeff = 	1.056345214
        if (JELEM	.eq.	1339	)		coeff = 	0.939370562
        if (JELEM	.eq.	1340	)		coeff = 	0.894076557
        if (JELEM	.eq.	1341	)		coeff = 	0.8941742
        if (JELEM	.eq.	1342	)		coeff = 	0.961617646
        if (JELEM	.eq.	1343	)		coeff = 	1.011495305
        if (JELEM	.eq.	1344	)		coeff = 	0.976432846
        if (JELEM	.eq.	1345	)		coeff = 	1.019617026
        if (JELEM	.eq.	1346	)		coeff = 	0.908066582
        if (JELEM	.eq.	1347	)		coeff = 	0.77508883
        if (JELEM	.eq.	1348	)		coeff = 	1.021072525
        if (JELEM	.eq.	1349	)		coeff = 	0.970645944
        if (JELEM	.eq.	1350	)		coeff = 	0.993583515
        if (JELEM	.eq.	1351	)		coeff = 	1.013445658
        if (JELEM	.eq.	1352	)		coeff = 	0.991245685
        if (JELEM	.eq.	1353	)		coeff = 	0.815802995
        if (JELEM	.eq.	1354	)		coeff = 	0.738072248
        if (JELEM	.eq.	1355	)		coeff = 	1.184348074
        if (JELEM	.eq.	1356	)		coeff = 	0.934597625
        if (JELEM	.eq.	1357	)		coeff = 	0.867485479
        if (JELEM	.eq.	1358	)		coeff = 	0.922971254
        if (JELEM	.eq.	1359	)		coeff = 	1.079324262
        if (JELEM	.eq.	1360	)		coeff = 	0.89820028
        if (JELEM	.eq.	1361	)		coeff = 	1.088049704
        if (JELEM	.eq.	1362	)		coeff = 	0.939637951
        if (JELEM	.eq.	1363	)		coeff = 	0.835315238
        if (JELEM	.eq.	1364	)		coeff = 	1.03903221
        if (JELEM	.eq.	1365	)		coeff = 	1.050052051
        if (JELEM	.eq.	1366	)		coeff = 	1.125944989
        if (JELEM	.eq.	1367	)		coeff = 	0.846667144
        if (JELEM	.eq.	1368	)		coeff = 	1.135588962
        if (JELEM	.eq.	1369	)		coeff = 	1.000174211
        if (JELEM	.eq.	1370	)		coeff = 	1.057196761
        if (JELEM	.eq.	1371	)		coeff = 	0.962877909
        if (JELEM	.eq.	1372	)		coeff = 	0.965127378
        if (JELEM	.eq.	1373	)		coeff = 	0.824490298
        if (JELEM	.eq.	1374	)		coeff = 	0.943687572
        if (JELEM	.eq.	1375	)		coeff = 	1.058520041
        if (JELEM	.eq.	1376	)		coeff = 	0.791192886
        if (JELEM	.eq.	1377	)		coeff = 	0.914703592
        if (JELEM	.eq.	1378	)		coeff = 	0.963263672
        if (JELEM	.eq.	1379	)		coeff = 	0.849583359
        if (JELEM	.eq.	1380	)		coeff = 	1.022225019
        if (JELEM	.eq.	1381	)		coeff = 	0.81337428
        if (JELEM	.eq.	1382	)		coeff = 	0.95517574
        if (JELEM	.eq.	1383	)		coeff = 	1.090559099
        if (JELEM	.eq.	1384	)		coeff = 	0.968568368
        if (JELEM	.eq.	1385	)		coeff = 	1.012639686
        if (JELEM	.eq.	1386	)		coeff = 	0.951849718
        if (JELEM	.eq.	1387	)		coeff = 	0.73215407
        if (JELEM	.eq.	1388	)		coeff = 	1.000096575
        if (JELEM	.eq.	1389	)		coeff = 	1.137668473
        if (JELEM	.eq.	1390	)		coeff = 	0.887088886
        if (JELEM	.eq.	1391	)		coeff = 	1.011783533
        if (JELEM	.eq.	1392	)		coeff = 	0.914148503
        if (JELEM	.eq.	1393	)		coeff = 	0.964610918
        if (JELEM	.eq.	1394	)		coeff = 	0.829117307
        if (JELEM	.eq.	1395	)		coeff = 	0.883417554
        if (JELEM	.eq.	1396	)		coeff = 	0.99185967
        if (JELEM	.eq.	1397	)		coeff = 	0.960068844
        if (JELEM	.eq.	1398	)		coeff = 	0.939940037
        if (JELEM	.eq.	1399	)		coeff = 	0.961713369
        if (JELEM	.eq.	1400	)		coeff = 	0.889314058
        if (JELEM	.eq.	1401	)		coeff = 	0.891368187
        if (JELEM	.eq.	1402	)		coeff = 	1.050664243
        if (JELEM	.eq.	1403	)		coeff = 	0.976417917
        if (JELEM	.eq.	1404	)		coeff = 	0.958803234
        if (JELEM	.eq.	1405	)		coeff = 	1.066471389
        if (JELEM	.eq.	1406	)		coeff = 	0.869139545
        if (JELEM	.eq.	1407	)		coeff = 	1.060492696
        if (JELEM	.eq.	1408	)		coeff = 	0.965338855
        if (JELEM	.eq.	1409	)		coeff = 	1.030484519
        if (JELEM	.eq.	1410	)		coeff = 	0.928299102
        if (JELEM	.eq.	1411	)		coeff = 	0.975812997
        if (JELEM	.eq.	1412	)		coeff = 	1.031582361
        if (JELEM	.eq.	1413	)		coeff = 	0.954332462
        if (JELEM	.eq.	1414	)		coeff = 	1.091041861
        if (JELEM	.eq.	1415	)		coeff = 	0.946507285
        if (JELEM	.eq.	1416	)		coeff = 	0.832518953
        if (JELEM	.eq.	1417	)		coeff = 	0.891792398
        if (JELEM	.eq.	1418	)		coeff = 	0.9531057
        if (JELEM	.eq.	1419	)		coeff = 	0.91200162
        if (JELEM	.eq.	1420	)		coeff = 	0.827213644
        if (JELEM	.eq.	1421	)		coeff = 	0.949299745
        if (JELEM	.eq.	1422	)		coeff = 	0.789058411
        if (JELEM	.eq.	1423	)		coeff = 	0.947901177
        if (JELEM	.eq.	1424	)		coeff = 	0.828705574
        if (JELEM	.eq.	1425	)		coeff = 	0.976998121
        if (JELEM	.eq.	1426	)		coeff = 	1.062221734
        if (JELEM	.eq.	1427	)		coeff = 	1.037049962
        if (JELEM	.eq.	1428	)		coeff = 	0.931423654
        if (JELEM	.eq.	1429	)		coeff = 	0.98217067
        if (JELEM	.eq.	1430	)		coeff = 	0.932840729
        if (JELEM	.eq.	1431	)		coeff = 	1.037219969
        if (JELEM	.eq.	1432	)		coeff = 	0.765157817
        if (JELEM	.eq.	1433	)		coeff = 	0.90080681
        if (JELEM	.eq.	1434	)		coeff = 	0.910625685
        if (JELEM	.eq.	1435	)		coeff = 	0.993017553
        if (JELEM	.eq.	1436	)		coeff = 	0.920051296
        if (JELEM	.eq.	1437	)		coeff = 	0.97795285
        if (JELEM	.eq.	1438	)		coeff = 	0.631452937
        if (JELEM	.eq.	1439	)		coeff = 	1.090055376
        if (JELEM	.eq.	1440	)		coeff = 	0.972735743
        if (JELEM	.eq.	1441	)		coeff = 	0.904310141
        if (JELEM	.eq.	1442	)		coeff = 	0.754699149
        if (JELEM	.eq.	1443	)		coeff = 	0.925212891
        if (JELEM	.eq.	1444	)		coeff = 	0.929752973
        if (JELEM	.eq.	1445	)		coeff = 	0.921346541
        if (JELEM	.eq.	1446	)		coeff = 	1.065555615
        if (JELEM	.eq.	1447	)		coeff = 	1.034798297
        if (JELEM	.eq.	1448	)		coeff = 	0.949816006
        if (JELEM	.eq.	1449	)		coeff = 	1.035810404
        if (JELEM	.eq.	1450	)		coeff = 	0.962880465
        if (JELEM	.eq.	1451	)		coeff = 	0.860108127
        if (JELEM	.eq.	1452	)		coeff = 	0.91474483
        if (JELEM	.eq.	1453	)		coeff = 	1.126913584
        if (JELEM	.eq.	1454	)		coeff = 	1.046145315
        if (JELEM	.eq.	1455	)		coeff = 	0.990451293
        if (JELEM	.eq.	1456	)		coeff = 	0.781797175
        if (JELEM	.eq.	1457	)		coeff = 	0.932583146
        if (JELEM	.eq.	1458	)		coeff = 	1.042889242
        if (JELEM	.eq.	1459	)		coeff = 	1.00517057
        if (JELEM	.eq.	1460	)		coeff = 	1.02792365
        if (JELEM	.eq.	1461	)		coeff = 	0.995250781
        if (JELEM	.eq.	1462	)		coeff = 	0.81531593
        if (JELEM	.eq.	1463	)		coeff = 	0.948402173
        if (JELEM	.eq.	1464	)		coeff = 	1.011777135
        if (JELEM	.eq.	1465	)		coeff = 	0.919713618
        if (JELEM	.eq.	1466	)		coeff = 	0.931838357
        if (JELEM	.eq.	1467	)		coeff = 	1.092378794
        if (JELEM	.eq.	1468	)		coeff = 	1.000291974
        if (JELEM	.eq.	1469	)		coeff = 	0.757159393
        if (JELEM	.eq.	1470	)		coeff = 	0.777117618
        if (JELEM	.eq.	1471	)		coeff = 	1.029989686
        if (JELEM	.eq.	1472	)		coeff = 	1.028900334
        if (JELEM	.eq.	1473	)		coeff = 	0.869864702
        if (JELEM	.eq.	1474	)		coeff = 	0.990779974
        if (JELEM	.eq.	1475	)		coeff = 	0.809295846
        if (JELEM	.eq.	1476	)		coeff = 	0.977154896
        if (JELEM	.eq.	1477	)		coeff = 	0.964849148
        if (JELEM	.eq.	1478	)		coeff = 	0.998619364
        if (JELEM	.eq.	1479	)		coeff = 	0.996246772
        if (JELEM	.eq.	1480	)		coeff = 	0.975547382
        if (JELEM	.eq.	1481	)		coeff = 	0.895487491
        if (JELEM	.eq.	1482	)		coeff = 	0.890841352
        if (JELEM	.eq.	1483	)		coeff = 	1.008678248
        if (JELEM	.eq.	1484	)		coeff = 	0.973850945
        if (JELEM	.eq.	1485	)		coeff = 	0.997332192
        if (JELEM	.eq.	1486	)		coeff = 	1.106439939
        if (JELEM	.eq.	1487	)		coeff = 	0.66708022
        if (JELEM	.eq.	1488	)		coeff = 	0.936161919
        if (JELEM	.eq.	1489	)		coeff = 	0.97771296
        if (JELEM	.eq.	1490	)		coeff = 	0.700348767
        if (JELEM	.eq.	1491	)		coeff = 	0.906092876
        if (JELEM	.eq.	1492	)		coeff = 	0.967343269
        if (JELEM	.eq.	1493	)		coeff = 	1.002867659
        if (JELEM	.eq.	1494	)		coeff = 	0.982776093
        if (JELEM	.eq.	1495	)		coeff = 	1.046734165
        if (JELEM	.eq.	1496	)		coeff = 	0.951341383
        if (JELEM	.eq.	1497	)		coeff = 	0.730762681
        if (JELEM	.eq.	1498	)		coeff = 	1.053036534
        if (JELEM	.eq.	1499	)		coeff = 	0.874340293
        if (JELEM	.eq.	1500	)		coeff = 	0.841228809
        if (JELEM	.eq.	1501	)		coeff = 	0.791228457
        if (JELEM	.eq.	1502	)		coeff = 	0.925483201
        if (JELEM	.eq.	1503	)		coeff = 	0.93698004
        if (JELEM	.eq.	1504	)		coeff = 	0.873293757
        if (JELEM	.eq.	1505	)		coeff = 	1.026428121
        if (JELEM	.eq.	1506	)		coeff = 	1.169810397
        if (JELEM	.eq.	1507	)		coeff = 	1.027778208
        if (JELEM	.eq.	1508	)		coeff = 	0.95448803
        if (JELEM	.eq.	1509	)		coeff = 	1.040755165
        if (JELEM	.eq.	1510	)		coeff = 	0.787369634
        if (JELEM	.eq.	1511	)		coeff = 	1.081182693
        if (JELEM	.eq.	1512	)		coeff = 	1.028383567
        if (JELEM	.eq.	1513	)		coeff = 	0.969600698
        if (JELEM	.eq.	1514	)		coeff = 	0.915610169
        if (JELEM	.eq.	1515	)		coeff = 	0.84203737
        if (JELEM	.eq.	1516	)		coeff = 	0.924190069
        if (JELEM	.eq.	1517	)		coeff = 	0.857371152
        if (JELEM	.eq.	1518	)		coeff = 	0.943309508
        if (JELEM	.eq.	1519	)		coeff = 	0.96631225
        if (JELEM	.eq.	1520	)		coeff = 	0.859981219
        if (JELEM	.eq.	1521	)		coeff = 	0.980182283
        if (JELEM	.eq.	1522	)		coeff = 	1.116813142
        if (JELEM	.eq.	1523	)		coeff = 	0.645806595
        if (JELEM	.eq.	1524	)		coeff = 	0.883165634
        if (JELEM	.eq.	1525	)		coeff = 	0.941464639
        if (JELEM	.eq.	1526	)		coeff = 	1.00921951
        if (JELEM	.eq.	1527	)		coeff = 	0.937053408
        if (JELEM	.eq.	1528	)		coeff = 	0.79709264
        if (JELEM	.eq.	1529	)		coeff = 	0.838035397
        if (JELEM	.eq.	1530	)		coeff = 	0.822143935
        if (JELEM	.eq.	1531	)		coeff = 	1.085281625
        if (JELEM	.eq.	1532	)		coeff = 	0.866690125
        if (JELEM	.eq.	1533	)		coeff = 	0.989683884
        if (JELEM	.eq.	1534	)		coeff = 	1.016776642
        if (JELEM	.eq.	1535	)		coeff = 	1.137597934
        if (JELEM	.eq.	1536	)		coeff = 	0.898221321
        if (JELEM	.eq.	1537	)		coeff = 	0.837637611
        if (JELEM	.eq.	1538	)		coeff = 	1.11604795
        if (JELEM	.eq.	1539	)		coeff = 	1.034509149
        if (JELEM	.eq.	1540	)		coeff = 	0.742790063
        if (JELEM	.eq.	1541	)		coeff = 	1.082358837
        if (JELEM	.eq.	1542	)		coeff = 	0.729927524
        if (JELEM	.eq.	1543	)		coeff = 	1.00613697
        if (JELEM	.eq.	1544	)		coeff = 	0.813750168
        if (JELEM	.eq.	1545	)		coeff = 	1.02557212
        if (JELEM	.eq.	1546	)		coeff = 	0.896737837
        if (JELEM	.eq.	1547	)		coeff = 	1.116176498
        if (JELEM	.eq.	1548	)		coeff = 	1.083628215
        if (JELEM	.eq.	1549	)		coeff = 	0.828292133
        if (JELEM	.eq.	1550	)		coeff = 	0.984728414
        if (JELEM	.eq.	1551	)		coeff = 	0.946217188
        if (JELEM	.eq.	1552	)		coeff = 	0.889762291
        if (JELEM	.eq.	1553	)		coeff = 	0.724603066
        if (JELEM	.eq.	1554	)		coeff = 	1.10955007
        if (JELEM	.eq.	1555	)		coeff = 	1.078234244
        if (JELEM	.eq.	1556	)		coeff = 	1.073586157
        if (JELEM	.eq.	1557	)		coeff = 	0.859394612
        if (JELEM	.eq.	1558	)		coeff = 	1.062382794
        if (JELEM	.eq.	1559	)		coeff = 	0.956851546
        if (JELEM	.eq.	1560	)		coeff = 	0.898909129
        if (JELEM	.eq.	1561	)		coeff = 	0.628714811
        if (JELEM	.eq.	1562	)		coeff = 	0.931828599
        if (JELEM	.eq.	1563	)		coeff = 	1.02111894
        if (JELEM	.eq.	1564	)		coeff = 	0.94960275
        if (JELEM	.eq.	1565	)		coeff = 	0.952529441
        if (JELEM	.eq.	1566	)		coeff = 	0.798051209
        if (JELEM	.eq.	1567	)		coeff = 	1.102199759
        if (JELEM	.eq.	1568	)		coeff = 	0.662882091
        if (JELEM	.eq.	1569	)		coeff = 	0.7935734
        if (JELEM	.eq.	1570	)		coeff = 	0.758399378
        if (JELEM	.eq.	1571	)		coeff = 	0.911918597
        if (JELEM	.eq.	1572	)		coeff = 	0.892153031
        if (JELEM	.eq.	1573	)		coeff = 	1.090175911
        if (JELEM	.eq.	1574	)		coeff = 	1.071718787
        if (JELEM	.eq.	1575	)		coeff = 	0.8554742
        if (JELEM	.eq.	1576	)		coeff = 	0.836661243
        if (JELEM	.eq.	1577	)		coeff = 	1.08558958
        if (JELEM	.eq.	1578	)		coeff = 	0.880425525
        if (JELEM	.eq.	1579	)		coeff = 	1.012812737
        if (JELEM	.eq.	1580	)		coeff = 	1.036739533
        if (JELEM	.eq.	1581	)		coeff = 	0.839965677
        if (JELEM	.eq.	1582	)		coeff = 	0.985495849
        if (JELEM	.eq.	1583	)		coeff = 	0.915936247
        if (JELEM	.eq.	1584	)		coeff = 	1.127845661
        if (JELEM	.eq.	1585	)		coeff = 	0.886574197
        if (JELEM	.eq.	1586	)		coeff = 	0.547866168
        if (JELEM	.eq.	1587	)		coeff = 	0.995586287
        if (JELEM	.eq.	1588	)		coeff = 	0.955795434
        if (JELEM	.eq.	1589	)		coeff = 	0.868130603
        if (JELEM	.eq.	1590	)		coeff = 	1.061650001
        if (JELEM	.eq.	1591	)		coeff = 	0.797284384
        if (JELEM	.eq.	1592	)		coeff = 	1.015846347
        if (JELEM	.eq.	1593	)		coeff = 	0.909358282
        if (JELEM	.eq.	1594	)		coeff = 	0.975247339
        if (JELEM	.eq.	1595	)		coeff = 	0.914519837
        if (JELEM	.eq.	1596	)		coeff = 	0.908111402
        if (JELEM	.eq.	1597	)		coeff = 	1.035502352
        if (JELEM	.eq.	1598	)		coeff = 	0.951663251
        if (JELEM	.eq.	1599	)		coeff = 	1.060922225
        if (JELEM	.eq.	1600	)		coeff = 	0.79277883
        if (JELEM	.eq.	1601	)		coeff = 	0.943220379
        if (JELEM	.eq.	1602	)		coeff = 	0.836668196
        if (JELEM	.eq.	1603	)		coeff = 	0.690382469
        if (JELEM	.eq.	1604	)		coeff = 	0.862452554
        if (JELEM	.eq.	1605	)		coeff = 	0.948663491
        if (JELEM	.eq.	1606	)		coeff = 	0.818949908
        if (JELEM	.eq.	1607	)		coeff = 	1.028030395
        if (JELEM	.eq.	1608	)		coeff = 	0.755700427
        if (JELEM	.eq.	1609	)		coeff = 	0.629583047
        if (JELEM	.eq.	1610	)		coeff = 	0.78265953
        if (JELEM	.eq.	1611	)		coeff = 	1.126517507
        if (JELEM	.eq.	1612	)		coeff = 	0.939273911
        if (JELEM	.eq.	1613	)		coeff = 	0.985614142
        if (JELEM	.eq.	1614	)		coeff = 	1.119405272
        if (JELEM	.eq.	1615	)		coeff = 	0.774566662
        if (JELEM	.eq.	1616	)		coeff = 	0.882307014
        if (JELEM	.eq.	1617	)		coeff = 	1.108717708
        if (JELEM	.eq.	1618	)		coeff = 	0.893842577
        if (JELEM	.eq.	1619	)		coeff = 	1.022614524
        if (JELEM	.eq.	1620	)		coeff = 	0.717527036
        if (JELEM	.eq.	1621	)		coeff = 	1.066379299
        if (JELEM	.eq.	1622	)		coeff = 	0.935885362
        if (JELEM	.eq.	1623	)		coeff = 	0.973612858
        if (JELEM	.eq.	1624	)		coeff = 	0.988230186
        if (JELEM	.eq.	1625	)		coeff = 	0.813720572
        if (JELEM	.eq.	1626	)		coeff = 	1.040833603
        if (JELEM	.eq.	1627	)		coeff = 	0.951883688
        if (JELEM	.eq.	1628	)		coeff = 	0.74661126
        if (JELEM	.eq.	1629	)		coeff = 	0.918133717
        if (JELEM	.eq.	1630	)		coeff = 	1.12946545
        if (JELEM	.eq.	1631	)		coeff = 	0.994311283
        if (JELEM	.eq.	1632	)		coeff = 	0.818386372
        if (JELEM	.eq.	1633	)		coeff = 	0.994458445
        if (JELEM	.eq.	1634	)		coeff = 	0.884118987
        if (JELEM	.eq.	1635	)		coeff = 	0.886453301
        if (JELEM	.eq.	1636	)		coeff = 	1.03932337
        if (JELEM	.eq.	1637	)		coeff = 	0.776075592
        if (JELEM	.eq.	1638	)		coeff = 	1.05816133
        if (JELEM	.eq.	1639	)		coeff = 	0.845195478
        if (JELEM	.eq.	1640	)		coeff = 	0.988363208
        if (JELEM	.eq.	1641	)		coeff = 	0.874282872
        if (JELEM	.eq.	1642	)		coeff = 	1.007608661
        if (JELEM	.eq.	1643	)		coeff = 	0.805532143
        if (JELEM	.eq.	1644	)		coeff = 	0.830530606
        if (JELEM	.eq.	1645	)		coeff = 	1.073142154
        if (JELEM	.eq.	1646	)		coeff = 	0.908296838
        if (JELEM	.eq.	1647	)		coeff = 	1.066285966
        if (JELEM	.eq.	1648	)		coeff = 	0.988492087
        if (JELEM	.eq.	1649	)		coeff = 	0.949872225
        if (JELEM	.eq.	1650	)		coeff = 	0.949335933
        if (JELEM	.eq.	1651	)		coeff = 	0.873120364
        if (JELEM	.eq.	1652	)		coeff = 	0.886595327
        if (JELEM	.eq.	1653	)		coeff = 	0.791805019
        if (JELEM	.eq.	1654	)		coeff = 	0.974010219
        if (JELEM	.eq.	1655	)		coeff = 	0.991078723
        if (JELEM	.eq.	1656	)		coeff = 	0.751870563
        if (JELEM	.eq.	1657	)		coeff = 	1.089078601
        if (JELEM	.eq.	1658	)		coeff = 	0.947319744
        if (JELEM	.eq.	1659	)		coeff = 	1.013863124
        if (JELEM	.eq.	1660	)		coeff = 	0.750662692
        if (JELEM	.eq.	1661	)		coeff = 	0.929194831
        if (JELEM	.eq.	1662	)		coeff = 	0.973826914
        if (JELEM	.eq.	1663	)		coeff = 	0.814773539
        if (JELEM	.eq.	1664	)		coeff = 	0.88977767
        if (JELEM	.eq.	1665	)		coeff = 	0.938301829
        if (JELEM	.eq.	1666	)		coeff = 	1.086676452
        if (JELEM	.eq.	1667	)		coeff = 	1.143489109
        if (JELEM	.eq.	1668	)		coeff = 	1.094786917
        if (JELEM	.eq.	1669	)		coeff = 	0.976347039
        if (JELEM	.eq.	1670	)		coeff = 	1.05259611
        if (JELEM	.eq.	1671	)		coeff = 	1.099781348
        if (JELEM	.eq.	1672	)		coeff = 	1.010626642
        if (JELEM	.eq.	1673	)		coeff = 	1.040091128
        if (JELEM	.eq.	1674	)		coeff = 	1.008392357
        if (JELEM	.eq.	1675	)		coeff = 	0.949149322
        if (JELEM	.eq.	1676	)		coeff = 	0.709231072
        if (JELEM	.eq.	1677	)		coeff = 	1.040108826
        if (JELEM	.eq.	1678	)		coeff = 	0.991256547
        if (JELEM	.eq.	1679	)		coeff = 	0.626202075
        if (JELEM	.eq.	1680	)		coeff = 	0.974473844
        if (JELEM	.eq.	1681	)		coeff = 	0.659607924
        if (JELEM	.eq.	1682	)		coeff = 	0.964445868
        if (JELEM	.eq.	1683	)		coeff = 	0.642933087
        if (JELEM	.eq.	1684	)		coeff = 	0.953276723
        if (JELEM	.eq.	1685	)		coeff = 	1.072743992
        if (JELEM	.eq.	1686	)		coeff = 	1.007833602
        if (JELEM	.eq.	1687	)		coeff = 	0.891206503
        if (JELEM	.eq.	1688	)		coeff = 	0.987065138
        if (JELEM	.eq.	1689	)		coeff = 	0.775096018
        if (JELEM	.eq.	1690	)		coeff = 	1.102574962
        if (JELEM	.eq.	1691	)		coeff = 	0.996873982
        if (JELEM	.eq.	1692	)		coeff = 	0.945446244
        if (JELEM	.eq.	1693	)		coeff = 	0.880135917
        if (JELEM	.eq.	1694	)		coeff = 	0.889394419
        if (JELEM	.eq.	1695	)		coeff = 	0.996372695
        if (JELEM	.eq.	1696	)		coeff = 	0.794487005
        if (JELEM	.eq.	1697	)		coeff = 	1.051030897
        if (JELEM	.eq.	1698	)		coeff = 	0.806304753
        if (JELEM	.eq.	1699	)		coeff = 	1.030778978
        if (JELEM	.eq.	1700	)		coeff = 	0.987602522
        if (JELEM	.eq.	1701	)		coeff = 	1.190530311
        if (JELEM	.eq.	1702	)		coeff = 	0.951057624
        if (JELEM	.eq.	1703	)		coeff = 	1.006240822
        if (JELEM	.eq.	1704	)		coeff = 	0.920676605
        if (JELEM	.eq.	1705	)		coeff = 	0.76674019
        if (JELEM	.eq.	1706	)		coeff = 	0.938755154
        if (JELEM	.eq.	1707	)		coeff = 	1.030064413
        if (JELEM	.eq.	1708	)		coeff = 	0.798532473
        if (JELEM	.eq.	1709	)		coeff = 	1.150414454
        if (JELEM	.eq.	1710	)		coeff = 	1.03989536
        if (JELEM	.eq.	1711	)		coeff = 	0.895835791
        if (JELEM	.eq.	1712	)		coeff = 	1.087155917
        if (JELEM	.eq.	1713	)		coeff = 	0.895775965
        if (JELEM	.eq.	1714	)		coeff = 	1.036332546
        if (JELEM	.eq.	1715	)		coeff = 	0.987823822
        if (JELEM	.eq.	1716	)		coeff = 	1.021728325
        if (JELEM	.eq.	1717	)		coeff = 	1.065523749
        if (JELEM	.eq.	1718	)		coeff = 	1.015982445
        if (JELEM	.eq.	1719	)		coeff = 	0.924597562
        if (JELEM	.eq.	1720	)		coeff = 	1.14252456
        if (JELEM	.eq.	1721	)		coeff = 	0.927543413
        if (JELEM	.eq.	1722	)		coeff = 	1.048677503
        if (JELEM	.eq.	1723	)		coeff = 	0.992092142
        if (JELEM	.eq.	1724	)		coeff = 	0.973431833
        if (JELEM	.eq.	1725	)		coeff = 	0.861502603
        if (JELEM	.eq.	1726	)		coeff = 	1.042042442
        if (JELEM	.eq.	1727	)		coeff = 	0.996063476
        if (JELEM	.eq.	1728	)		coeff = 	0.792601422
        if (JELEM	.eq.	1729	)		coeff = 	1.076422827
        if (JELEM	.eq.	1730	)		coeff = 	1.119782423
        if (JELEM	.eq.	1731	)		coeff = 	0.905757666
        if (JELEM	.eq.	1732	)		coeff = 	1.019266442
        if (JELEM	.eq.	1733	)		coeff = 	0.567187238
        if (JELEM	.eq.	1734	)		coeff = 	1.054072815
        if (JELEM	.eq.	1735	)		coeff = 	0.948366771
        if (JELEM	.eq.	1736	)		coeff = 	0.937673943
        if (JELEM	.eq.	1737	)		coeff = 	0.964280076
        if (JELEM	.eq.	1738	)		coeff = 	0.974326001
        if (JELEM	.eq.	1739	)		coeff = 	1.065095849
        if (JELEM	.eq.	1740	)		coeff = 	1.030666041
        if (JELEM	.eq.	1741	)		coeff = 	0.893748864
        if (JELEM	.eq.	1742	)		coeff = 	0.570610129
        if (JELEM	.eq.	1743	)		coeff = 	0.948319289
        if (JELEM	.eq.	1744	)		coeff = 	0.921297842
        if (JELEM	.eq.	1745	)		coeff = 	0.908839132
        if (JELEM	.eq.	1746	)		coeff = 	0.738636588
        if (JELEM	.eq.	1747	)		coeff = 	1.000402497
        if (JELEM	.eq.	1748	)		coeff = 	0.972503582
        if (JELEM	.eq.	1749	)		coeff = 	1.074502724
        if (JELEM	.eq.	1750	)		coeff = 	1.097185971
        if (JELEM	.eq.	1751	)		coeff = 	0.951583338
        if (JELEM	.eq.	1752	)		coeff = 	0.911766129
        if (JELEM	.eq.	1753	)		coeff = 	0.837708021
        if (JELEM	.eq.	1754	)		coeff = 	0.935544239
        if (JELEM	.eq.	1755	)		coeff = 	1.107422492
        if (JELEM	.eq.	1756	)		coeff = 	0.887828374
        if (JELEM	.eq.	1757	)		coeff = 	1.064793344
        if (JELEM	.eq.	1758	)		coeff = 	1.094362691
        if (JELEM	.eq.	1759	)		coeff = 	1.116709326
        if (JELEM	.eq.	1760	)		coeff = 	1.096862375
        if (JELEM	.eq.	1761	)		coeff = 	1.044469266
        if (JELEM	.eq.	1762	)		coeff = 	0.994397188
        if (JELEM	.eq.	1763	)		coeff = 	1.022143594
        if (JELEM	.eq.	1764	)		coeff = 	1.06623652
        if (JELEM	.eq.	1765	)		coeff = 	0.924476782
        if (JELEM	.eq.	1766	)		coeff = 	1.080765276
        if (JELEM	.eq.	1767	)		coeff = 	0.889074727
        if (JELEM	.eq.	1768	)		coeff = 	0.806558391
        if (JELEM	.eq.	1769	)		coeff = 	0.979677507
        if (JELEM	.eq.	1770	)		coeff = 	0.989109095
        if (JELEM	.eq.	1771	)		coeff = 	0.98754941
        if (JELEM	.eq.	1772	)		coeff = 	0.815010253
        if (JELEM	.eq.	1773	)		coeff = 	0.956374086
        if (JELEM	.eq.	1774	)		coeff = 	1.074827599
        if (JELEM	.eq.	1775	)		coeff = 	1.045026244
        if (JELEM	.eq.	1776	)		coeff = 	1.059494081
        if (JELEM	.eq.	1777	)		coeff = 	1.003583164
        if (JELEM	.eq.	1778	)		coeff = 	1.000903095
        if (JELEM	.eq.	1779	)		coeff = 	0.949372362
        if (JELEM	.eq.	1780	)		coeff = 	1.073367818
        if (JELEM	.eq.	1781	)		coeff = 	0.898246523
        if (JELEM	.eq.	1782	)		coeff = 	1.003511183
        if (JELEM	.eq.	1783	)		coeff = 	0.92018971
        if (JELEM	.eq.	1784	)		coeff = 	0.984456729
        if (JELEM	.eq.	1785	)		coeff = 	0.879825611
        if (JELEM	.eq.	1786	)		coeff = 	1.035310863
        if (JELEM	.eq.	1787	)		coeff = 	1.052554082
        if (JELEM	.eq.	1788	)		coeff = 	1.065108167
        if (JELEM	.eq.	1789	)		coeff = 	0.623119783
        if (JELEM	.eq.	1790	)		coeff = 	0.69316078
        if (JELEM	.eq.	1791	)		coeff = 	0.777861444
        if (JELEM	.eq.	1792	)		coeff = 	0.977395659
        if (JELEM	.eq.	1793	)		coeff = 	1.04815335
        if (JELEM	.eq.	1794	)		coeff = 	0.997912292
        if (JELEM	.eq.	1795	)		coeff = 	0.885628476
        if (JELEM	.eq.	1796	)		coeff = 	0.797555648
        if (JELEM	.eq.	1797	)		coeff = 	0.827360391
        if (JELEM	.eq.	1798	)		coeff = 	0.697036679
        if (JELEM	.eq.	1799	)		coeff = 	1.0655659
        if (JELEM	.eq.	1800	)		coeff = 	1.07672221



             
             
       DISP = 0.0
       E_EQV = 0.0
       DO I = 1, NNODE
        DISP(2*I-1,1) = U(3*I-2)
        DISP(2*I,1) = U(3*I-1)
							END DO
							DO I = 1, 4
        E_EQV(I,1) = U(3*I)
       END DO
C     ==================================================================
C     Saving time and increment no.
C     ==================================================================
       TIMEZ=USRVAR(JELEM,9,1)
       IF (TIMEZ.LT.TIME(2)) THEN
        USRVAR(JELEM,9,1)=TIME(2)
        USRVAR(JELEM,10,1)=0.0
       ELSE
        USRVAR(JELEM,10,1)=USRVAR(JELEM,10,1) + 1.0
       ENDIF
       STEPITER=USRVAR(JELEM,10,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
       EMOD = PROPS(1)
       ENU = PROPS(2)
       C = PROPS(3)
       AKI = PROPS(4) * coeff 
       beta = PROPS(5)
	 alpha = 0.99 
	 AK = 10.
	 h = EMOD/1D9
	 R = 0.05
	 Anta = 5.
       c1 = (AK-1.)/(2.*AK*(1.-2.*ENU))
       c2 = (AK-1.)**2./(1.-2.*ENU)**2.
       c3 = 2.*AK/(1.+ENU)**2.
	 del_oper(1:3) = (/1.0,1.0,0.0/)
	 dI1_de(1:3) = (/1.0,1.0,0.0/)
                              
       thickness = 250.0d0
C     ==================================================================
C     Calculating materials constitutive matrix (plane strain)
C     ==================================================================
       DMAT = 0.
       DMAT(1,1) = EMOD*(1.-ENU)/((1.+ENU)*(1.-2.*ENU))
       DMAT(2,2) = EMOD*(1.-ENU)/((1.+ENU)*(1.-2.*ENU))
       DMAT(3,3) = EMOD*(0.5-ENU)/((1.+ENU)*(1.-2.*ENU))
       DMAT(1,2) = EMOD*ENU/((1.+ENU)*(1.-2.*ENU))
       DMAT(2,1) = EMOD*ENU/((1.+ENU)*(1.-2.*ENU))
C     ==================================================================
C     Initialisation
C     ==================================================================
       RHS=0.;AMATRX=0.;AK_UU=0.;AK_EE=0.;AK_UE=0.;AK_EU=0.;F_EE=0.
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
C      Three Point Gauss Quadrature
       GW(1,1) = 5./9.; GW(2,1) = 5./9.
       GW(1,2) = 8./9.; GW(2,2) = 5./9.
       GW(1,3) = 5./9.; GW(2,3) = 5./9.
       GW(1,4) = 5./9.; GW(2,4) = 8./9.
       GW(1,5) = 8./9.; GW(2,5) = 8./9.
       GW(1,6) = 5./9.; GW(2,6) = 8./9.
       GW(1,7) = 5./9.; GW(2,7) = 5./9.
       GW(1,8) = 8./9.; GW(2,8) = 5./9.
       GW(1,9) = 5./9.; GW(2,9) = 5./9.
       GP(1,1) = -SQRT(3./5.); GP(2,1) = -SQRT(3./5.)
       GP(1,2) =  0.000000 ;   GP(2,2) = -SQRT(3./5.)
       GP(1,3) =  SQRT(3./5.); GP(2,3) = -SQRT(3./5.)
       GP(1,4) = -SQRT(3./5.); GP(2,4) =  0.0000000
       GP(1,5) =  0.000000 ;   GP(2,5) =  0.0000000
       GP(1,6) =  SQRT(3./5.); GP(2,6) =  0.0000000
       GP(1,7) = -SQRT(3./5.); GP(2,7) =  SQRT(3./5.)
       GP(1,8) =  0.000000 ;   GP(2,8) =  SQRT(3./5.)
       GP(1,9) =  SQRT(3./5.); GP(2,9) =  SQRT(3./5.)
C
      DO INPT=1,NGP!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GP LOOP BEGIN
C     ==================================================================
C      Shape Functions and derivatives
C     ==================================================================
C      Local coordinates of the integration point
       XI(1) = GP(1,INPT)
       XI(2) = GP(2,INPT) 
C      Shape functions and local derivatives
       CALL SHAPEFUN(AN,dNdxi,AH,dHdxi,XI)
C      Jacobian
       AJACOB = 0.0
       DO I = 1,2
        DO J = 1,2
         DO K = 1,NNODE
          AJACOB(I,J) = AJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
         END DO
        END DO
       END DO
C        
       DTM = 0.0
       DTM = AJACOB(1,1)*AJACOB(2,2)-AJACOB(1,2)*AJACOB(2,1)
       IF (DTM.LT.0.0) THEN
        WRITE(7,*) 'Negative Jacobian',DTM
        CALL XIT	
       END IF
C      Inverse of Jacobian
       AJABOBINV(1,1)=AJACOB(2,2)/DTM
       AJABOBINV(1,2)=-AJACOB(1,2)/DTM
       AJABOBINV(2,1)=-AJACOB(2,1)/DTM
       AJABOBINV(2,2)=AJACOB(1,1)/DTM       
C      Derivatives of shape functions Q8
       dNdx = 0.0
       DO K = 1,NNODE
        DO I = 1,2
         DO J = 1,2
          dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*AJABOBINV(J,I)
         END DO
        END DO
       END DO
C      Derivatives of shape functions Q4
       dHdx = 0.0
       DO K = 1,4
        DO I = 1,2
         DO J = 1,2
          dHdx(K,I) = dHdx(K,I) + dHdxi(K,J)*AJABOBINV(J,I)
         END DO
        END DO
       END DO
C
C      Calculating B matrix for disp and e_eqv
       BU = 0.0
       DO I = 1,NNODE
        BU(1,2*I-1) = dNdx(I,1)
        BU(2,2*I) = dNdx(I,2)
        BU(3,2*I-1) = dNdx(I,2)
        BU(3,2*I) = dNdx(I,1)
						 END DO
C
       BE = 0.0
						 DO I = 1,4
        BE(1,I) = dHdx(I,1)
        BE(2,I) = dHdx(I,2)
						 END DO
C
C     ==================================================================
C      Updating SVARS and USRVAR
C     ==================================================================
C      Transfering maximum of last LOAD increment K_GP->K0
       IF (STEPITER .EQ. 0.0) THEN
        AK0(1,INPT) = USRVAR(JELEM,6,INPT)
       ELSE
        AK0(1,INPT) = USRVAR(JELEM,7,INPT)
       ENDIF
C      Strains{1ST-3RD SDV}.............................................
       E(1,INPT) = DOT_PRODUCT(BU(1,:),DISP(:,1))
       E(2,INPT) = DOT_PRODUCT(BU(2,:),DISP(:,1))
       E(3,INPT) = DOT_PRODUCT(BU(3,:),DISP(:,1))
       SVARS(NSDV*(INPT-1)+1) = E(1,INPT)
       SVARS(NSDV*(INPT-1)+2) = E(2,INPT)
       SVARS(NSDV*(INPT-1)+3) = E(3,INPT)
       USRVAR(JELEM,1,INPT) = E(1,INPT)
       USRVAR(JELEM,2,INPT) = E(2,INPT)
       USRVAR(JELEM,3,INPT) = E(3,INPT)
C      ALOC_E_EQV{4TH SDV}..............................................
       ALOC_E_EQV(1,INPT) = c1*(E(1,INPT) + E(2,INPT)) + (1./(2.*AK))*
     &  SQRT((c2*(E(1,INPT)+E(2,INPT))**2.)+(c3*(E(1,INPT)**2. + 
     &   E(2,INPT)**2 - E(1,INPT)*E(2,INPT) + 3.*(E(3,INPT)**2.))))
					  SVARS(NSDV*(INPT-1)+4) = ALOC_E_EQV(1,INPT)
							USRVAR(JELEM,4,INPT) = ALOC_E_EQV(1,INPT)
C      GP_E_EQV(5TH SDV)................................................
       GP_E_EQV(1,INPT) = DOT_PRODUCT(AH(1,:),E_EQV(:,1))
       SVARS(NSDV*(INPT-1)+5) = GP_E_EQV(1,INPT)
       USRVAR(JELEM,5,INPT) = GP_E_EQV(1,INPT)
C      AK_GP(6TH SDV)...................................................
							IF (GP_E_EQV(1,INPT) .LT. AK0(1,INPT)) THEN
							  AK_GP(1,INPT) = AK0(1,INPT)
							ELSE
         AK_GP(1,INPT) = GP_E_EQV(1,INPT)
       END IF
							SVARS(NSDV*(INPT-1)+6) = AK_GP(1,INPT)
       USRVAR(JELEM,6,INPT) = AK_GP(1,INPT)
C      AK0(7th SDV).....................................................
							SVARS(NSDV*(INPT-1)+7) = AK0(1,INPT)!AK0_OUT
       USRVAR(JELEM,7,INPT) = AK0(1,INPT)!AK0_OUT
C      DAMAGE(8TH SDV)..................................................	
       D_OLD = SVARS(NSDV*(INPT-1)+8)
       
       IF (AK_GP(1,INPT) .LT. AKI) THEN
							  D = 0.
	 ELSE
         D = 1.-(AKI/AK_GP(1,INPT))*(1.-alpha+alpha*
     &           DEXP(-beta*(AK_GP(1,INPT)-AKI)))
       END IF
	 SVARS(NSDV*(INPT-1)+8) = D
       USRVAR(JELEM,8,INPT)   = D
       
       DELTA_D = D - D_OLD
       
       IF (DELTA_D .GE. 0.1D0) THEN
           pnewdt = 0.25
       END IF
C     ==================================================================
C     Computing parameters for the integration point
C     ==================================================================
      AI1 = E(1,INPT) + E(2,INPT)
C						
						AJ2 = 2.*(E(1,INPT)**2. + E(2,INPT)**2. + 3.*E(3,INPT) - 
     &						  E(1,INPT)*E(2,INPT))
C
					 dJ2_de(1) = 4.*E(1,INPT) - 2.*E(2,INPT)
						dJ2_de(2) = 4.*E(2,INPT) - 2.*E(1,INPT)
						dJ2_de(3) = 12.*E(3,INPT)
C
      A_term=(((AK-1.)*AI1/(1.-2.*ENU))**2.)+(2.*AK/(1. + ENU)**2.)*AJ2
C
      IF (A_term .LT. 1D-10) THEN
						  H_mat = 0.
						ELSE
						  H_mat = RESHAPE((/2.,-1.,0.,-1.,2.,0.,0.,0.,1.5/), (/3, 3/))
						  H_mat = H_mat*(h/(((1. + ENU)**2.)*SQRT(A_term)))
      END IF								
      !Derivative of damage wrt history parameter 
						IF (AK_GP(1,INPT) .LT. AKI) THEN
						 dD_dK = 0.
      ELSEIF (AK_GP(1,INPT) .GT. AK0(1,INPT)) THEN
       dD_dk = ((AKI/(AK_GP(1,INPT)**2.))*
     & (1.-alpha+alpha*DEXP(-beta*(AK_GP(1,INPT)-AKI)))) +
     & ((alpha*beta*AKI*DEXP(-beta*(AK_GP(1,INPT)-AKI)))/AK_GP(1,INPT))
      ELSE
       dD_dk = 0.
      END IF
						! Parameter G
						G = (((1. - R)*DEXP(-Anta*D)) + R - DEXP(-Anta)) /
     &    (1. - DEXP(-Anta))
						dG_dk = dD_dk*(Anta*(R - 1.)*DEXP(-Anta*D))/
     &        (1. - DEXP(-Anta))
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
      T_INT = GW(1,INPT)*GW(2,INPT)*DTM * thickness 
						del_E = MATMUL(BE,E_EQV)
C     K_UU..............................................................
       DO K=1,2*NNODE
        DO L=1,2*NNODE
         DO I=1,3
          DO J=1,3
           AK_UU(K,L)=AK_UU(K,L)+BU(I,K)*DMAT(I,J)*BU(J,L)*(1.-D)*T_INT
          END DO
         END DO
        END DO
       END DO
C     K_EE..............................................................
       DO K=1,4
        DO L=1,4
         AK_EE(K,L) = AK_EE(K,L) + AH(1,K)*AH(1,L)*h*T_INT
         DO I=1,2
          AK_EE(K,L) = AK_EE(K,L) + BE(I,K)*del_E(I,1)*
     &     AH(1,L)*h*C*dG_dk*T_INT
          AK_EE(K,L) = AK_EE(K,L) + BE(I,K)*BE(I,L)*G*h*C*T_INT
         END DO
        END DO
       END DO
C     K_UE..............................................................
       DO K = 1, 2*NNODE
        DO L = 1, 4
         DO I = 1, 3
										AK_UE(K,L) = AK_UE(K,L) - 
     &                 BU(I,K)*h*del_oper(I)*AH(1,L)*T_INT
          DO J = 1, 3
									  AK_UE(K,L) = AK_UE(K,L) - 
     &                  BU(I,K)*H_mat(I,J)*E(J,INPT)*AH(1,L)*T_INT
									  AK_UE(K,L) = AK_UE(K,L) - 
     &                  BU(I,K)*DMAT(I,J)*E(J,INPT)*AH(1,L)*dD_dK*T_INT
					     END DO
         END DO
        END DO
       END DO
C     K_EU..............................................................
       Denominator = c2*(E(1,INPT) + E(2,INPT))**2. + 
     &     (c3*(E(1,INPT)**2. + E(2,INPT)**2. - 
     &     E(1,INPT)*E(2,INPT)+3.*(E(3,INPT)**2.)))
       IF (Denominator .LT. 1D-10) THEN ! for avoiding NaN
								DEEQV_DE(:) = 0.0
							ELSE
C      Computing derivatives of principal strains
       DEEQV_DE(1) = c1 + (c3*(2.*E(1,INPT) - E(2,INPT))
     &     + 2.*c2*(E(1,INPT) + E(2,INPT))) / 
     &     (4.*AK*(Denominator)**0.5)
       DEEQV_DE(2) = c1 + (c3*(2.*E(2,INPT) - E(1,INPT))
     &     + 2.*c2*(E(1,INPT) + E(2,INPT))) / 
     &     (4.*AK*(Denominator)**0.5)
       DEEQV_DE(3) = (3.*c3*E(3,INPT)) /
     &     (2.*AK*(Denominator)**0.5)
							END IF
							DO K=1,4
        DO L=1,2*NNODE
         DO J=1,3
          AK_EU(K,L) = AK_EU(K,L) - AH(1,K)*h*DEEQV_DE(J)*BU(J,L)*T_INT
         END DO
        END DO
       END DO
C     F_EE..............................................................	 
       DO K = 1,4
        F_EE(K,1) = F_EE(K,1) - (AH(1,K)*h*(GP_E_EQV(1,INPT)
     &   -ALOC_E_EQV(1,INPT)))*T_INT
					   DO I = 1, 2
					    F_EE(K,1) = F_EE(K,1) - BE(I,K)*del_E(I,1)*G*h*C*T_INT
					   END DO
       END DO
C
C       WRITE(6,201) D
       
       
       
      USRVAR(JELEM,11,INPT) = coeff
      USRVAR(JELEM,12,INPT) = DELTA_D
              
              
              
      END DO!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GP LOOP END
C     Assembly of K into AMATRX........................................
      !Assemble K_UU
      DO I=1,NNODE
       DO J=1,NNODE
        AMATRX(3*I-2,3*J-2) = AK_UU(2*I-1,2*J-1)
        AMATRX(3*I-1,3*J-2) = AK_UU(2*I,2*J-1)
        AMATRX(3*I-2,3*J-1) = AK_UU(2*I-1,2*J)
        AMATRX(3*I-1,3*J-1) = AK_UU(2*I,2*J)
							END DO
						END DO
		    !Assemble K_UE
						DO I=1,NNODE
							DO J=1,4
        AMATRX(3*I-2,3*J) = AK_UE(2*I-1,J)
        AMATRX(3*I-1,3*J) = AK_UE(2*I,J)
							END DO
						END DO
		    !Assemble K_EU
						DO I=1,4
							DO J=1,NNODE
        AMATRX(3*I,3*J-2) = AK_EU(I,2*J-1)
        AMATRX(3*I,3*J-1) = AK_EU(I,2*J)
							END DO
						END DO
		    !Assemble K_EE
						DO I=1,4
							DO J=1,4
        AMATRX(3*I,3*J) = AK_EE(I,J)
       END DO
      END DO
						!Conditioning
						AMATRX(15,15) = 1.0; AMATRX(18,18) = 1.0
						AMATRX(21,21) = 1.0; AMATRX(24,24) = 1.0
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
       RHS = 0.0
	      !FINT_U
       DO I=1,NNODE
        DO J=1,2*NNODE
         RHS(3*I-2,1) = RHS(3*I-2,1)-AK_UU(2*I-1,J)*(DISP(J,1))
         RHS(3*I-1,1) = RHS(3*I-1,1)-AK_UU(2*I,J)*(DISP(J,1))
        END DO
							END DO
	      !FINT_E
       DO I=1,4
        RHS(3*I,1) = RHS(3*I,1) + F_EE(I,1)
       END DO
C
      RETURN
      END
C***********************************************************************
C     ==================================================================
C     SHAPE FUNCTION SUBROUTINE
C     ==================================================================
      SUBROUTINE SHAPEFUN(AN,dNdxi,AH,dHdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(1,8),dNdxi(8,2)
						Real*8 AH(1,4),dHdxi(4,2)
      Real*8 XI(2)
C
C     Values of shape functions for Q8
						AN(1,1) = -0.25*(1.-XI(1))*(1.-XI(2))*(1.+XI(1)+XI(2))
      AN(1,2) = -0.25*(1.+XI(1))*(1.-XI(2))*(1.-XI(1)+XI(2))
      AN(1,3) = -0.25*(1.+XI(1))*(1.+XI(2))*(1.-XI(1)-XI(2))
      AN(1,4) = -0.25*(1.-XI(1))*(1.+XI(2))*(1.+XI(1)-XI(2))
      AN(1,5) =  0.5*(1.-XI(1))*(1.+XI(1))*(1.-XI(2))
      AN(1,6) =  0.5*(1.+XI(1))*(1.+XI(2))*(1.-XI(2))
      AN(1,7) =  0.5*(1.-XI(1))*(1.+XI(1))*(1.+XI(2))
      AN(1,8) =  0.5*(1.-XI(1))*(1.+XI(2))*(1.-XI(2))
C
C     Values of shape functions for Q4
      AH(1,1) = 0.25*(1.-XI(1))*(1.-XI(2))
      AH(1,2) = 0.25*(1.+XI(1))*(1.-XI(2))
      AH(1,3) = 0.25*(1.+XI(1))*(1.+XI(2))
      AH(1,4) = 0.25*(1.-XI(1))*(1.+XI(2))
C
C     Derivatives of shape functions Q8
      dNdxi = 0.0
      dNdxi(1,1) = -0.25*(-1.+XI(2))*(2.*XI(1)+XI(2))
      dNdxi(1,2) = -0.25*(-1.+XI(1))*(XI(1)+2.*XI(2))
      dNdxi(2,1) =  0.25*(-1.+XI(2))*(XI(2)-2.*XI(1))
      dNdxi(2,2) =  0.25*(1.+XI(1))*(2.*XI(2)-XI(1))
      dNdxi(3,1) =  0.25*(1.+XI(2))*(2.*XI(1)+XI(2))
      dNdxi(3,2) =  0.25*(1.+XI(1))*(XI(1)+2.*XI(2))
      dNdxi(4,1) = -0.25*(1.+XI(2))*(XI(2)-2.*XI(1))
      dNdxi(4,2) = -0.25*(-1.+XI(1))*(2.*XI(2)-XI(1))
						dNdxi(5,1) =  XI(1)*(-1.+XI(2))  
      dNdxi(5,2) =  0.5*(1.+XI(1))*(-1.+XI(1))
      dNdxi(6,1) = -0.5*(1.+XI(2))*(-1.+XI(2)) 
      dNdxi(6,2) = -XI(2)*(1.+XI(1))
      dNdxi(7,1) = -XI(1)*(1.+XI(2)) 
      dNdxi(7,2) = -0.5*(1.+XI(1))*(-1.+XI(1))
      dNdxi(8,1) =  0.5*(1.+XI(2))*(-1.+XI(2)) 
      dNdxi(8,2) =  XI(2)*(-1.+XI(1))
C
C     Derivatives of shape functions Q4
      dHdxi = 0.0
      dHdxi(1,1) = -0.25*(1.-XI(2))
      dHdxi(1,2) = -0.25*(1.-XI(1))
      dHdxi(2,1) =  0.25*(1.-XI(2))
      dHdxi(2,2) = -0.25*(1.+XI(1))
      dHdxi(3,1) =  0.25*(1.+XI(2))
      dHdxi(3,2) =  0.25*(1.+XI(1))
      dHdxi(4,1) = -0.25*(1.+XI(2))
      dHdxi(4,2) =  0.25*(1.-XI(1))
      RETURN
      END
C***********************************************************************      
C Subroutine UMAT  : 
C Dummy material
C
C ==============================================================
C !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
C ==============================================================
C
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=10000,NSDV=120)
       DATA NEWTON,TOLER/40,1.D-2/ 
C       
       COMMON/KUSER/USRVAR(N_ELEM,NSDV,9)
C 
C ----------------------------------------------------------- 
C          Material properties
C ----------------------------------------------------------- 
C          PROPS(1) - Young's modulus 
C          PROPS(2) - Poisson ratio 
C ----------------------------------------------------------- 
C
C	Elastic properties
C
       EMOD=PROPS(1)
       ENU=PROPS(2)
       EG=EMOD/(TWO*(ONE+ENU))
       EG2=EG*TWO
       ELAM=EG2*ENU/(ONE-TWO*ENU)
C
C	Stiffness tensor
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         DDSDDE(K2, K1)=0.0
        END DO
       END DO
C
       DO K1=1, NDI
        DO K2=1, NDI
         DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
       END DO 
C
       DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
       END DO
C
C	Calculate Stresses
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
       END DO 
C
       NELEMAN=NOEL-N_ELEM

       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
C       
       RETURN
       END      