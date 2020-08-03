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
       
        if (JELEM	.eq.	1	)		coeff = 	0.958232231
        if (JELEM	.eq.	2	)		coeff = 	0.843382496
        if (JELEM	.eq.	3	)		coeff = 	1.084378428
        if (JELEM	.eq.	4	)		coeff = 	0.960185612
        if (JELEM	.eq.	5	)		coeff = 	1.128663003
        if (JELEM	.eq.	6	)		coeff = 	0.836163973
        if (JELEM	.eq.	7	)		coeff = 	0.991247058
        if (JELEM	.eq.	8	)		coeff = 	1.046929199
        if (JELEM	.eq.	9	)		coeff = 	1.122396102
        if (JELEM	.eq.	10	)		coeff = 	1.123483604
        if (JELEM	.eq.	11	)		coeff = 	0.723815338
        if (JELEM	.eq.	12	)		coeff = 	1.05935742
        if (JELEM	.eq.	13	)		coeff = 	0.677502564
        if (JELEM	.eq.	14	)		coeff = 	1.009616145
        if (JELEM	.eq.	15	)		coeff = 	0.670514522
        if (JELEM	.eq.	16	)		coeff = 	1.008445952
        if (JELEM	.eq.	17	)		coeff = 	1.135201431
        if (JELEM	.eq.	18	)		coeff = 	0.914749542
        if (JELEM	.eq.	19	)		coeff = 	0.935713406
        if (JELEM	.eq.	20	)		coeff = 	0.74821197
        if (JELEM	.eq.	21	)		coeff = 	1.060786203
        if (JELEM	.eq.	22	)		coeff = 	0.836025083
        if (JELEM	.eq.	23	)		coeff = 	1.006877211
        if (JELEM	.eq.	24	)		coeff = 	0.640274717
        if (JELEM	.eq.	25	)		coeff = 	0.806540206
        if (JELEM	.eq.	26	)		coeff = 	0.873830042
        if (JELEM	.eq.	27	)		coeff = 	1.055248837
        if (JELEM	.eq.	28	)		coeff = 	1.19542101
        if (JELEM	.eq.	29	)		coeff = 	1.05719688
        if (JELEM	.eq.	30	)		coeff = 	0.942516601
        if (JELEM	.eq.	31	)		coeff = 	1.104913144
        if (JELEM	.eq.	32	)		coeff = 	0.794086065
        if (JELEM	.eq.	33	)		coeff = 	0.743497135
        if (JELEM	.eq.	34	)		coeff = 	1.011690007
        if (JELEM	.eq.	35	)		coeff = 	0.965775142
        if (JELEM	.eq.	36	)		coeff = 	1.246923925
        if (JELEM	.eq.	37	)		coeff = 	1.019629324
        if (JELEM	.eq.	38	)		coeff = 	0.960292243
        if (JELEM	.eq.	39	)		coeff = 	1.037142992
        if (JELEM	.eq.	40	)		coeff = 	0.886399281
        if (JELEM	.eq.	41	)		coeff = 	0.743130955
        if (JELEM	.eq.	42	)		coeff = 	0.930343882
        if (JELEM	.eq.	43	)		coeff = 	1.035995625
        if (JELEM	.eq.	44	)		coeff = 	0.722456416
        if (JELEM	.eq.	45	)		coeff = 	1.134537579
        if (JELEM	.eq.	46	)		coeff = 	0.909435312
        if (JELEM	.eq.	47	)		coeff = 	0.692857286
        if (JELEM	.eq.	48	)		coeff = 	0.949545664
        if (JELEM	.eq.	49	)		coeff = 	0.677256911
        if (JELEM	.eq.	50	)		coeff = 	1.033525421
        if (JELEM	.eq.	51	)		coeff = 	0.945847527
        if (JELEM	.eq.	52	)		coeff = 	0.743123602
        if (JELEM	.eq.	53	)		coeff = 	1.05533442
        if (JELEM	.eq.	54	)		coeff = 	0.812651298
        if (JELEM	.eq.	55	)		coeff = 	1.089009941
        if (JELEM	.eq.	56	)		coeff = 	0.976403532
        if (JELEM	.eq.	57	)		coeff = 	0.777260225
        if (JELEM	.eq.	58	)		coeff = 	0.877184963
        if (JELEM	.eq.	59	)		coeff = 	0.95624464
        if (JELEM	.eq.	60	)		coeff = 	0.915842893
        if (JELEM	.eq.	61	)		coeff = 	1.124311475
        if (JELEM	.eq.	62	)		coeff = 	1.044913993
        if (JELEM	.eq.	63	)		coeff = 	1.175994769
        if (JELEM	.eq.	64	)		coeff = 	0.799342559
        if (JELEM	.eq.	65	)		coeff = 	1.196798399
        if (JELEM	.eq.	66	)		coeff = 	0.937571005
        if (JELEM	.eq.	67	)		coeff = 	1.098674247
        if (JELEM	.eq.	68	)		coeff = 	0.923950183
        if (JELEM	.eq.	69	)		coeff = 	0.911017655
        if (JELEM	.eq.	70	)		coeff = 	0.857256528
        if (JELEM	.eq.	71	)		coeff = 	1.019973401
        if (JELEM	.eq.	72	)		coeff = 	0.737152837
        if (JELEM	.eq.	73	)		coeff = 	0.839015161
        if (JELEM	.eq.	74	)		coeff = 	1.109253722
        if (JELEM	.eq.	75	)		coeff = 	1.03022937
        if (JELEM	.eq.	76	)		coeff = 	0.747098292
        if (JELEM	.eq.	77	)		coeff = 	1.103416903
        if (JELEM	.eq.	78	)		coeff = 	0.834838919
        if (JELEM	.eq.	79	)		coeff = 	1.108518483
        if (JELEM	.eq.	80	)		coeff = 	0.846658805
        if (JELEM	.eq.	81	)		coeff = 	0.895190825
        if (JELEM	.eq.	82	)		coeff = 	0.84394747
        if (JELEM	.eq.	83	)		coeff = 	0.717760014
        if (JELEM	.eq.	84	)		coeff = 	1.071091552
        if (JELEM	.eq.	85	)		coeff = 	1.185113118
        if (JELEM	.eq.	86	)		coeff = 	0.808156509
        if (JELEM	.eq.	87	)		coeff = 	0.964397088
        if (JELEM	.eq.	88	)		coeff = 	1.104805317
        if (JELEM	.eq.	89	)		coeff = 	1.116717812
        if (JELEM	.eq.	90	)		coeff = 	1.101157271
        if (JELEM	.eq.	91	)		coeff = 	1.075035788
        if (JELEM	.eq.	92	)		coeff = 	0.63922099
        if (JELEM	.eq.	93	)		coeff = 	1.028465512
        if (JELEM	.eq.	94	)		coeff = 	0.919828449
        if (JELEM	.eq.	95	)		coeff = 	0.795222819
        if (JELEM	.eq.	96	)		coeff = 	1.043641417
        if (JELEM	.eq.	97	)		coeff = 	0.866772737
        if (JELEM	.eq.	98	)		coeff = 	0.988865066
        if (JELEM	.eq.	99	)		coeff = 	1.019737308
        if (JELEM	.eq.	100	)		coeff = 	0.898551323
        if (JELEM	.eq.	101	)		coeff = 	0.99189098
        if (JELEM	.eq.	102	)		coeff = 	1.09885555
        if (JELEM	.eq.	103	)		coeff = 	0.619506965
        if (JELEM	.eq.	104	)		coeff = 	0.948898945
        if (JELEM	.eq.	105	)		coeff = 	0.952537398
        if (JELEM	.eq.	106	)		coeff = 	1.014703506
        if (JELEM	.eq.	107	)		coeff = 	1.062418449
        if (JELEM	.eq.	108	)		coeff = 	0.787201689
        if (JELEM	.eq.	109	)		coeff = 	1.035044286
        if (JELEM	.eq.	110	)		coeff = 	0.941735735
        if (JELEM	.eq.	111	)		coeff = 	0.937196519
        if (JELEM	.eq.	112	)		coeff = 	0.663223357
        if (JELEM	.eq.	113	)		coeff = 	1.034418964
        if (JELEM	.eq.	114	)		coeff = 	1.035833901
        if (JELEM	.eq.	115	)		coeff = 	1.019951737
        if (JELEM	.eq.	116	)		coeff = 	0.935109643
        if (JELEM	.eq.	117	)		coeff = 	0.996502541
        if (JELEM	.eq.	118	)		coeff = 	1.043270389
        if (JELEM	.eq.	119	)		coeff = 	0.889330839
        if (JELEM	.eq.	120	)		coeff = 	1.037464937
        if (JELEM	.eq.	121	)		coeff = 	0.667213053
        if (JELEM	.eq.	122	)		coeff = 	0.956053853
        if (JELEM	.eq.	123	)		coeff = 	0.979769603
        if (JELEM	.eq.	124	)		coeff = 	0.931326699
        if (JELEM	.eq.	125	)		coeff = 	0.914534229
        if (JELEM	.eq.	126	)		coeff = 	0.824922629
        if (JELEM	.eq.	127	)		coeff = 	0.772593976
        if (JELEM	.eq.	128	)		coeff = 	1.078665845
        if (JELEM	.eq.	129	)		coeff = 	1.116638423
        if (JELEM	.eq.	130	)		coeff = 	1.018848801
        if (JELEM	.eq.	131	)		coeff = 	0.916693161
        if (JELEM	.eq.	132	)		coeff = 	0.996316547
        if (JELEM	.eq.	133	)		coeff = 	0.726231431
        if (JELEM	.eq.	134	)		coeff = 	1.140592434
        if (JELEM	.eq.	135	)		coeff = 	0.5275791
        if (JELEM	.eq.	136	)		coeff = 	0.893173118
        if (JELEM	.eq.	137	)		coeff = 	0.912963213
        if (JELEM	.eq.	138	)		coeff = 	1.084342641
        if (JELEM	.eq.	139	)		coeff = 	1.171209198
        if (JELEM	.eq.	140	)		coeff = 	0.929628152
        if (JELEM	.eq.	141	)		coeff = 	1.169821621
        if (JELEM	.eq.	142	)		coeff = 	0.873952574
        if (JELEM	.eq.	143	)		coeff = 	1.055592914
        if (JELEM	.eq.	144	)		coeff = 	1.039732813
        if (JELEM	.eq.	145	)		coeff = 	0.673407626
        if (JELEM	.eq.	146	)		coeff = 	1.068919135
        if (JELEM	.eq.	147	)		coeff = 	0.829928788
        if (JELEM	.eq.	148	)		coeff = 	1.157297598
        if (JELEM	.eq.	149	)		coeff = 	1.020243445
        if (JELEM	.eq.	150	)		coeff = 	0.653984803
        if (JELEM	.eq.	151	)		coeff = 	1.035366436
        if (JELEM	.eq.	152	)		coeff = 	1.097489028
        if (JELEM	.eq.	153	)		coeff = 	0.916683681
        if (JELEM	.eq.	154	)		coeff = 	0.901943107
        if (JELEM	.eq.	155	)		coeff = 	0.777973288
        if (JELEM	.eq.	156	)		coeff = 	0.757434078
        if (JELEM	.eq.	157	)		coeff = 	1.144022884
        if (JELEM	.eq.	158	)		coeff = 	0.851004574
        if (JELEM	.eq.	159	)		coeff = 	0.790251804
        if (JELEM	.eq.	160	)		coeff = 	0.961488351
        if (JELEM	.eq.	161	)		coeff = 	1.126266369
        if (JELEM	.eq.	162	)		coeff = 	1.115643884
        if (JELEM	.eq.	163	)		coeff = 	1.061201575
        if (JELEM	.eq.	164	)		coeff = 	0.841173728
        if (JELEM	.eq.	165	)		coeff = 	0.986676948
        if (JELEM	.eq.	166	)		coeff = 	0.812485625
        if (JELEM	.eq.	167	)		coeff = 	1.084303232
        if (JELEM	.eq.	168	)		coeff = 	0.864066515
        if (JELEM	.eq.	169	)		coeff = 	1.083384218
        if (JELEM	.eq.	170	)		coeff = 	0.804929904
        if (JELEM	.eq.	171	)		coeff = 	0.939367704
        if (JELEM	.eq.	172	)		coeff = 	0.667913193
        if (JELEM	.eq.	173	)		coeff = 	1.125978125
        if (JELEM	.eq.	174	)		coeff = 	0.715199861
        if (JELEM	.eq.	175	)		coeff = 	1.002877392
        if (JELEM	.eq.	176	)		coeff = 	1.132407886
        if (JELEM	.eq.	177	)		coeff = 	0.656869391
        if (JELEM	.eq.	178	)		coeff = 	0.900302214
        if (JELEM	.eq.	179	)		coeff = 	0.982181628
        if (JELEM	.eq.	180	)		coeff = 	0.920794529
        if (JELEM	.eq.	181	)		coeff = 	1.007093334
        if (JELEM	.eq.	182	)		coeff = 	1.054163124
        if (JELEM	.eq.	183	)		coeff = 	1.159820307
        if (JELEM	.eq.	184	)		coeff = 	1.081712521
        if (JELEM	.eq.	185	)		coeff = 	0.998131908
        if (JELEM	.eq.	186	)		coeff = 	0.775977065
        if (JELEM	.eq.	187	)		coeff = 	0.975256802
        if (JELEM	.eq.	188	)		coeff = 	0.904950367
        if (JELEM	.eq.	189	)		coeff = 	0.993248641
        if (JELEM	.eq.	190	)		coeff = 	0.80576801
        if (JELEM	.eq.	191	)		coeff = 	0.847505732
        if (JELEM	.eq.	192	)		coeff = 	1.095522904
        if (JELEM	.eq.	193	)		coeff = 	1.114409692
        if (JELEM	.eq.	194	)		coeff = 	0.895901719
        if (JELEM	.eq.	195	)		coeff = 	0.996678321
        if (JELEM	.eq.	196	)		coeff = 	1.105295492
        if (JELEM	.eq.	197	)		coeff = 	0.910399805
        if (JELEM	.eq.	198	)		coeff = 	0.911252999
        if (JELEM	.eq.	199	)		coeff = 	0.742645807
        if (JELEM	.eq.	200	)		coeff = 	0.917154779
        if (JELEM	.eq.	201	)		coeff = 	0.960250046
        if (JELEM	.eq.	202	)		coeff = 	0.928070687
        if (JELEM	.eq.	203	)		coeff = 	0.605469959
        if (JELEM	.eq.	204	)		coeff = 	0.948788731
        if (JELEM	.eq.	205	)		coeff = 	1.006540595
        if (JELEM	.eq.	206	)		coeff = 	0.60210772
        if (JELEM	.eq.	207	)		coeff = 	0.547058815
        if (JELEM	.eq.	208	)		coeff = 	0.758570097
        if (JELEM	.eq.	209	)		coeff = 	0.859911532
        if (JELEM	.eq.	210	)		coeff = 	0.894089467
        if (JELEM	.eq.	211	)		coeff = 	0.863853957
        if (JELEM	.eq.	212	)		coeff = 	0.977022817
        if (JELEM	.eq.	213	)		coeff = 	0.98925588
        if (JELEM	.eq.	214	)		coeff = 	1.054058097
        if (JELEM	.eq.	215	)		coeff = 	0.771860675
        if (JELEM	.eq.	216	)		coeff = 	0.703031869
        if (JELEM	.eq.	217	)		coeff = 	0.691688004
        if (JELEM	.eq.	218	)		coeff = 	1.106074984
        if (JELEM	.eq.	219	)		coeff = 	1.139757608
        if (JELEM	.eq.	220	)		coeff = 	1.076953687
        if (JELEM	.eq.	221	)		coeff = 	0.945678025
        if (JELEM	.eq.	222	)		coeff = 	1.130234729
        if (JELEM	.eq.	223	)		coeff = 	1.209859304
        if (JELEM	.eq.	224	)		coeff = 	1.115031818
        if (JELEM	.eq.	225	)		coeff = 	0.916662031
        if (JELEM	.eq.	226	)		coeff = 	1.141407204
        if (JELEM	.eq.	227	)		coeff = 	0.89680633
        if (JELEM	.eq.	228	)		coeff = 	0.928407169
        if (JELEM	.eq.	229	)		coeff = 	1.127184695
        if (JELEM	.eq.	230	)		coeff = 	1.044862017
        if (JELEM	.eq.	231	)		coeff = 	1.005687193
        if (JELEM	.eq.	232	)		coeff = 	1.094917793
        if (JELEM	.eq.	233	)		coeff = 	0.860961452
        if (JELEM	.eq.	234	)		coeff = 	0.973726069
        if (JELEM	.eq.	235	)		coeff = 	1.190657559
        if (JELEM	.eq.	236	)		coeff = 	0.886995861
        if (JELEM	.eq.	237	)		coeff = 	0.772551069
        if (JELEM	.eq.	238	)		coeff = 	1.088269268
        if (JELEM	.eq.	239	)		coeff = 	0.8226051
        if (JELEM	.eq.	240	)		coeff = 	1.04060211
        if (JELEM	.eq.	241	)		coeff = 	1.178220934
        if (JELEM	.eq.	242	)		coeff = 	0.944538192
        if (JELEM	.eq.	243	)		coeff = 	0.713762352
        if (JELEM	.eq.	244	)		coeff = 	0.845934605
        if (JELEM	.eq.	245	)		coeff = 	0.987529026
        if (JELEM	.eq.	246	)		coeff = 	0.871100932
        if (JELEM	.eq.	247	)		coeff = 	0.961625169
        if (JELEM	.eq.	248	)		coeff = 	0.975247522
        if (JELEM	.eq.	249	)		coeff = 	0.702575769
        if (JELEM	.eq.	250	)		coeff = 	0.927325281
        if (JELEM	.eq.	251	)		coeff = 	0.838279908
        if (JELEM	.eq.	252	)		coeff = 	0.87247062
        if (JELEM	.eq.	253	)		coeff = 	1.087114883
        if (JELEM	.eq.	254	)		coeff = 	0.882185684
        if (JELEM	.eq.	255	)		coeff = 	0.991275136
        if (JELEM	.eq.	256	)		coeff = 	1.094381008
        if (JELEM	.eq.	257	)		coeff = 	0.903711373
        if (JELEM	.eq.	258	)		coeff = 	0.956077335
        if (JELEM	.eq.	259	)		coeff = 	1.134756816
        if (JELEM	.eq.	260	)		coeff = 	0.947055407
        if (JELEM	.eq.	261	)		coeff = 	1.018356846
        if (JELEM	.eq.	262	)		coeff = 	0.79325241
        if (JELEM	.eq.	263	)		coeff = 	1.213708503
        if (JELEM	.eq.	264	)		coeff = 	1.060228049
        if (JELEM	.eq.	265	)		coeff = 	0.7727486
        if (JELEM	.eq.	266	)		coeff = 	1.053995211
        if (JELEM	.eq.	267	)		coeff = 	1.017256069
        if (JELEM	.eq.	268	)		coeff = 	1.008165201
        if (JELEM	.eq.	269	)		coeff = 	0.828990401
        if (JELEM	.eq.	270	)		coeff = 	0.930132132
        if (JELEM	.eq.	271	)		coeff = 	0.821156011
        if (JELEM	.eq.	272	)		coeff = 	1.118600578
        if (JELEM	.eq.	273	)		coeff = 	1.064560851
        if (JELEM	.eq.	274	)		coeff = 	0.908191241
        if (JELEM	.eq.	275	)		coeff = 	1.134844825
        if (JELEM	.eq.	276	)		coeff = 	0.859240324
        if (JELEM	.eq.	277	)		coeff = 	0.918186049
        if (JELEM	.eq.	278	)		coeff = 	0.829935909
        if (JELEM	.eq.	279	)		coeff = 	1.078382642
        if (JELEM	.eq.	280	)		coeff = 	0.861672573
        if (JELEM	.eq.	281	)		coeff = 	0.852259632
        if (JELEM	.eq.	282	)		coeff = 	0.772137396
        if (JELEM	.eq.	283	)		coeff = 	0.979274048
        if (JELEM	.eq.	284	)		coeff = 	1.066526488
        if (JELEM	.eq.	285	)		coeff = 	1.022402321
        if (JELEM	.eq.	286	)		coeff = 	1.016186591
        if (JELEM	.eq.	287	)		coeff = 	0.974966931
        if (JELEM	.eq.	288	)		coeff = 	1.002039069
        if (JELEM	.eq.	289	)		coeff = 	0.909707338
        if (JELEM	.eq.	290	)		coeff = 	0.707986177
        if (JELEM	.eq.	291	)		coeff = 	1.02814234
        if (JELEM	.eq.	292	)		coeff = 	1.043711506
        if (JELEM	.eq.	293	)		coeff = 	0.732025748
        if (JELEM	.eq.	294	)		coeff = 	0.725374068
        if (JELEM	.eq.	295	)		coeff = 	1.019968787
        if (JELEM	.eq.	296	)		coeff = 	1.125273925
        if (JELEM	.eq.	297	)		coeff = 	0.839771253
        if (JELEM	.eq.	298	)		coeff = 	1.092841929
        if (JELEM	.eq.	299	)		coeff = 	0.968376825
        if (JELEM	.eq.	300	)		coeff = 	0.975587403
        if (JELEM	.eq.	301	)		coeff = 	0.917581345
        if (JELEM	.eq.	302	)		coeff = 	1.094770538
        if (JELEM	.eq.	303	)		coeff = 	0.919820869
        if (JELEM	.eq.	304	)		coeff = 	1.092135086
        if (JELEM	.eq.	305	)		coeff = 	1.101042848
        if (JELEM	.eq.	306	)		coeff = 	0.975282116
        if (JELEM	.eq.	307	)		coeff = 	1.0338619
        if (JELEM	.eq.	308	)		coeff = 	0.699799197
        if (JELEM	.eq.	309	)		coeff = 	0.949458835
        if (JELEM	.eq.	310	)		coeff = 	0.973857129
        if (JELEM	.eq.	311	)		coeff = 	1.020444224
        if (JELEM	.eq.	312	)		coeff = 	0.965620669
        if (JELEM	.eq.	313	)		coeff = 	0.840549256
        if (JELEM	.eq.	314	)		coeff = 	1.074659566
        if (JELEM	.eq.	315	)		coeff = 	1.181467576
        if (JELEM	.eq.	316	)		coeff = 	0.875214022
        if (JELEM	.eq.	317	)		coeff = 	0.917525553
        if (JELEM	.eq.	318	)		coeff = 	0.984519513
        if (JELEM	.eq.	319	)		coeff = 	0.739135329
        if (JELEM	.eq.	320	)		coeff = 	0.963038577
        if (JELEM	.eq.	321	)		coeff = 	0.737477983
        if (JELEM	.eq.	322	)		coeff = 	0.95628434
        if (JELEM	.eq.	323	)		coeff = 	0.723654021
        if (JELEM	.eq.	324	)		coeff = 	0.877478463
        if (JELEM	.eq.	325	)		coeff = 	0.881436537
        if (JELEM	.eq.	326	)		coeff = 	1.23418129
        if (JELEM	.eq.	327	)		coeff = 	0.941191956
        if (JELEM	.eq.	328	)		coeff = 	0.984758699
        if (JELEM	.eq.	329	)		coeff = 	0.908725224
        if (JELEM	.eq.	330	)		coeff = 	0.94001753
        if (JELEM	.eq.	331	)		coeff = 	1.0548108
        if (JELEM	.eq.	332	)		coeff = 	0.936968281
        if (JELEM	.eq.	333	)		coeff = 	0.787313041
        if (JELEM	.eq.	334	)		coeff = 	1.232089351
        if (JELEM	.eq.	335	)		coeff = 	1.095786047
        if (JELEM	.eq.	336	)		coeff = 	0.842898234
        if (JELEM	.eq.	337	)		coeff = 	1.185179171
        if (JELEM	.eq.	338	)		coeff = 	1.086886463
        if (JELEM	.eq.	339	)		coeff = 	1.01630723
        if (JELEM	.eq.	340	)		coeff = 	0.749984697
        if (JELEM	.eq.	341	)		coeff = 	0.9141082
        if (JELEM	.eq.	342	)		coeff = 	1.194388407
        if (JELEM	.eq.	343	)		coeff = 	1.019458768
        if (JELEM	.eq.	344	)		coeff = 	1.064294647
        if (JELEM	.eq.	345	)		coeff = 	1.032270623
        if (JELEM	.eq.	346	)		coeff = 	1.174570708
        if (JELEM	.eq.	347	)		coeff = 	0.810140111
        if (JELEM	.eq.	348	)		coeff = 	0.848561613
        if (JELEM	.eq.	349	)		coeff = 	0.869564015
        if (JELEM	.eq.	350	)		coeff = 	0.815255987
        if (JELEM	.eq.	351	)		coeff = 	0.848424408
        if (JELEM	.eq.	352	)		coeff = 	1.09819106
        if (JELEM	.eq.	353	)		coeff = 	0.994423563
        if (JELEM	.eq.	354	)		coeff = 	1.172761094
        if (JELEM	.eq.	355	)		coeff = 	1.267030119
        if (JELEM	.eq.	356	)		coeff = 	1.221683097
        if (JELEM	.eq.	357	)		coeff = 	1.075120397
        if (JELEM	.eq.	358	)		coeff = 	1.083125507
        if (JELEM	.eq.	359	)		coeff = 	0.949858372
        if (JELEM	.eq.	360	)		coeff = 	0.982766149
        if (JELEM	.eq.	361	)		coeff = 	0.921774443
        if (JELEM	.eq.	362	)		coeff = 	1.031721436
        if (JELEM	.eq.	363	)		coeff = 	0.922141153
        if (JELEM	.eq.	364	)		coeff = 	0.924987459
        if (JELEM	.eq.	365	)		coeff = 	1.097641675
        if (JELEM	.eq.	366	)		coeff = 	1.123772023
        if (JELEM	.eq.	367	)		coeff = 	1.166995111
        if (JELEM	.eq.	368	)		coeff = 	1.074948122
        if (JELEM	.eq.	369	)		coeff = 	0.94683383
        if (JELEM	.eq.	370	)		coeff = 	0.8029149
        if (JELEM	.eq.	371	)		coeff = 	1.021209968
        if (JELEM	.eq.	372	)		coeff = 	1.064971521
        if (JELEM	.eq.	373	)		coeff = 	1.133837456
        if (JELEM	.eq.	374	)		coeff = 	0.696688582
        if (JELEM	.eq.	375	)		coeff = 	1.043648346
        if (JELEM	.eq.	376	)		coeff = 	0.909783234
        if (JELEM	.eq.	377	)		coeff = 	1.058560556
        if (JELEM	.eq.	378	)		coeff = 	1.044795859
        if (JELEM	.eq.	379	)		coeff = 	0.859021104
        if (JELEM	.eq.	380	)		coeff = 	1.080505684
        if (JELEM	.eq.	381	)		coeff = 	0.814966957
        if (JELEM	.eq.	382	)		coeff = 	1.113576531
        if (JELEM	.eq.	383	)		coeff = 	0.614791813
        if (JELEM	.eq.	384	)		coeff = 	0.898785022
        if (JELEM	.eq.	385	)		coeff = 	1.029144643
        if (JELEM	.eq.	386	)		coeff = 	0.783803104
        if (JELEM	.eq.	387	)		coeff = 	1.169601938
        if (JELEM	.eq.	388	)		coeff = 	0.850277453
        if (JELEM	.eq.	389	)		coeff = 	0.758551243
        if (JELEM	.eq.	390	)		coeff = 	0.912301678
        if (JELEM	.eq.	391	)		coeff = 	0.86975374
        if (JELEM	.eq.	392	)		coeff = 	0.967573281
        if (JELEM	.eq.	393	)		coeff = 	0.851323423
        if (JELEM	.eq.	394	)		coeff = 	1.221055127
        if (JELEM	.eq.	395	)		coeff = 	1.063995653
        if (JELEM	.eq.	396	)		coeff = 	1.115177156
        if (JELEM	.eq.	397	)		coeff = 	1.075053847
        if (JELEM	.eq.	398	)		coeff = 	0.702752564
        if (JELEM	.eq.	399	)		coeff = 	0.797233863
        if (JELEM	.eq.	400	)		coeff = 	1.089530733
        if (JELEM	.eq.	401	)		coeff = 	0.895735293
        if (JELEM	.eq.	402	)		coeff = 	1.15906636
        if (JELEM	.eq.	403	)		coeff = 	1.084209741
        if (JELEM	.eq.	404	)		coeff = 	1.067685321
        if (JELEM	.eq.	405	)		coeff = 	0.901939318
        if (JELEM	.eq.	406	)		coeff = 	0.792017829
        if (JELEM	.eq.	407	)		coeff = 	0.794369612
        if (JELEM	.eq.	408	)		coeff = 	0.957270392
        if (JELEM	.eq.	409	)		coeff = 	0.839771534
        if (JELEM	.eq.	410	)		coeff = 	1.032566565
        if (JELEM	.eq.	411	)		coeff = 	0.653256327
        if (JELEM	.eq.	412	)		coeff = 	0.558863376
        if (JELEM	.eq.	413	)		coeff = 	1.119753491
        if (JELEM	.eq.	414	)		coeff = 	1.013049298
        if (JELEM	.eq.	415	)		coeff = 	1.034886487
        if (JELEM	.eq.	416	)		coeff = 	1.054840753
        if (JELEM	.eq.	417	)		coeff = 	0.847139437
        if (JELEM	.eq.	418	)		coeff = 	1.000960364
        if (JELEM	.eq.	419	)		coeff = 	1.081434555
        if (JELEM	.eq.	420	)		coeff = 	1.076365085
        if (JELEM	.eq.	421	)		coeff = 	0.761237126
        if (JELEM	.eq.	422	)		coeff = 	0.697053205
        if (JELEM	.eq.	423	)		coeff = 	0.949172008
        if (JELEM	.eq.	424	)		coeff = 	1.189090068
        if (JELEM	.eq.	425	)		coeff = 	0.98262348
        if (JELEM	.eq.	426	)		coeff = 	1.097376996
        if (JELEM	.eq.	427	)		coeff = 	0.686532301
        if (JELEM	.eq.	428	)		coeff = 	1.034197742
        if (JELEM	.eq.	429	)		coeff = 	0.760963595
        if (JELEM	.eq.	430	)		coeff = 	0.759498609
        if (JELEM	.eq.	431	)		coeff = 	0.84265191
        if (JELEM	.eq.	432	)		coeff = 	0.745695164
        if (JELEM	.eq.	433	)		coeff = 	1.192647616
        if (JELEM	.eq.	434	)		coeff = 	0.671952539
        if (JELEM	.eq.	435	)		coeff = 	1.126153386
        if (JELEM	.eq.	436	)		coeff = 	1.028778131
        if (JELEM	.eq.	437	)		coeff = 	0.959301283
        if (JELEM	.eq.	438	)		coeff = 	0.644769279
        if (JELEM	.eq.	439	)		coeff = 	0.931210322
        if (JELEM	.eq.	440	)		coeff = 	1.032815219
        if (JELEM	.eq.	441	)		coeff = 	1.093698466
        if (JELEM	.eq.	442	)		coeff = 	1.080884484
        if (JELEM	.eq.	443	)		coeff = 	1.034698523
        if (JELEM	.eq.	444	)		coeff = 	0.977706647
        if (JELEM	.eq.	445	)		coeff = 	0.990514075
        if (JELEM	.eq.	446	)		coeff = 	0.980621122
        if (JELEM	.eq.	447	)		coeff = 	0.914019006
        if (JELEM	.eq.	448	)		coeff = 	0.973345661
        if (JELEM	.eq.	449	)		coeff = 	1.204173059
        if (JELEM	.eq.	450	)		coeff = 	1.164151298
        if (JELEM	.eq.	451	)		coeff = 	0.989546604
        if (JELEM	.eq.	452	)		coeff = 	0.616117996
        if (JELEM	.eq.	453	)		coeff = 	0.940289345
        if (JELEM	.eq.	454	)		coeff = 	0.741806835
        if (JELEM	.eq.	455	)		coeff = 	1.015003209
        if (JELEM	.eq.	456	)		coeff = 	1.026848738
        if (JELEM	.eq.	457	)		coeff = 	0.774147598
        if (JELEM	.eq.	458	)		coeff = 	0.840710918
        if (JELEM	.eq.	459	)		coeff = 	0.991017487
        if (JELEM	.eq.	460	)		coeff = 	0.771254686
        if (JELEM	.eq.	461	)		coeff = 	0.980478869
        if (JELEM	.eq.	462	)		coeff = 	0.966714888
        if (JELEM	.eq.	463	)		coeff = 	0.980727412
        if (JELEM	.eq.	464	)		coeff = 	0.970502344
        if (JELEM	.eq.	465	)		coeff = 	0.881095786
        if (JELEM	.eq.	466	)		coeff = 	1.112584163
        if (JELEM	.eq.	467	)		coeff = 	0.863807602
        if (JELEM	.eq.	468	)		coeff = 	1.093216207
        if (JELEM	.eq.	469	)		coeff = 	0.95641437
        if (JELEM	.eq.	470	)		coeff = 	1.038070528
        if (JELEM	.eq.	471	)		coeff = 	0.922572392
        if (JELEM	.eq.	472	)		coeff = 	1.050233479
        if (JELEM	.eq.	473	)		coeff = 	0.680379232
        if (JELEM	.eq.	474	)		coeff = 	1.135537498
        if (JELEM	.eq.	475	)		coeff = 	1.068335483
        if (JELEM	.eq.	476	)		coeff = 	0.910233133
        if (JELEM	.eq.	477	)		coeff = 	1.026531032
        if (JELEM	.eq.	478	)		coeff = 	1.089578528
        if (JELEM	.eq.	479	)		coeff = 	1.065911422
        if (JELEM	.eq.	480	)		coeff = 	0.677429807
        if (JELEM	.eq.	481	)		coeff = 	1.024676156
        if (JELEM	.eq.	482	)		coeff = 	1.046529666
        if (JELEM	.eq.	483	)		coeff = 	0.8959674
        if (JELEM	.eq.	484	)		coeff = 	1.236090713
        if (JELEM	.eq.	485	)		coeff = 	1.124071312
        if (JELEM	.eq.	486	)		coeff = 	1.020094127
        if (JELEM	.eq.	487	)		coeff = 	1.051659521
        if (JELEM	.eq.	488	)		coeff = 	0.881154709
        if (JELEM	.eq.	489	)		coeff = 	0.878989402
        if (JELEM	.eq.	490	)		coeff = 	0.795510037
        if (JELEM	.eq.	491	)		coeff = 	0.934477051
        if (JELEM	.eq.	492	)		coeff = 	0.930405344
        if (JELEM	.eq.	493	)		coeff = 	0.944469098
        if (JELEM	.eq.	494	)		coeff = 	0.790265126
        if (JELEM	.eq.	495	)		coeff = 	0.973500987
        if (JELEM	.eq.	496	)		coeff = 	0.994745031
        if (JELEM	.eq.	497	)		coeff = 	0.921198319
        if (JELEM	.eq.	498	)		coeff = 	1.021693494
        if (JELEM	.eq.	499	)		coeff = 	1.109217208
        if (JELEM	.eq.	500	)		coeff = 	0.99862142
        if (JELEM	.eq.	501	)		coeff = 	1.052372864
        if (JELEM	.eq.	502	)		coeff = 	0.902789697
        if (JELEM	.eq.	503	)		coeff = 	0.842050229
        if (JELEM	.eq.	504	)		coeff = 	0.979006977
        if (JELEM	.eq.	505	)		coeff = 	0.911460676
        if (JELEM	.eq.	506	)		coeff = 	1.030962984
        if (JELEM	.eq.	507	)		coeff = 	0.962742902
        if (JELEM	.eq.	508	)		coeff = 	1.018286264
        if (JELEM	.eq.	509	)		coeff = 	0.812402662
        if (JELEM	.eq.	510	)		coeff = 	0.859143497
        if (JELEM	.eq.	511	)		coeff = 	1.086068184
        if (JELEM	.eq.	512	)		coeff = 	1.131370963
        if (JELEM	.eq.	513	)		coeff = 	0.944041766
        if (JELEM	.eq.	514	)		coeff = 	0.858756552
        if (JELEM	.eq.	515	)		coeff = 	0.806754743
        if (JELEM	.eq.	516	)		coeff = 	0.869791746
        if (JELEM	.eq.	517	)		coeff = 	0.76654641
        if (JELEM	.eq.	518	)		coeff = 	1.077349805
        if (JELEM	.eq.	519	)		coeff = 	1.06078024
        if (JELEM	.eq.	520	)		coeff = 	0.771572358
        if (JELEM	.eq.	521	)		coeff = 	0.96821498
        if (JELEM	.eq.	522	)		coeff = 	1.066710489
        if (JELEM	.eq.	523	)		coeff = 	0.867838977
        if (JELEM	.eq.	524	)		coeff = 	0.770179606
        if (JELEM	.eq.	525	)		coeff = 	0.830185659
        if (JELEM	.eq.	526	)		coeff = 	0.801803054
        if (JELEM	.eq.	527	)		coeff = 	0.941959831
        if (JELEM	.eq.	528	)		coeff = 	1.057380268
        if (JELEM	.eq.	529	)		coeff = 	0.944635849
        if (JELEM	.eq.	530	)		coeff = 	0.801573132
        if (JELEM	.eq.	531	)		coeff = 	1.055072951
        if (JELEM	.eq.	532	)		coeff = 	0.871596762
        if (JELEM	.eq.	533	)		coeff = 	0.787375232
        if (JELEM	.eq.	534	)		coeff = 	0.965053722
        if (JELEM	.eq.	535	)		coeff = 	0.968708437
        if (JELEM	.eq.	536	)		coeff = 	0.702312688
        if (JELEM	.eq.	537	)		coeff = 	0.908642505
        if (JELEM	.eq.	538	)		coeff = 	0.825283683
        if (JELEM	.eq.	539	)		coeff = 	0.988612261
        if (JELEM	.eq.	540	)		coeff = 	1.015964694
        if (JELEM	.eq.	541	)		coeff = 	0.996293637
        if (JELEM	.eq.	542	)		coeff = 	0.87057015
        if (JELEM	.eq.	543	)		coeff = 	1.064172265
        if (JELEM	.eq.	544	)		coeff = 	0.958660598
        if (JELEM	.eq.	545	)		coeff = 	1.09808143
        if (JELEM	.eq.	546	)		coeff = 	1.088220386
        if (JELEM	.eq.	547	)		coeff = 	0.6947041
        if (JELEM	.eq.	548	)		coeff = 	0.907628905
        if (JELEM	.eq.	549	)		coeff = 	0.97609833
        if (JELEM	.eq.	550	)		coeff = 	0.875647366
        if (JELEM	.eq.	551	)		coeff = 	0.944916317
        if (JELEM	.eq.	552	)		coeff = 	0.899105431
        if (JELEM	.eq.	553	)		coeff = 	0.841228919
        if (JELEM	.eq.	554	)		coeff = 	0.912630497
        if (JELEM	.eq.	555	)		coeff = 	0.954499966
        if (JELEM	.eq.	556	)		coeff = 	0.960240844
        if (JELEM	.eq.	557	)		coeff = 	1.103947975
        if (JELEM	.eq.	558	)		coeff = 	0.737100075
        if (JELEM	.eq.	559	)		coeff = 	0.978135328
        if (JELEM	.eq.	560	)		coeff = 	0.980030569
        if (JELEM	.eq.	561	)		coeff = 	1.170797259
        if (JELEM	.eq.	562	)		coeff = 	0.951617138
        if (JELEM	.eq.	563	)		coeff = 	0.662192867
        if (JELEM	.eq.	564	)		coeff = 	1.032313391
        if (JELEM	.eq.	565	)		coeff = 	0.859021668
        if (JELEM	.eq.	566	)		coeff = 	1.100235977
        if (JELEM	.eq.	567	)		coeff = 	0.845933549
        if (JELEM	.eq.	568	)		coeff = 	0.847815943
        if (JELEM	.eq.	569	)		coeff = 	0.940998476
        if (JELEM	.eq.	570	)		coeff = 	1.074433677
        if (JELEM	.eq.	571	)		coeff = 	0.720274723
        if (JELEM	.eq.	572	)		coeff = 	1.057938238
        if (JELEM	.eq.	573	)		coeff = 	1.186361364
        if (JELEM	.eq.	574	)		coeff = 	0.918857872
        if (JELEM	.eq.	575	)		coeff = 	0.916686228
        if (JELEM	.eq.	576	)		coeff = 	0.965873169
        if (JELEM	.eq.	577	)		coeff = 	0.920842004
        if (JELEM	.eq.	578	)		coeff = 	0.868499989
        if (JELEM	.eq.	579	)		coeff = 	0.986396775
        if (JELEM	.eq.	580	)		coeff = 	0.920356155
        if (JELEM	.eq.	581	)		coeff = 	0.980464996
        if (JELEM	.eq.	582	)		coeff = 	1.097507187
        if (JELEM	.eq.	583	)		coeff = 	1.006112074
        if (JELEM	.eq.	584	)		coeff = 	0.938877843
        if (JELEM	.eq.	585	)		coeff = 	1.094729505
        if (JELEM	.eq.	586	)		coeff = 	1.121756247
        if (JELEM	.eq.	587	)		coeff = 	1.027227586
        if (JELEM	.eq.	588	)		coeff = 	0.889063001
        if (JELEM	.eq.	589	)		coeff = 	0.969881778
        if (JELEM	.eq.	590	)		coeff = 	0.955489872
        if (JELEM	.eq.	591	)		coeff = 	1.139324917
        if (JELEM	.eq.	592	)		coeff = 	1.12735273
        if (JELEM	.eq.	593	)		coeff = 	1.010258283
        if (JELEM	.eq.	594	)		coeff = 	1.010880554
        if (JELEM	.eq.	595	)		coeff = 	0.592335845
        if (JELEM	.eq.	596	)		coeff = 	0.85495704
        if (JELEM	.eq.	597	)		coeff = 	0.901472547
        if (JELEM	.eq.	598	)		coeff = 	0.770783074
        if (JELEM	.eq.	599	)		coeff = 	0.951109745
        if (JELEM	.eq.	600	)		coeff = 	0.997175597
        if (JELEM	.eq.	601	)		coeff = 	0.769995085
        if (JELEM	.eq.	602	)		coeff = 	1.043271613
        if (JELEM	.eq.	603	)		coeff = 	1.076826825
        if (JELEM	.eq.	604	)		coeff = 	0.896601289
        if (JELEM	.eq.	605	)		coeff = 	1.039572042
        if (JELEM	.eq.	606	)		coeff = 	1.079762365
        if (JELEM	.eq.	607	)		coeff = 	1.031857518
        if (JELEM	.eq.	608	)		coeff = 	1.023932269
        if (JELEM	.eq.	609	)		coeff = 	1.226836162
        if (JELEM	.eq.	610	)		coeff = 	1.049830999
        if (JELEM	.eq.	611	)		coeff = 	0.8918636
        if (JELEM	.eq.	612	)		coeff = 	0.607309161
        if (JELEM	.eq.	613	)		coeff = 	1.039374391
        if (JELEM	.eq.	614	)		coeff = 	0.91218722
        if (JELEM	.eq.	615	)		coeff = 	0.641206799
        if (JELEM	.eq.	616	)		coeff = 	0.957911324
        if (JELEM	.eq.	617	)		coeff = 	1.023111063
        if (JELEM	.eq.	618	)		coeff = 	0.920549376
        if (JELEM	.eq.	619	)		coeff = 	0.940071739
        if (JELEM	.eq.	620	)		coeff = 	1.0933434
        if (JELEM	.eq.	621	)		coeff = 	0.936927375
        if (JELEM	.eq.	622	)		coeff = 	1.083778769
        if (JELEM	.eq.	623	)		coeff = 	1.038482069
        if (JELEM	.eq.	624	)		coeff = 	0.945244108
        if (JELEM	.eq.	625	)		coeff = 	0.821071586
        if (JELEM	.eq.	626	)		coeff = 	0.870769718
        if (JELEM	.eq.	627	)		coeff = 	0.954937174
        if (JELEM	.eq.	628	)		coeff = 	0.760652615
        if (JELEM	.eq.	629	)		coeff = 	0.688578568
        if (JELEM	.eq.	630	)		coeff = 	0.968080478
        if (JELEM	.eq.	631	)		coeff = 	1.027893692
        if (JELEM	.eq.	632	)		coeff = 	0.725181662
        if (JELEM	.eq.	633	)		coeff = 	0.956483325
        if (JELEM	.eq.	634	)		coeff = 	0.775542943
        if (JELEM	.eq.	635	)		coeff = 	1.055913255
        if (JELEM	.eq.	636	)		coeff = 	1.101560836
        if (JELEM	.eq.	637	)		coeff = 	1.096191668
        if (JELEM	.eq.	638	)		coeff = 	0.929033518
        if (JELEM	.eq.	639	)		coeff = 	0.966716842
        if (JELEM	.eq.	640	)		coeff = 	0.997909292
        if (JELEM	.eq.	641	)		coeff = 	0.853481158
        if (JELEM	.eq.	642	)		coeff = 	0.970369574
        if (JELEM	.eq.	643	)		coeff = 	0.876421433
        if (JELEM	.eq.	644	)		coeff = 	0.721090266
        if (JELEM	.eq.	645	)		coeff = 	0.816678471
        if (JELEM	.eq.	646	)		coeff = 	1.157517087
        if (JELEM	.eq.	647	)		coeff = 	0.656241643
        if (JELEM	.eq.	648	)		coeff = 	0.766233022
        if (JELEM	.eq.	649	)		coeff = 	0.904554157
        if (JELEM	.eq.	650	)		coeff = 	1.207315272
        if (JELEM	.eq.	651	)		coeff = 	0.624058253
        if (JELEM	.eq.	652	)		coeff = 	0.994018869
        if (JELEM	.eq.	653	)		coeff = 	0.98924689
        if (JELEM	.eq.	654	)		coeff = 	0.95166505
        if (JELEM	.eq.	655	)		coeff = 	0.9371377
        if (JELEM	.eq.	656	)		coeff = 	0.672781824
        if (JELEM	.eq.	657	)		coeff = 	0.792874968
        if (JELEM	.eq.	658	)		coeff = 	0.963200643
        if (JELEM	.eq.	659	)		coeff = 	0.775132385
        if (JELEM	.eq.	660	)		coeff = 	0.992548772
        if (JELEM	.eq.	661	)		coeff = 	1.147805002
        if (JELEM	.eq.	662	)		coeff = 	1.02084989
        if (JELEM	.eq.	663	)		coeff = 	1.016377271
        if (JELEM	.eq.	664	)		coeff = 	0.968439994
        if (JELEM	.eq.	665	)		coeff = 	0.959295648
        if (JELEM	.eq.	666	)		coeff = 	0.876103348
        if (JELEM	.eq.	667	)		coeff = 	1.000197912
        if (JELEM	.eq.	668	)		coeff = 	0.704012597
        if (JELEM	.eq.	669	)		coeff = 	1.034927258
        if (JELEM	.eq.	670	)		coeff = 	0.858052602
        if (JELEM	.eq.	671	)		coeff = 	0.745424107
        if (JELEM	.eq.	672	)		coeff = 	0.925204094
        if (JELEM	.eq.	673	)		coeff = 	0.748155988
        if (JELEM	.eq.	674	)		coeff = 	1.140881471
        if (JELEM	.eq.	675	)		coeff = 	0.935310897
        if (JELEM	.eq.	676	)		coeff = 	1.246335649
        if (JELEM	.eq.	677	)		coeff = 	1.081839131
        if (JELEM	.eq.	678	)		coeff = 	0.98660557
        if (JELEM	.eq.	679	)		coeff = 	1.047974222
        if (JELEM	.eq.	680	)		coeff = 	0.942923605
        if (JELEM	.eq.	681	)		coeff = 	1.019384493
        if (JELEM	.eq.	682	)		coeff = 	0.875256497
        if (JELEM	.eq.	683	)		coeff = 	0.924278669
        if (JELEM	.eq.	684	)		coeff = 	1.045768785
        if (JELEM	.eq.	685	)		coeff = 	1.134687165
        if (JELEM	.eq.	686	)		coeff = 	1.242849177
        if (JELEM	.eq.	687	)		coeff = 	1.065852446
        if (JELEM	.eq.	688	)		coeff = 	0.958243294
        if (JELEM	.eq.	689	)		coeff = 	0.913224049
        if (JELEM	.eq.	690	)		coeff = 	1.0970085
        if (JELEM	.eq.	691	)		coeff = 	1.091673547
        if (JELEM	.eq.	692	)		coeff = 	1.018482481
        if (JELEM	.eq.	693	)		coeff = 	1.117498415
        if (JELEM	.eq.	694	)		coeff = 	0.986166641
        if (JELEM	.eq.	695	)		coeff = 	1.029791211
        if (JELEM	.eq.	696	)		coeff = 	0.907766699
        if (JELEM	.eq.	697	)		coeff = 	1.05331964
        if (JELEM	.eq.	698	)		coeff = 	0.898367454
        if (JELEM	.eq.	699	)		coeff = 	1.064002951
        if (JELEM	.eq.	700	)		coeff = 	1.113192356
        if (JELEM	.eq.	701	)		coeff = 	0.734992776
        if (JELEM	.eq.	702	)		coeff = 	0.934771578
        if (JELEM	.eq.	703	)		coeff = 	0.983184479
        if (JELEM	.eq.	704	)		coeff = 	1.103779065
        if (JELEM	.eq.	705	)		coeff = 	0.936739812
        if (JELEM	.eq.	706	)		coeff = 	1.103307677
        if (JELEM	.eq.	707	)		coeff = 	1.046982553
        if (JELEM	.eq.	708	)		coeff = 	1.088847667
        if (JELEM	.eq.	709	)		coeff = 	1.013587271
        if (JELEM	.eq.	710	)		coeff = 	1.007314904
        if (JELEM	.eq.	711	)		coeff = 	0.722082441
        if (JELEM	.eq.	712	)		coeff = 	0.765916356
        if (JELEM	.eq.	713	)		coeff = 	0.834712893
        if (JELEM	.eq.	714	)		coeff = 	1.126721035
        if (JELEM	.eq.	715	)		coeff = 	0.98245784
        if (JELEM	.eq.	716	)		coeff = 	0.990237024
        if (JELEM	.eq.	717	)		coeff = 	0.984772096
        if (JELEM	.eq.	718	)		coeff = 	0.634866036
        if (JELEM	.eq.	719	)		coeff = 	1.046264469
        if (JELEM	.eq.	720	)		coeff = 	1.104845343
        if (JELEM	.eq.	721	)		coeff = 	0.82901844
        if (JELEM	.eq.	722	)		coeff = 	1.152158845
        if (JELEM	.eq.	723	)		coeff = 	0.670536247
        if (JELEM	.eq.	724	)		coeff = 	0.828077947
        if (JELEM	.eq.	725	)		coeff = 	0.669160395
        if (JELEM	.eq.	726	)		coeff = 	0.847800814
        if (JELEM	.eq.	727	)		coeff = 	0.916360278
        if (JELEM	.eq.	728	)		coeff = 	0.569203091
        if (JELEM	.eq.	729	)		coeff = 	0.933218714
        if (JELEM	.eq.	730	)		coeff = 	0.956634663
        if (JELEM	.eq.	731	)		coeff = 	1.040534169
        if (JELEM	.eq.	732	)		coeff = 	1.062962818
        if (JELEM	.eq.	733	)		coeff = 	0.957717997
        if (JELEM	.eq.	734	)		coeff = 	0.874245461
        if (JELEM	.eq.	735	)		coeff = 	0.578007746
        if (JELEM	.eq.	736	)		coeff = 	1.01138126
        if (JELEM	.eq.	737	)		coeff = 	1.05315766
        if (JELEM	.eq.	738	)		coeff = 	1.007713335
        if (JELEM	.eq.	739	)		coeff = 	1.115078761
        if (JELEM	.eq.	740	)		coeff = 	0.907371995
        if (JELEM	.eq.	741	)		coeff = 	0.996502906
        if (JELEM	.eq.	742	)		coeff = 	0.766781339
        if (JELEM	.eq.	743	)		coeff = 	0.980189294
        if (JELEM	.eq.	744	)		coeff = 	0.763659963
        if (JELEM	.eq.	745	)		coeff = 	1.051558
        if (JELEM	.eq.	746	)		coeff = 	0.830807862
        if (JELEM	.eq.	747	)		coeff = 	0.790721864
        if (JELEM	.eq.	748	)		coeff = 	0.979997529
        if (JELEM	.eq.	749	)		coeff = 	0.923613242
        if (JELEM	.eq.	750	)		coeff = 	1.115019436
        if (JELEM	.eq.	751	)		coeff = 	1.041455355
        if (JELEM	.eq.	752	)		coeff = 	0.921789784
        if (JELEM	.eq.	753	)		coeff = 	0.948087875
        if (JELEM	.eq.	754	)		coeff = 	0.739944556
        if (JELEM	.eq.	755	)		coeff = 	0.894683198
        if (JELEM	.eq.	756	)		coeff = 	0.955339201
        if (JELEM	.eq.	757	)		coeff = 	0.906117886
        if (JELEM	.eq.	758	)		coeff = 	1.118480935
        if (JELEM	.eq.	759	)		coeff = 	1.156445125
        if (JELEM	.eq.	760	)		coeff = 	0.855006912
        if (JELEM	.eq.	761	)		coeff = 	1.015143005
        if (JELEM	.eq.	762	)		coeff = 	1.132752582
        if (JELEM	.eq.	763	)		coeff = 	0.729092839
        if (JELEM	.eq.	764	)		coeff = 	1.142947007
        if (JELEM	.eq.	765	)		coeff = 	1.131620256
        if (JELEM	.eq.	766	)		coeff = 	0.929240723
        if (JELEM	.eq.	767	)		coeff = 	0.975915072
        if (JELEM	.eq.	768	)		coeff = 	1.173580753
        if (JELEM	.eq.	769	)		coeff = 	1.147287005
        if (JELEM	.eq.	770	)		coeff = 	1.046012804
        if (JELEM	.eq.	771	)		coeff = 	1.09462403
        if (JELEM	.eq.	772	)		coeff = 	0.907868757
        if (JELEM	.eq.	773	)		coeff = 	0.963299932
        if (JELEM	.eq.	774	)		coeff = 	1.104034591
        if (JELEM	.eq.	775	)		coeff = 	0.977186332
        if (JELEM	.eq.	776	)		coeff = 	0.752921404
        if (JELEM	.eq.	777	)		coeff = 	0.932035336
        if (JELEM	.eq.	778	)		coeff = 	0.865140643
        if (JELEM	.eq.	779	)		coeff = 	0.610209733
        if (JELEM	.eq.	780	)		coeff = 	0.923854103
        if (JELEM	.eq.	781	)		coeff = 	0.954129403
        if (JELEM	.eq.	782	)		coeff = 	0.944232547
        if (JELEM	.eq.	783	)		coeff = 	0.917096045
        if (JELEM	.eq.	784	)		coeff = 	0.848541549
        if (JELEM	.eq.	785	)		coeff = 	0.919228561
        if (JELEM	.eq.	786	)		coeff = 	0.868904209
        if (JELEM	.eq.	787	)		coeff = 	0.977259465
        if (JELEM	.eq.	788	)		coeff = 	0.901475101
        if (JELEM	.eq.	789	)		coeff = 	1.071174151
        if (JELEM	.eq.	790	)		coeff = 	0.776261958
        if (JELEM	.eq.	791	)		coeff = 	1.041040063
        if (JELEM	.eq.	792	)		coeff = 	1.003880082
        if (JELEM	.eq.	793	)		coeff = 	1.022291342
        if (JELEM	.eq.	794	)		coeff = 	0.652284551
        if (JELEM	.eq.	795	)		coeff = 	1.064496173
        if (JELEM	.eq.	796	)		coeff = 	1.111177753
        if (JELEM	.eq.	797	)		coeff = 	0.881807327
        if (JELEM	.eq.	798	)		coeff = 	1.112308
        if (JELEM	.eq.	799	)		coeff = 	1.061543577
        if (JELEM	.eq.	800	)		coeff = 	0.916805805
        if (JELEM	.eq.	801	)		coeff = 	0.874773341
        if (JELEM	.eq.	802	)		coeff = 	0.763888129
        if (JELEM	.eq.	803	)		coeff = 	0.918470524
        if (JELEM	.eq.	804	)		coeff = 	0.78754112
        if (JELEM	.eq.	805	)		coeff = 	0.889424058
        if (JELEM	.eq.	806	)		coeff = 	0.837516551
        if (JELEM	.eq.	807	)		coeff = 	0.942453202
        if (JELEM	.eq.	808	)		coeff = 	1.235105607
        if (JELEM	.eq.	809	)		coeff = 	1.061863817
        if (JELEM	.eq.	810	)		coeff = 	0.976075037
        if (JELEM	.eq.	811	)		coeff = 	0.93475303
        if (JELEM	.eq.	812	)		coeff = 	0.994153718
        if (JELEM	.eq.	813	)		coeff = 	0.901341596
        if (JELEM	.eq.	814	)		coeff = 	0.924790559
        if (JELEM	.eq.	815	)		coeff = 	0.84749859
        if (JELEM	.eq.	816	)		coeff = 	1.072077691
        if (JELEM	.eq.	817	)		coeff = 	0.954145543
        if (JELEM	.eq.	818	)		coeff = 	0.994983989
        if (JELEM	.eq.	819	)		coeff = 	0.624267715
        if (JELEM	.eq.	820	)		coeff = 	0.712007375
        if (JELEM	.eq.	821	)		coeff = 	0.761828001
        if (JELEM	.eq.	822	)		coeff = 	0.97445254
        if (JELEM	.eq.	823	)		coeff = 	0.855025424
        if (JELEM	.eq.	824	)		coeff = 	0.745068942
        if (JELEM	.eq.	825	)		coeff = 	0.998153755
        if (JELEM	.eq.	826	)		coeff = 	1.100517056
        if (JELEM	.eq.	827	)		coeff = 	0.946501702
        if (JELEM	.eq.	828	)		coeff = 	1.009140415
        if (JELEM	.eq.	829	)		coeff = 	0.797254993
        if (JELEM	.eq.	830	)		coeff = 	1.067058817
        if (JELEM	.eq.	831	)		coeff = 	1.057721223
        if (JELEM	.eq.	832	)		coeff = 	0.944640951
        if (JELEM	.eq.	833	)		coeff = 	1.119863424
        if (JELEM	.eq.	834	)		coeff = 	1.045166556
        if (JELEM	.eq.	835	)		coeff = 	0.90621603
        if (JELEM	.eq.	836	)		coeff = 	1.103423272
        if (JELEM	.eq.	837	)		coeff = 	0.916696213
        if (JELEM	.eq.	838	)		coeff = 	0.764598157
        if (JELEM	.eq.	839	)		coeff = 	0.958420963
        if (JELEM	.eq.	840	)		coeff = 	1.018718123
        if (JELEM	.eq.	841	)		coeff = 	0.84966185
        if (JELEM	.eq.	842	)		coeff = 	0.681519541
        if (JELEM	.eq.	843	)		coeff = 	0.901731071
        if (JELEM	.eq.	844	)		coeff = 	0.960258946
        if (JELEM	.eq.	845	)		coeff = 	0.976258257
        if (JELEM	.eq.	846	)		coeff = 	1.066994692
        if (JELEM	.eq.	847	)		coeff = 	0.974586088
        if (JELEM	.eq.	848	)		coeff = 	0.95930037
        if (JELEM	.eq.	849	)		coeff = 	0.788557
        if (JELEM	.eq.	850	)		coeff = 	0.871719954
        if (JELEM	.eq.	851	)		coeff = 	0.911741243
        if (JELEM	.eq.	852	)		coeff = 	1.116252849
        if (JELEM	.eq.	853	)		coeff = 	1.103036592
        if (JELEM	.eq.	854	)		coeff = 	0.928491363
        if (JELEM	.eq.	855	)		coeff = 	0.968317246
        if (JELEM	.eq.	856	)		coeff = 	0.939338446
        if (JELEM	.eq.	857	)		coeff = 	1.159620601
        if (JELEM	.eq.	858	)		coeff = 	0.927139255
        if (JELEM	.eq.	859	)		coeff = 	0.908395261
        if (JELEM	.eq.	860	)		coeff = 	0.9614186
        if (JELEM	.eq.	861	)		coeff = 	1.078413748
        if (JELEM	.eq.	862	)		coeff = 	0.853599763
        if (JELEM	.eq.	863	)		coeff = 	0.8327671
        if (JELEM	.eq.	864	)		coeff = 	1.026544074
        if (JELEM	.eq.	865	)		coeff = 	0.952592877
        if (JELEM	.eq.	866	)		coeff = 	1.046227229
        if (JELEM	.eq.	867	)		coeff = 	1.064331283
        if (JELEM	.eq.	868	)		coeff = 	0.823904054
        if (JELEM	.eq.	869	)		coeff = 	1.132609236
        if (JELEM	.eq.	870	)		coeff = 	1.063003903
        if (JELEM	.eq.	871	)		coeff = 	0.632274866
        if (JELEM	.eq.	872	)		coeff = 	0.92867342
        if (JELEM	.eq.	873	)		coeff = 	0.925697125
        if (JELEM	.eq.	874	)		coeff = 	0.689816295
        if (JELEM	.eq.	875	)		coeff = 	1.109858772
        if (JELEM	.eq.	876	)		coeff = 	0.774352289
        if (JELEM	.eq.	877	)		coeff = 	0.763031519
        if (JELEM	.eq.	878	)		coeff = 	1.055593542
        if (JELEM	.eq.	879	)		coeff = 	1.075494075
        if (JELEM	.eq.	880	)		coeff = 	0.901634556
        if (JELEM	.eq.	881	)		coeff = 	1.069587297
        if (JELEM	.eq.	882	)		coeff = 	1.100374551
        if (JELEM	.eq.	883	)		coeff = 	1.000865268
        if (JELEM	.eq.	884	)		coeff = 	0.988330681
        if (JELEM	.eq.	885	)		coeff = 	1.116229475
        if (JELEM	.eq.	886	)		coeff = 	0.886119809
        if (JELEM	.eq.	887	)		coeff = 	0.938824996
        if (JELEM	.eq.	888	)		coeff = 	0.587231085
        if (JELEM	.eq.	889	)		coeff = 	1.056108958
        if (JELEM	.eq.	890	)		coeff = 	0.734249891
        if (JELEM	.eq.	891	)		coeff = 	0.893914415
        if (JELEM	.eq.	892	)		coeff = 	0.858989105
        if (JELEM	.eq.	893	)		coeff = 	0.996411723
        if (JELEM	.eq.	894	)		coeff = 	0.963860784
        if (JELEM	.eq.	895	)		coeff = 	1.106180054
        if (JELEM	.eq.	896	)		coeff = 	0.76997571
        if (JELEM	.eq.	897	)		coeff = 	0.973920907
        if (JELEM	.eq.	898	)		coeff = 	0.89264754
        if (JELEM	.eq.	899	)		coeff = 	0.866282282
        if (JELEM	.eq.	900	)		coeff = 	1.028594933
        if (JELEM	.eq.	901	)		coeff = 	0.9218213
        if (JELEM	.eq.	902	)		coeff = 	0.999099694
        if (JELEM	.eq.	903	)		coeff = 	0.831547849
        if (JELEM	.eq.	904	)		coeff = 	1.127405317
        if (JELEM	.eq.	905	)		coeff = 	0.540545817
        if (JELEM	.eq.	906	)		coeff = 	0.842386005
        if (JELEM	.eq.	907	)		coeff = 	1.220079788
        if (JELEM	.eq.	908	)		coeff = 	0.765340099
        if (JELEM	.eq.	909	)		coeff = 	0.786245527
        if (JELEM	.eq.	910	)		coeff = 	0.880147907
        if (JELEM	.eq.	911	)		coeff = 	0.768906435
        if (JELEM	.eq.	912	)		coeff = 	0.966964
        if (JELEM	.eq.	913	)		coeff = 	0.998750687
        if (JELEM	.eq.	914	)		coeff = 	1.120716005
        if (JELEM	.eq.	915	)		coeff = 	0.479302123
        if (JELEM	.eq.	916	)		coeff = 	0.781727326
        if (JELEM	.eq.	917	)		coeff = 	1.326099053
        if (JELEM	.eq.	918	)		coeff = 	0.887729684
        if (JELEM	.eq.	919	)		coeff = 	0.898113945
        if (JELEM	.eq.	920	)		coeff = 	0.985148779
        if (JELEM	.eq.	921	)		coeff = 	0.963315544
        if (JELEM	.eq.	922	)		coeff = 	0.993091717
        if (JELEM	.eq.	923	)		coeff = 	0.971569632
        if (JELEM	.eq.	924	)		coeff = 	1.174495514
        if (JELEM	.eq.	925	)		coeff = 	0.980607673
        if (JELEM	.eq.	926	)		coeff = 	1.106445328
        if (JELEM	.eq.	927	)		coeff = 	1.153944219
        if (JELEM	.eq.	928	)		coeff = 	0.935531101
        if (JELEM	.eq.	929	)		coeff = 	0.915200919
        if (JELEM	.eq.	930	)		coeff = 	0.861239123
        if (JELEM	.eq.	931	)		coeff = 	1.03908761
        if (JELEM	.eq.	932	)		coeff = 	0.867197601
        if (JELEM	.eq.	933	)		coeff = 	0.873175627
        if (JELEM	.eq.	934	)		coeff = 	0.953921243
        if (JELEM	.eq.	935	)		coeff = 	0.76179732
        if (JELEM	.eq.	936	)		coeff = 	0.840580609
        if (JELEM	.eq.	937	)		coeff = 	0.956987959
        if (JELEM	.eq.	938	)		coeff = 	0.84998642
        if (JELEM	.eq.	939	)		coeff = 	0.893137311
        if (JELEM	.eq.	940	)		coeff = 	0.621047528
        if (JELEM	.eq.	941	)		coeff = 	1.100088957
        if (JELEM	.eq.	942	)		coeff = 	0.782134792
        if (JELEM	.eq.	943	)		coeff = 	0.931303826
        if (JELEM	.eq.	944	)		coeff = 	1.098415031
        if (JELEM	.eq.	945	)		coeff = 	0.645440562
        if (JELEM	.eq.	946	)		coeff = 	0.803879885
        if (JELEM	.eq.	947	)		coeff = 	1.136104228
        if (JELEM	.eq.	948	)		coeff = 	1.031325767
        if (JELEM	.eq.	949	)		coeff = 	0.985884614
        if (JELEM	.eq.	950	)		coeff = 	0.909238538
        if (JELEM	.eq.	951	)		coeff = 	0.734116743
        if (JELEM	.eq.	952	)		coeff = 	0.984023529
        if (JELEM	.eq.	953	)		coeff = 	0.873967829
        if (JELEM	.eq.	954	)		coeff = 	0.643923133
        if (JELEM	.eq.	955	)		coeff = 	0.977172727
        if (JELEM	.eq.	956	)		coeff = 	0.905890343
        if (JELEM	.eq.	957	)		coeff = 	0.753247805
        if (JELEM	.eq.	958	)		coeff = 	0.780306963
        if (JELEM	.eq.	959	)		coeff = 	1.0868458
        if (JELEM	.eq.	960	)		coeff = 	1.061044399
        if (JELEM	.eq.	961	)		coeff = 	0.911512097
        if (JELEM	.eq.	962	)		coeff = 	0.842920386
        if (JELEM	.eq.	963	)		coeff = 	1.042875225
        if (JELEM	.eq.	964	)		coeff = 	0.999465737
        if (JELEM	.eq.	965	)		coeff = 	1.009149933
        if (JELEM	.eq.	966	)		coeff = 	0.970102888
        if (JELEM	.eq.	967	)		coeff = 	1.073762799
        if (JELEM	.eq.	968	)		coeff = 	0.857949535
        if (JELEM	.eq.	969	)		coeff = 	1.11110847
        if (JELEM	.eq.	970	)		coeff = 	0.764126598
        if (JELEM	.eq.	971	)		coeff = 	1.12134713
        if (JELEM	.eq.	972	)		coeff = 	0.839420498
        if (JELEM	.eq.	973	)		coeff = 	0.797775161
        if (JELEM	.eq.	974	)		coeff = 	0.97196875
        if (JELEM	.eq.	975	)		coeff = 	0.917459308
        if (JELEM	.eq.	976	)		coeff = 	1.004567735
        if (JELEM	.eq.	977	)		coeff = 	0.851217152
        if (JELEM	.eq.	978	)		coeff = 	0.923398006
        if (JELEM	.eq.	979	)		coeff = 	0.915093327
        if (JELEM	.eq.	980	)		coeff = 	0.801664945
        if (JELEM	.eq.	981	)		coeff = 	0.969900178
        if (JELEM	.eq.	982	)		coeff = 	0.770849307
        if (JELEM	.eq.	983	)		coeff = 	0.957430239
        if (JELEM	.eq.	984	)		coeff = 	0.910770793
        if (JELEM	.eq.	985	)		coeff = 	0.933359627
        if (JELEM	.eq.	986	)		coeff = 	0.86933164
        if (JELEM	.eq.	987	)		coeff = 	0.938869204
        if (JELEM	.eq.	988	)		coeff = 	0.837823179
        if (JELEM	.eq.	989	)		coeff = 	0.973498697
        if (JELEM	.eq.	990	)		coeff = 	0.742336658
        if (JELEM	.eq.	991	)		coeff = 	1.118331307
        if (JELEM	.eq.	992	)		coeff = 	0.424939744
        if (JELEM	.eq.	993	)		coeff = 	0.866579638
        if (JELEM	.eq.	994	)		coeff = 	0.934711584
        if (JELEM	.eq.	995	)		coeff = 	1.086792906
        if (JELEM	.eq.	996	)		coeff = 	0.805469194
        if (JELEM	.eq.	997	)		coeff = 	0.867669093
        if (JELEM	.eq.	998	)		coeff = 	0.983272639
        if (JELEM	.eq.	999	)		coeff = 	1.033106396
        if (JELEM	.eq.	1000	)		coeff = 	1.018298163
        if (JELEM	.eq.	1001	)		coeff = 	0.926255074
        if (JELEM	.eq.	1002	)		coeff = 	0.883724745
        if (JELEM	.eq.	1003	)		coeff = 	0.786445795
        if (JELEM	.eq.	1004	)		coeff = 	1.05686855
        if (JELEM	.eq.	1005	)		coeff = 	0.952889185
        if (JELEM	.eq.	1006	)		coeff = 	0.980687447
        if (JELEM	.eq.	1007	)		coeff = 	1.045574465
        if (JELEM	.eq.	1008	)		coeff = 	0.8520149
        if (JELEM	.eq.	1009	)		coeff = 	0.833663013
        if (JELEM	.eq.	1010	)		coeff = 	0.993899981
        if (JELEM	.eq.	1011	)		coeff = 	1.01712879
        if (JELEM	.eq.	1012	)		coeff = 	0.961760909
        if (JELEM	.eq.	1013	)		coeff = 	1.03116883
        if (JELEM	.eq.	1014	)		coeff = 	1.079561841
        if (JELEM	.eq.	1015	)		coeff = 	0.990681304
        if (JELEM	.eq.	1016	)		coeff = 	1.028662205
        if (JELEM	.eq.	1017	)		coeff = 	1.044298876
        if (JELEM	.eq.	1018	)		coeff = 	0.903070633
        if (JELEM	.eq.	1019	)		coeff = 	0.676815873
        if (JELEM	.eq.	1020	)		coeff = 	1.067014513
        if (JELEM	.eq.	1021	)		coeff = 	0.744946196
        if (JELEM	.eq.	1022	)		coeff = 	0.930982997
        if (JELEM	.eq.	1023	)		coeff = 	1.06398284
        if (JELEM	.eq.	1024	)		coeff = 	0.877217762
        if (JELEM	.eq.	1025	)		coeff = 	0.899447042
        if (JELEM	.eq.	1026	)		coeff = 	1.012321671
        if (JELEM	.eq.	1027	)		coeff = 	0.947742769
        if (JELEM	.eq.	1028	)		coeff = 	1.177610159
        if (JELEM	.eq.	1029	)		coeff = 	1.181777012
        if (JELEM	.eq.	1030	)		coeff = 	0.829889677
        if (JELEM	.eq.	1031	)		coeff = 	1.105284036
        if (JELEM	.eq.	1032	)		coeff = 	0.883962067
        if (JELEM	.eq.	1033	)		coeff = 	0.842922252
        if (JELEM	.eq.	1034	)		coeff = 	0.833898586
        if (JELEM	.eq.	1035	)		coeff = 	0.874882011
        if (JELEM	.eq.	1036	)		coeff = 	0.695809413
        if (JELEM	.eq.	1037	)		coeff = 	1.02417666
        if (JELEM	.eq.	1038	)		coeff = 	0.504366518
        if (JELEM	.eq.	1039	)		coeff = 	0.936130909
        if (JELEM	.eq.	1040	)		coeff = 	0.89762408
        if (JELEM	.eq.	1041	)		coeff = 	0.697779693
        if (JELEM	.eq.	1042	)		coeff = 	1.05139822
        if (JELEM	.eq.	1043	)		coeff = 	0.607923461
        if (JELEM	.eq.	1044	)		coeff = 	1.073852696
        if (JELEM	.eq.	1045	)		coeff = 	0.704598419
        if (JELEM	.eq.	1046	)		coeff = 	0.676742996
        if (JELEM	.eq.	1047	)		coeff = 	1.071910675
        if (JELEM	.eq.	1048	)		coeff = 	1.14299281
        if (JELEM	.eq.	1049	)		coeff = 	1.034919588
        if (JELEM	.eq.	1050	)		coeff = 	0.939139389
        if (JELEM	.eq.	1051	)		coeff = 	0.842836459
        if (JELEM	.eq.	1052	)		coeff = 	0.918049706
        if (JELEM	.eq.	1053	)		coeff = 	1.017812618
        if (JELEM	.eq.	1054	)		coeff = 	0.977887755
        if (JELEM	.eq.	1055	)		coeff = 	0.801937489
        if (JELEM	.eq.	1056	)		coeff = 	0.90950084
        if (JELEM	.eq.	1057	)		coeff = 	1.033333307
        if (JELEM	.eq.	1058	)		coeff = 	0.950191819
        if (JELEM	.eq.	1059	)		coeff = 	0.765840024
        if (JELEM	.eq.	1060	)		coeff = 	1.072226737
        if (JELEM	.eq.	1061	)		coeff = 	0.816928124
        if (JELEM	.eq.	1062	)		coeff = 	1.069220553
        if (JELEM	.eq.	1063	)		coeff = 	1.111545313
        if (JELEM	.eq.	1064	)		coeff = 	1.167997662
        if (JELEM	.eq.	1065	)		coeff = 	1.209761657
        if (JELEM	.eq.	1066	)		coeff = 	0.653085345
        if (JELEM	.eq.	1067	)		coeff = 	0.766170111
        if (JELEM	.eq.	1068	)		coeff = 	1.147300404
        if (JELEM	.eq.	1069	)		coeff = 	0.891826552
        if (JELEM	.eq.	1070	)		coeff = 	0.626243645
        if (JELEM	.eq.	1071	)		coeff = 	1.053117799
        if (JELEM	.eq.	1072	)		coeff = 	1.120015714
        if (JELEM	.eq.	1073	)		coeff = 	1.003068001
        if (JELEM	.eq.	1074	)		coeff = 	0.863939172
        if (JELEM	.eq.	1075	)		coeff = 	0.918657541
        if (JELEM	.eq.	1076	)		coeff = 	0.971801821
        if (JELEM	.eq.	1077	)		coeff = 	0.741344175
        if (JELEM	.eq.	1078	)		coeff = 	1.244442558
        if (JELEM	.eq.	1079	)		coeff = 	0.951737245
        if (JELEM	.eq.	1080	)		coeff = 	1.085995441
        if (JELEM	.eq.	1081	)		coeff = 	0.95361552
        if (JELEM	.eq.	1082	)		coeff = 	0.931173567
        if (JELEM	.eq.	1083	)		coeff = 	0.903209379
        if (JELEM	.eq.	1084	)		coeff = 	1.058454874
        if (JELEM	.eq.	1085	)		coeff = 	1.009462997
        if (JELEM	.eq.	1086	)		coeff = 	1.013852658
        if (JELEM	.eq.	1087	)		coeff = 	0.913799192
        if (JELEM	.eq.	1088	)		coeff = 	0.922749475
        if (JELEM	.eq.	1089	)		coeff = 	0.884341674
        if (JELEM	.eq.	1090	)		coeff = 	1.075725789
        if (JELEM	.eq.	1091	)		coeff = 	0.899185871
        if (JELEM	.eq.	1092	)		coeff = 	0.846296318
        if (JELEM	.eq.	1093	)		coeff = 	0.990560346
        if (JELEM	.eq.	1094	)		coeff = 	0.990392335
        if (JELEM	.eq.	1095	)		coeff = 	0.882531729
        if (JELEM	.eq.	1096	)		coeff = 	1.046349542
        if (JELEM	.eq.	1097	)		coeff = 	0.691813393
        if (JELEM	.eq.	1098	)		coeff = 	1.090297873
        if (JELEM	.eq.	1099	)		coeff = 	0.951239601
        if (JELEM	.eq.	1100	)		coeff = 	1.1557584
        if (JELEM	.eq.	1101	)		coeff = 	1.019102828
        if (JELEM	.eq.	1102	)		coeff = 	1.015833918
        if (JELEM	.eq.	1103	)		coeff = 	0.6927823
        if (JELEM	.eq.	1104	)		coeff = 	1.058289596
        if (JELEM	.eq.	1105	)		coeff = 	0.958881422
        if (JELEM	.eq.	1106	)		coeff = 	0.974065509
        if (JELEM	.eq.	1107	)		coeff = 	0.816437672
        if (JELEM	.eq.	1108	)		coeff = 	1.190933284
        if (JELEM	.eq.	1109	)		coeff = 	0.835567569
        if (JELEM	.eq.	1110	)		coeff = 	0.940990009
        if (JELEM	.eq.	1111	)		coeff = 	1.151475404
        if (JELEM	.eq.	1112	)		coeff = 	0.950429375
        if (JELEM	.eq.	1113	)		coeff = 	0.973744958
        if (JELEM	.eq.	1114	)		coeff = 	0.958187983
        if (JELEM	.eq.	1115	)		coeff = 	0.896076876
        if (JELEM	.eq.	1116	)		coeff = 	0.953366673
        if (JELEM	.eq.	1117	)		coeff = 	0.859776043
        if (JELEM	.eq.	1118	)		coeff = 	0.700212864
        if (JELEM	.eq.	1119	)		coeff = 	0.829426639
        if (JELEM	.eq.	1120	)		coeff = 	0.976411977
        if (JELEM	.eq.	1121	)		coeff = 	1.1159059
        if (JELEM	.eq.	1122	)		coeff = 	0.715570329
        if (JELEM	.eq.	1123	)		coeff = 	0.996439196
        if (JELEM	.eq.	1124	)		coeff = 	0.851305174
        if (JELEM	.eq.	1125	)		coeff = 	0.826817784
        if (JELEM	.eq.	1126	)		coeff = 	0.826556849
        if (JELEM	.eq.	1127	)		coeff = 	0.921036099
        if (JELEM	.eq.	1128	)		coeff = 	1.014715623
        if (JELEM	.eq.	1129	)		coeff = 	0.856601132
        if (JELEM	.eq.	1130	)		coeff = 	0.963280291
        if (JELEM	.eq.	1131	)		coeff = 	0.949824791
        if (JELEM	.eq.	1132	)		coeff = 	0.962044461
        if (JELEM	.eq.	1133	)		coeff = 	0.981035626
        if (JELEM	.eq.	1134	)		coeff = 	1.071470785
        if (JELEM	.eq.	1135	)		coeff = 	1.03042724
        if (JELEM	.eq.	1136	)		coeff = 	0.775254545
        if (JELEM	.eq.	1137	)		coeff = 	0.778444139
        if (JELEM	.eq.	1138	)		coeff = 	0.958797742
        if (JELEM	.eq.	1139	)		coeff = 	0.917098155
        if (JELEM	.eq.	1140	)		coeff = 	0.871157482
        if (JELEM	.eq.	1141	)		coeff = 	1.08816918
        if (JELEM	.eq.	1142	)		coeff = 	0.925756018
        if (JELEM	.eq.	1143	)		coeff = 	1.063363132
        if (JELEM	.eq.	1144	)		coeff = 	1.135270515
        if (JELEM	.eq.	1145	)		coeff = 	1.134813974
        if (JELEM	.eq.	1146	)		coeff = 	0.706422563
        if (JELEM	.eq.	1147	)		coeff = 	0.948225599
        if (JELEM	.eq.	1148	)		coeff = 	1.083277497
        if (JELEM	.eq.	1149	)		coeff = 	0.901448068
        if (JELEM	.eq.	1150	)		coeff = 	0.925487854
        if (JELEM	.eq.	1151	)		coeff = 	1.019129401
        if (JELEM	.eq.	1152	)		coeff = 	0.89101994
        if (JELEM	.eq.	1153	)		coeff = 	0.990125033
        if (JELEM	.eq.	1154	)		coeff = 	1.053414423
        if (JELEM	.eq.	1155	)		coeff = 	1.045147046
        if (JELEM	.eq.	1156	)		coeff = 	0.993732302
        if (JELEM	.eq.	1157	)		coeff = 	1.165128206
        if (JELEM	.eq.	1158	)		coeff = 	1.007218676
        if (JELEM	.eq.	1159	)		coeff = 	0.936990716
        if (JELEM	.eq.	1160	)		coeff = 	0.934699602
        if (JELEM	.eq.	1161	)		coeff = 	1.132498382
        if (JELEM	.eq.	1162	)		coeff = 	0.926165533
        if (JELEM	.eq.	1163	)		coeff = 	0.963489595
        if (JELEM	.eq.	1164	)		coeff = 	0.89524011
        if (JELEM	.eq.	1165	)		coeff = 	0.836909937
        if (JELEM	.eq.	1166	)		coeff = 	1.065448921
        if (JELEM	.eq.	1167	)		coeff = 	0.667163788
        if (JELEM	.eq.	1168	)		coeff = 	0.943012538
        if (JELEM	.eq.	1169	)		coeff = 	0.739421515
        if (JELEM	.eq.	1170	)		coeff = 	0.901818818
        if (JELEM	.eq.	1171	)		coeff = 	1.081691848
        if (JELEM	.eq.	1172	)		coeff = 	1.144941388
        if (JELEM	.eq.	1173	)		coeff = 	1.196351084
        if (JELEM	.eq.	1174	)		coeff = 	1.062795918
        if (JELEM	.eq.	1175	)		coeff = 	0.945594506
        if (JELEM	.eq.	1176	)		coeff = 	1.08444522
        if (JELEM	.eq.	1177	)		coeff = 	1.127700144
        if (JELEM	.eq.	1178	)		coeff = 	0.935990688
        if (JELEM	.eq.	1179	)		coeff = 	0.608894973
        if (JELEM	.eq.	1180	)		coeff = 	1.077346274
        if (JELEM	.eq.	1181	)		coeff = 	1.080228996
        if (JELEM	.eq.	1182	)		coeff = 	1.076977923
        if (JELEM	.eq.	1183	)		coeff = 	1.123845142
        if (JELEM	.eq.	1184	)		coeff = 	0.870830894
        if (JELEM	.eq.	1185	)		coeff = 	0.97352821
        if (JELEM	.eq.	1186	)		coeff = 	0.95172355
        if (JELEM	.eq.	1187	)		coeff = 	0.917527864
        if (JELEM	.eq.	1188	)		coeff = 	0.770602409
        if (JELEM	.eq.	1189	)		coeff = 	1.097467645
        if (JELEM	.eq.	1190	)		coeff = 	1.026942135
        if (JELEM	.eq.	1191	)		coeff = 	1.136083374
        if (JELEM	.eq.	1192	)		coeff = 	1.086446876
        if (JELEM	.eq.	1193	)		coeff = 	1.131602627
        if (JELEM	.eq.	1194	)		coeff = 	1.067547418
        if (JELEM	.eq.	1195	)		coeff = 	0.846275553
        if (JELEM	.eq.	1196	)		coeff = 	0.757094746
        if (JELEM	.eq.	1197	)		coeff = 	0.948892004
        if (JELEM	.eq.	1198	)		coeff = 	0.800643204
        if (JELEM	.eq.	1199	)		coeff = 	1.049790325
        if (JELEM	.eq.	1200	)		coeff = 	0.979567054
        if (JELEM	.eq.	1201	)		coeff = 	0.99986012
        if (JELEM	.eq.	1202	)		coeff = 	0.880823862
        if (JELEM	.eq.	1203	)		coeff = 	0.69935729
        if (JELEM	.eq.	1204	)		coeff = 	1.054060066
        if (JELEM	.eq.	1205	)		coeff = 	0.854288477
        if (JELEM	.eq.	1206	)		coeff = 	0.948346172
        if (JELEM	.eq.	1207	)		coeff = 	0.850718553
        if (JELEM	.eq.	1208	)		coeff = 	0.772519599
        if (JELEM	.eq.	1209	)		coeff = 	1.046424171
        if (JELEM	.eq.	1210	)		coeff = 	0.706683589
        if (JELEM	.eq.	1211	)		coeff = 	1.173904569
        if (JELEM	.eq.	1212	)		coeff = 	0.933594081
        if (JELEM	.eq.	1213	)		coeff = 	0.763957682
        if (JELEM	.eq.	1214	)		coeff = 	1.075476682
        if (JELEM	.eq.	1215	)		coeff = 	0.980622375
        if (JELEM	.eq.	1216	)		coeff = 	1.01405312
        if (JELEM	.eq.	1217	)		coeff = 	0.931399545
        if (JELEM	.eq.	1218	)		coeff = 	0.812520856
        if (JELEM	.eq.	1219	)		coeff = 	0.991939507
        if (JELEM	.eq.	1220	)		coeff = 	0.911349249
        if (JELEM	.eq.	1221	)		coeff = 	0.658767622
        if (JELEM	.eq.	1222	)		coeff = 	0.811558668
        if (JELEM	.eq.	1223	)		coeff = 	0.821055968
        if (JELEM	.eq.	1224	)		coeff = 	0.68042414
        if (JELEM	.eq.	1225	)		coeff = 	1.068390217
        if (JELEM	.eq.	1226	)		coeff = 	0.891734513
        if (JELEM	.eq.	1227	)		coeff = 	0.862790394
        if (JELEM	.eq.	1228	)		coeff = 	0.817252281
        if (JELEM	.eq.	1229	)		coeff = 	0.887288679
        if (JELEM	.eq.	1230	)		coeff = 	1.01295745
        if (JELEM	.eq.	1231	)		coeff = 	1.163623391
        if (JELEM	.eq.	1232	)		coeff = 	1.098394697
        if (JELEM	.eq.	1233	)		coeff = 	1.029288892
        if (JELEM	.eq.	1234	)		coeff = 	0.805062205
        if (JELEM	.eq.	1235	)		coeff = 	0.753451587
        if (JELEM	.eq.	1236	)		coeff = 	0.972283379
        if (JELEM	.eq.	1237	)		coeff = 	0.572437467
        if (JELEM	.eq.	1238	)		coeff = 	1.16476131
        if (JELEM	.eq.	1239	)		coeff = 	0.789837976
        if (JELEM	.eq.	1240	)		coeff = 	0.945298477
        if (JELEM	.eq.	1241	)		coeff = 	0.630101377
        if (JELEM	.eq.	1242	)		coeff = 	0.775634551
        if (JELEM	.eq.	1243	)		coeff = 	1.078043273
        if (JELEM	.eq.	1244	)		coeff = 	1.004941763
        if (JELEM	.eq.	1245	)		coeff = 	1.139730439
        if (JELEM	.eq.	1246	)		coeff = 	0.974254448
        if (JELEM	.eq.	1247	)		coeff = 	0.97527528
        if (JELEM	.eq.	1248	)		coeff = 	1.01412095
        if (JELEM	.eq.	1249	)		coeff = 	0.739556252
        if (JELEM	.eq.	1250	)		coeff = 	0.95709236
        if (JELEM	.eq.	1251	)		coeff = 	0.874536659
        if (JELEM	.eq.	1252	)		coeff = 	0.89412869
        if (JELEM	.eq.	1253	)		coeff = 	0.809388365
        if (JELEM	.eq.	1254	)		coeff = 	1.038672747
        if (JELEM	.eq.	1255	)		coeff = 	0.845980232
        if (JELEM	.eq.	1256	)		coeff = 	0.985568665
        if (JELEM	.eq.	1257	)		coeff = 	1.141786818
        if (JELEM	.eq.	1258	)		coeff = 	0.781065762
        if (JELEM	.eq.	1259	)		coeff = 	0.639838026
        if (JELEM	.eq.	1260	)		coeff = 	1.100140799
        if (JELEM	.eq.	1261	)		coeff = 	0.952349025
        if (JELEM	.eq.	1262	)		coeff = 	0.853002999
        if (JELEM	.eq.	1263	)		coeff = 	1.104040582
        if (JELEM	.eq.	1264	)		coeff = 	1.164666311
        if (JELEM	.eq.	1265	)		coeff = 	0.9752739
        if (JELEM	.eq.	1266	)		coeff = 	1.0933368
        if (JELEM	.eq.	1267	)		coeff = 	0.64185982
        if (JELEM	.eq.	1268	)		coeff = 	0.705986723
        if (JELEM	.eq.	1269	)		coeff = 	0.716193538
        if (JELEM	.eq.	1270	)		coeff = 	0.598098761
        if (JELEM	.eq.	1271	)		coeff = 	0.882134805
        if (JELEM	.eq.	1272	)		coeff = 	0.877999156
        if (JELEM	.eq.	1273	)		coeff = 	0.998834268
        if (JELEM	.eq.	1274	)		coeff = 	1.134649288
        if (JELEM	.eq.	1275	)		coeff = 	0.974478853
        if (JELEM	.eq.	1276	)		coeff = 	0.940211928
        if (JELEM	.eq.	1277	)		coeff = 	0.941212195
        if (JELEM	.eq.	1278	)		coeff = 	0.836806806
        if (JELEM	.eq.	1279	)		coeff = 	0.918933678
        if (JELEM	.eq.	1280	)		coeff = 	0.695010926
        if (JELEM	.eq.	1281	)		coeff = 	0.907521484
        if (JELEM	.eq.	1282	)		coeff = 	0.774263923
        if (JELEM	.eq.	1283	)		coeff = 	1.097056282
        if (JELEM	.eq.	1284	)		coeff = 	0.82096677
        if (JELEM	.eq.	1285	)		coeff = 	0.622850901
        if (JELEM	.eq.	1286	)		coeff = 	1.048825448
        if (JELEM	.eq.	1287	)		coeff = 	0.741419092
        if (JELEM	.eq.	1288	)		coeff = 	0.889251566
        if (JELEM	.eq.	1289	)		coeff = 	0.867753501
        if (JELEM	.eq.	1290	)		coeff = 	0.885256697
        if (JELEM	.eq.	1291	)		coeff = 	0.830928752
        if (JELEM	.eq.	1292	)		coeff = 	0.944923123
        if (JELEM	.eq.	1293	)		coeff = 	0.853469127
        if (JELEM	.eq.	1294	)		coeff = 	0.953993849
        if (JELEM	.eq.	1295	)		coeff = 	0.747590241
        if (JELEM	.eq.	1296	)		coeff = 	1.096565145
        if (JELEM	.eq.	1297	)		coeff = 	0.795211346
        if (JELEM	.eq.	1298	)		coeff = 	1.071897056
        if (JELEM	.eq.	1299	)		coeff = 	0.848619552
        if (JELEM	.eq.	1300	)		coeff = 	0.969105681
        if (JELEM	.eq.	1301	)		coeff = 	0.917180496
        if (JELEM	.eq.	1302	)		coeff = 	1.001305621
        if (JELEM	.eq.	1303	)		coeff = 	0.988475708
        if (JELEM	.eq.	1304	)		coeff = 	0.865408794
        if (JELEM	.eq.	1305	)		coeff = 	0.680358686
        if (JELEM	.eq.	1306	)		coeff = 	0.774296896
        if (JELEM	.eq.	1307	)		coeff = 	0.915791693
        if (JELEM	.eq.	1308	)		coeff = 	0.65226881
        if (JELEM	.eq.	1309	)		coeff = 	0.906378658
        if (JELEM	.eq.	1310	)		coeff = 	0.841723757
        if (JELEM	.eq.	1311	)		coeff = 	0.756284674
        if (JELEM	.eq.	1312	)		coeff = 	1.050644034
        if (JELEM	.eq.	1313	)		coeff = 	1.109122192
        if (JELEM	.eq.	1314	)		coeff = 	0.811696231
        if (JELEM	.eq.	1315	)		coeff = 	0.980882877
        if (JELEM	.eq.	1316	)		coeff = 	0.86680652
        if (JELEM	.eq.	1317	)		coeff = 	0.711046714
        if (JELEM	.eq.	1318	)		coeff = 	1.024704819
        if (JELEM	.eq.	1319	)		coeff = 	1.019176235
        if (JELEM	.eq.	1320	)		coeff = 	0.937827522
        if (JELEM	.eq.	1321	)		coeff = 	0.936797975
        if (JELEM	.eq.	1322	)		coeff = 	0.903910954
        if (JELEM	.eq.	1323	)		coeff = 	0.920716698
        if (JELEM	.eq.	1324	)		coeff = 	1.078116751
        if (JELEM	.eq.	1325	)		coeff = 	0.983094434
        if (JELEM	.eq.	1326	)		coeff = 	0.999600504
        if (JELEM	.eq.	1327	)		coeff = 	1.122605873
        if (JELEM	.eq.	1328	)		coeff = 	1.046872967
        if (JELEM	.eq.	1329	)		coeff = 	0.523911699
        if (JELEM	.eq.	1330	)		coeff = 	0.745175967
        if (JELEM	.eq.	1331	)		coeff = 	0.834296328
        if (JELEM	.eq.	1332	)		coeff = 	1.090602967
        if (JELEM	.eq.	1333	)		coeff = 	0.937598355
        if (JELEM	.eq.	1334	)		coeff = 	1.143419883
        if (JELEM	.eq.	1335	)		coeff = 	0.88917349
        if (JELEM	.eq.	1336	)		coeff = 	0.815532978
        if (JELEM	.eq.	1337	)		coeff = 	0.742893545
        if (JELEM	.eq.	1338	)		coeff = 	1.08466948
        if (JELEM	.eq.	1339	)		coeff = 	1.021756235
        if (JELEM	.eq.	1340	)		coeff = 	0.840317219
        if (JELEM	.eq.	1341	)		coeff = 	0.959605635
        if (JELEM	.eq.	1342	)		coeff = 	0.938013579
        if (JELEM	.eq.	1343	)		coeff = 	0.774792253
        if (JELEM	.eq.	1344	)		coeff = 	0.573964068
        if (JELEM	.eq.	1345	)		coeff = 	1.187538996
        if (JELEM	.eq.	1346	)		coeff = 	1.005643156
        if (JELEM	.eq.	1347	)		coeff = 	0.943158543
        if (JELEM	.eq.	1348	)		coeff = 	0.856862618
        if (JELEM	.eq.	1349	)		coeff = 	0.646309609
        if (JELEM	.eq.	1350	)		coeff = 	0.881218473
        if (JELEM	.eq.	1351	)		coeff = 	1.033070772
        if (JELEM	.eq.	1352	)		coeff = 	0.687821203
        if (JELEM	.eq.	1353	)		coeff = 	0.746922752
        if (JELEM	.eq.	1354	)		coeff = 	0.969752256
        if (JELEM	.eq.	1355	)		coeff = 	0.863944937
        if (JELEM	.eq.	1356	)		coeff = 	0.934771408
        if (JELEM	.eq.	1357	)		coeff = 	0.528910714
        if (JELEM	.eq.	1358	)		coeff = 	0.992743892
        if (JELEM	.eq.	1359	)		coeff = 	0.771814471
        if (JELEM	.eq.	1360	)		coeff = 	0.761326268
        if (JELEM	.eq.	1361	)		coeff = 	0.993859427
        if (JELEM	.eq.	1362	)		coeff = 	0.866593832
        if (JELEM	.eq.	1363	)		coeff = 	0.741515376
        if (JELEM	.eq.	1364	)		coeff = 	0.750627863
        if (JELEM	.eq.	1365	)		coeff = 	1.121786611
        if (JELEM	.eq.	1366	)		coeff = 	1.165332991
        if (JELEM	.eq.	1367	)		coeff = 	0.929650305
        if (JELEM	.eq.	1368	)		coeff = 	1.212737293
        if (JELEM	.eq.	1369	)		coeff = 	1.061980756
        if (JELEM	.eq.	1370	)		coeff = 	1.184633782
        if (JELEM	.eq.	1371	)		coeff = 	0.833460204
        if (JELEM	.eq.	1372	)		coeff = 	1.065110941
        if (JELEM	.eq.	1373	)		coeff = 	0.686102494
        if (JELEM	.eq.	1374	)		coeff = 	0.845463408
        if (JELEM	.eq.	1375	)		coeff = 	0.819106165
        if (JELEM	.eq.	1376	)		coeff = 	0.980350216
        if (JELEM	.eq.	1377	)		coeff = 	1.023397373
        if (JELEM	.eq.	1378	)		coeff = 	0.887596895
        if (JELEM	.eq.	1379	)		coeff = 	0.790545098
        if (JELEM	.eq.	1380	)		coeff = 	0.780136476
        if (JELEM	.eq.	1381	)		coeff = 	0.757693571
        if (JELEM	.eq.	1382	)		coeff = 	0.801374747
        if (JELEM	.eq.	1383	)		coeff = 	0.960316532
        if (JELEM	.eq.	1384	)		coeff = 	0.95915625
        if (JELEM	.eq.	1385	)		coeff = 	1.133666696
        if (JELEM	.eq.	1386	)		coeff = 	0.969902505
        if (JELEM	.eq.	1387	)		coeff = 	0.921124998
        if (JELEM	.eq.	1388	)		coeff = 	1.103985567
        if (JELEM	.eq.	1389	)		coeff = 	1.008226802
        if (JELEM	.eq.	1390	)		coeff = 	1.036442157
        if (JELEM	.eq.	1391	)		coeff = 	0.926811447
        if (JELEM	.eq.	1392	)		coeff = 	0.977200036
        if (JELEM	.eq.	1393	)		coeff = 	0.984040888
        if (JELEM	.eq.	1394	)		coeff = 	1.143080291
        if (JELEM	.eq.	1395	)		coeff = 	1.038056766
        if (JELEM	.eq.	1396	)		coeff = 	1.043871181
        if (JELEM	.eq.	1397	)		coeff = 	0.875049881
        if (JELEM	.eq.	1398	)		coeff = 	0.896500673
        if (JELEM	.eq.	1399	)		coeff = 	1.098939338
        if (JELEM	.eq.	1400	)		coeff = 	1.09362825
        if (JELEM	.eq.	1401	)		coeff = 	1.124212366
        if (JELEM	.eq.	1402	)		coeff = 	0.862960486
        if (JELEM	.eq.	1403	)		coeff = 	0.935521142
        if (JELEM	.eq.	1404	)		coeff = 	0.807668012
        if (JELEM	.eq.	1405	)		coeff = 	1.106259752
        if (JELEM	.eq.	1406	)		coeff = 	1.03922074
        if (JELEM	.eq.	1407	)		coeff = 	1.143160795
        if (JELEM	.eq.	1408	)		coeff = 	0.897435155
        if (JELEM	.eq.	1409	)		coeff = 	1.031556218
        if (JELEM	.eq.	1410	)		coeff = 	0.875034078
        if (JELEM	.eq.	1411	)		coeff = 	1.093073071
        if (JELEM	.eq.	1412	)		coeff = 	0.946894619
        if (JELEM	.eq.	1413	)		coeff = 	1.151053128
        if (JELEM	.eq.	1414	)		coeff = 	1.109098847
        if (JELEM	.eq.	1415	)		coeff = 	0.965845907
        if (JELEM	.eq.	1416	)		coeff = 	0.874352609
        if (JELEM	.eq.	1417	)		coeff = 	0.918957338
        if (JELEM	.eq.	1418	)		coeff = 	1.05927521
        if (JELEM	.eq.	1419	)		coeff = 	0.912335516
        if (JELEM	.eq.	1420	)		coeff = 	0.879465213
        if (JELEM	.eq.	1421	)		coeff = 	1.018255801
        if (JELEM	.eq.	1422	)		coeff = 	0.826230791
        if (JELEM	.eq.	1423	)		coeff = 	1.070595088
        if (JELEM	.eq.	1424	)		coeff = 	1.044962334
        if (JELEM	.eq.	1425	)		coeff = 	0.921432736
        if (JELEM	.eq.	1426	)		coeff = 	0.783007553
        if (JELEM	.eq.	1427	)		coeff = 	0.81780437
        if (JELEM	.eq.	1428	)		coeff = 	0.92315961
        if (JELEM	.eq.	1429	)		coeff = 	1.253110175
        if (JELEM	.eq.	1430	)		coeff = 	0.699112169
        if (JELEM	.eq.	1431	)		coeff = 	0.833344029
        if (JELEM	.eq.	1432	)		coeff = 	0.989358352
        if (JELEM	.eq.	1433	)		coeff = 	0.778778973
        if (JELEM	.eq.	1434	)		coeff = 	0.890078664
        if (JELEM	.eq.	1435	)		coeff = 	1.044559124
        if (JELEM	.eq.	1436	)		coeff = 	1.074824333
        if (JELEM	.eq.	1437	)		coeff = 	0.788562098
        if (JELEM	.eq.	1438	)		coeff = 	1.014141598
        if (JELEM	.eq.	1439	)		coeff = 	1.014658432
        if (JELEM	.eq.	1440	)		coeff = 	0.989073292
        if (JELEM	.eq.	1441	)		coeff = 	0.848959773
        if (JELEM	.eq.	1442	)		coeff = 	1.076797542
        if (JELEM	.eq.	1443	)		coeff = 	0.884184292
        if (JELEM	.eq.	1444	)		coeff = 	0.902753915
        if (JELEM	.eq.	1445	)		coeff = 	0.816119198
        if (JELEM	.eq.	1446	)		coeff = 	0.996381403
        if (JELEM	.eq.	1447	)		coeff = 	1.102958261
        if (JELEM	.eq.	1448	)		coeff = 	0.789813056
        if (JELEM	.eq.	1449	)		coeff = 	0.873207246
        if (JELEM	.eq.	1450	)		coeff = 	0.928702653
        if (JELEM	.eq.	1451	)		coeff = 	0.860125915
        if (JELEM	.eq.	1452	)		coeff = 	0.85395477
        if (JELEM	.eq.	1453	)		coeff = 	0.935901522
        if (JELEM	.eq.	1454	)		coeff = 	0.828504878
        if (JELEM	.eq.	1455	)		coeff = 	0.891768746
        if (JELEM	.eq.	1456	)		coeff = 	0.942996285
        if (JELEM	.eq.	1457	)		coeff = 	0.885198278
        if (JELEM	.eq.	1458	)		coeff = 	0.680250923
        if (JELEM	.eq.	1459	)		coeff = 	0.908872043
        if (JELEM	.eq.	1460	)		coeff = 	0.817514745
        if (JELEM	.eq.	1461	)		coeff = 	1.065897935
        if (JELEM	.eq.	1462	)		coeff = 	0.855050894
        if (JELEM	.eq.	1463	)		coeff = 	1.052599239
        if (JELEM	.eq.	1464	)		coeff = 	0.983586029
        if (JELEM	.eq.	1465	)		coeff = 	0.834247596
        if (JELEM	.eq.	1466	)		coeff = 	0.989144659
        if (JELEM	.eq.	1467	)		coeff = 	1.15185253
        if (JELEM	.eq.	1468	)		coeff = 	1.060860194
        if (JELEM	.eq.	1469	)		coeff = 	1.154209149
        if (JELEM	.eq.	1470	)		coeff = 	0.904400323
        if (JELEM	.eq.	1471	)		coeff = 	1.041014388
        if (JELEM	.eq.	1472	)		coeff = 	0.982428633
        if (JELEM	.eq.	1473	)		coeff = 	1.073313264
        if (JELEM	.eq.	1474	)		coeff = 	0.769382194
        if (JELEM	.eq.	1475	)		coeff = 	1.00917859
        if (JELEM	.eq.	1476	)		coeff = 	1.059191008
        if (JELEM	.eq.	1477	)		coeff = 	0.831691085
        if (JELEM	.eq.	1478	)		coeff = 	1.17808397
        if (JELEM	.eq.	1479	)		coeff = 	1.157593703
        if (JELEM	.eq.	1480	)		coeff = 	0.94485587
        if (JELEM	.eq.	1481	)		coeff = 	0.71585013
        if (JELEM	.eq.	1482	)		coeff = 	0.918631078
        if (JELEM	.eq.	1483	)		coeff = 	0.919021062
        if (JELEM	.eq.	1484	)		coeff = 	0.902535275
        if (JELEM	.eq.	1485	)		coeff = 	0.730720916
        if (JELEM	.eq.	1486	)		coeff = 	1.05351672
        if (JELEM	.eq.	1487	)		coeff = 	0.768984366
        if (JELEM	.eq.	1488	)		coeff = 	1.148135698
        if (JELEM	.eq.	1489	)		coeff = 	0.821886136
        if (JELEM	.eq.	1490	)		coeff = 	1.166304309
        if (JELEM	.eq.	1491	)		coeff = 	0.98406034
        if (JELEM	.eq.	1492	)		coeff = 	0.921456368
        if (JELEM	.eq.	1493	)		coeff = 	1.062616169
        if (JELEM	.eq.	1494	)		coeff = 	0.831288435
        if (JELEM	.eq.	1495	)		coeff = 	1.056816436
        if (JELEM	.eq.	1496	)		coeff = 	1.072174991
        if (JELEM	.eq.	1497	)		coeff = 	0.914227937
        if (JELEM	.eq.	1498	)		coeff = 	0.70460417
        if (JELEM	.eq.	1499	)		coeff = 	0.895391755
        if (JELEM	.eq.	1500	)		coeff = 	1.055249292
        if (JELEM	.eq.	1501	)		coeff = 	0.740466134
        if (JELEM	.eq.	1502	)		coeff = 	0.995822734
        if (JELEM	.eq.	1503	)		coeff = 	1.163262675
        if (JELEM	.eq.	1504	)		coeff = 	1.050339917
        if (JELEM	.eq.	1505	)		coeff = 	0.747752534
        if (JELEM	.eq.	1506	)		coeff = 	0.782403702
        if (JELEM	.eq.	1507	)		coeff = 	0.92597112
        if (JELEM	.eq.	1508	)		coeff = 	0.93385615
        if (JELEM	.eq.	1509	)		coeff = 	1.271292838
        if (JELEM	.eq.	1510	)		coeff = 	1.094853279
        if (JELEM	.eq.	1511	)		coeff = 	1.055352076
        if (JELEM	.eq.	1512	)		coeff = 	0.812629122
        if (JELEM	.eq.	1513	)		coeff = 	0.99331367
        if (JELEM	.eq.	1514	)		coeff = 	1.017121428
        if (JELEM	.eq.	1515	)		coeff = 	1.03728593
        if (JELEM	.eq.	1516	)		coeff = 	1.024707819
        if (JELEM	.eq.	1517	)		coeff = 	0.979632803
        if (JELEM	.eq.	1518	)		coeff = 	1.064448017
        if (JELEM	.eq.	1519	)		coeff = 	0.931440331
        if (JELEM	.eq.	1520	)		coeff = 	0.896983249
        if (JELEM	.eq.	1521	)		coeff = 	0.885346574
        if (JELEM	.eq.	1522	)		coeff = 	1.004817
        if (JELEM	.eq.	1523	)		coeff = 	0.677082992
        if (JELEM	.eq.	1524	)		coeff = 	1.044037809
        if (JELEM	.eq.	1525	)		coeff = 	0.856638742
        if (JELEM	.eq.	1526	)		coeff = 	0.992430636
        if (JELEM	.eq.	1527	)		coeff = 	0.800555678
        if (JELEM	.eq.	1528	)		coeff = 	0.998606282
        if (JELEM	.eq.	1529	)		coeff = 	0.875716989
        if (JELEM	.eq.	1530	)		coeff = 	1.007584726
        if (JELEM	.eq.	1531	)		coeff = 	0.775688205
        if (JELEM	.eq.	1532	)		coeff = 	0.855521951
        if (JELEM	.eq.	1533	)		coeff = 	0.977525862
        if (JELEM	.eq.	1534	)		coeff = 	0.988908921
        if (JELEM	.eq.	1535	)		coeff = 	0.99828265
        if (JELEM	.eq.	1536	)		coeff = 	0.974546905
        if (JELEM	.eq.	1537	)		coeff = 	0.995785272
        if (JELEM	.eq.	1538	)		coeff = 	1.128011184
        if (JELEM	.eq.	1539	)		coeff = 	0.975549725
        if (JELEM	.eq.	1540	)		coeff = 	0.917005512
        if (JELEM	.eq.	1541	)		coeff = 	1.021783708
        if (JELEM	.eq.	1542	)		coeff = 	0.986891301
        if (JELEM	.eq.	1543	)		coeff = 	1.109274658
        if (JELEM	.eq.	1544	)		coeff = 	0.96314108
        if (JELEM	.eq.	1545	)		coeff = 	1.026172466
        if (JELEM	.eq.	1546	)		coeff = 	1.193524462
        if (JELEM	.eq.	1547	)		coeff = 	0.769417748
        if (JELEM	.eq.	1548	)		coeff = 	1.05708928
        if (JELEM	.eq.	1549	)		coeff = 	0.886172549
        if (JELEM	.eq.	1550	)		coeff = 	1.081259652
        if (JELEM	.eq.	1551	)		coeff = 	0.633742503
        if (JELEM	.eq.	1552	)		coeff = 	0.954301619
        if (JELEM	.eq.	1553	)		coeff = 	0.763071014
        if (JELEM	.eq.	1554	)		coeff = 	0.86708924
        if (JELEM	.eq.	1555	)		coeff = 	0.956190472
        if (JELEM	.eq.	1556	)		coeff = 	0.916745558
        if (JELEM	.eq.	1557	)		coeff = 	1.006339557
        if (JELEM	.eq.	1558	)		coeff = 	0.946119053
        if (JELEM	.eq.	1559	)		coeff = 	0.830816116
        if (JELEM	.eq.	1560	)		coeff = 	1.032922391
        if (JELEM	.eq.	1561	)		coeff = 	1.026074993
        if (JELEM	.eq.	1562	)		coeff = 	1.00672893
        if (JELEM	.eq.	1563	)		coeff = 	0.77668383
        if (JELEM	.eq.	1564	)		coeff = 	1.113994777
        if (JELEM	.eq.	1565	)		coeff = 	1.058557752
        if (JELEM	.eq.	1566	)		coeff = 	1.022498578
        if (JELEM	.eq.	1567	)		coeff = 	1.00557522
        if (JELEM	.eq.	1568	)		coeff = 	0.866969479
        if (JELEM	.eq.	1569	)		coeff = 	0.767715812
        if (JELEM	.eq.	1570	)		coeff = 	1.145179389
        if (JELEM	.eq.	1571	)		coeff = 	1.123796226
        if (JELEM	.eq.	1572	)		coeff = 	1.197708546
        if (JELEM	.eq.	1573	)		coeff = 	1.064056648
        if (JELEM	.eq.	1574	)		coeff = 	1.175213307
        if (JELEM	.eq.	1575	)		coeff = 	0.965903442
        if (JELEM	.eq.	1576	)		coeff = 	1.044172858
        if (JELEM	.eq.	1577	)		coeff = 	0.733149575
        if (JELEM	.eq.	1578	)		coeff = 	0.87200357
        if (JELEM	.eq.	1579	)		coeff = 	0.959120199
        if (JELEM	.eq.	1580	)		coeff = 	1.124137777
        if (JELEM	.eq.	1581	)		coeff = 	0.944769016
        if (JELEM	.eq.	1582	)		coeff = 	0.893547974
        if (JELEM	.eq.	1583	)		coeff = 	0.91655697
        if (JELEM	.eq.	1584	)		coeff = 	0.982540685
        if (JELEM	.eq.	1585	)		coeff = 	1.103580613
        if (JELEM	.eq.	1586	)		coeff = 	0.81908085
        if (JELEM	.eq.	1587	)		coeff = 	1.041935665
        if (JELEM	.eq.	1588	)		coeff = 	0.952049436
        if (JELEM	.eq.	1589	)		coeff = 	0.826617932
        if (JELEM	.eq.	1590	)		coeff = 	0.951823748
        if (JELEM	.eq.	1591	)		coeff = 	1.022519458
        if (JELEM	.eq.	1592	)		coeff = 	0.787946846
        if (JELEM	.eq.	1593	)		coeff = 	1.140260732
        if (JELEM	.eq.	1594	)		coeff = 	0.865087242
        if (JELEM	.eq.	1595	)		coeff = 	0.753926504
        if (JELEM	.eq.	1596	)		coeff = 	0.967507646
        if (JELEM	.eq.	1597	)		coeff = 	0.96570009
        if (JELEM	.eq.	1598	)		coeff = 	1.081529901
        if (JELEM	.eq.	1599	)		coeff = 	0.666529004
        if (JELEM	.eq.	1600	)		coeff = 	0.776447753
        if (JELEM	.eq.	1601	)		coeff = 	0.959125838
        if (JELEM	.eq.	1602	)		coeff = 	0.986745396
        if (JELEM	.eq.	1603	)		coeff = 	1.095002707
        if (JELEM	.eq.	1604	)		coeff = 	0.726380258
        if (JELEM	.eq.	1605	)		coeff = 	1.228510349
        if (JELEM	.eq.	1606	)		coeff = 	1.066991988
        if (JELEM	.eq.	1607	)		coeff = 	1.015037256
        if (JELEM	.eq.	1608	)		coeff = 	1.146820182
        if (JELEM	.eq.	1609	)		coeff = 	1.085958563
        if (JELEM	.eq.	1610	)		coeff = 	0.865740645
        if (JELEM	.eq.	1611	)		coeff = 	0.961298265
        if (JELEM	.eq.	1612	)		coeff = 	1.010193541
        if (JELEM	.eq.	1613	)		coeff = 	1.046696775
        if (JELEM	.eq.	1614	)		coeff = 	0.971973545
        if (JELEM	.eq.	1615	)		coeff = 	1.067391767
        if (JELEM	.eq.	1616	)		coeff = 	1.014969929
        if (JELEM	.eq.	1617	)		coeff = 	1.03647667
        if (JELEM	.eq.	1618	)		coeff = 	0.810509538
        if (JELEM	.eq.	1619	)		coeff = 	0.880687515
        if (JELEM	.eq.	1620	)		coeff = 	1.011754878
        if (JELEM	.eq.	1621	)		coeff = 	0.926792318
        if (JELEM	.eq.	1622	)		coeff = 	1.027807358
        if (JELEM	.eq.	1623	)		coeff = 	1.036476973
        if (JELEM	.eq.	1624	)		coeff = 	1.037980565
        if (JELEM	.eq.	1625	)		coeff = 	0.888952005
        if (JELEM	.eq.	1626	)		coeff = 	0.948350337
        if (JELEM	.eq.	1627	)		coeff = 	1.125074999
        if (JELEM	.eq.	1628	)		coeff = 	1.141619661
        if (JELEM	.eq.	1629	)		coeff = 	1.038415476
        if (JELEM	.eq.	1630	)		coeff = 	0.975669022
        if (JELEM	.eq.	1631	)		coeff = 	1.029081747
        if (JELEM	.eq.	1632	)		coeff = 	0.888240592
        if (JELEM	.eq.	1633	)		coeff = 	0.690596032
        if (JELEM	.eq.	1634	)		coeff = 	0.84349239
        if (JELEM	.eq.	1635	)		coeff = 	0.905592683
        if (JELEM	.eq.	1636	)		coeff = 	0.853999468
        if (JELEM	.eq.	1637	)		coeff = 	0.857355402
        if (JELEM	.eq.	1638	)		coeff = 	0.924658808
        if (JELEM	.eq.	1639	)		coeff = 	0.843969545
        if (JELEM	.eq.	1640	)		coeff = 	0.991667058
        if (JELEM	.eq.	1641	)		coeff = 	0.917472573
        if (JELEM	.eq.	1642	)		coeff = 	1.042640582
        if (JELEM	.eq.	1643	)		coeff = 	1.026961106
        if (JELEM	.eq.	1644	)		coeff = 	1.18718697
        if (JELEM	.eq.	1645	)		coeff = 	1.00701751
        if (JELEM	.eq.	1646	)		coeff = 	1.08729652
        if (JELEM	.eq.	1647	)		coeff = 	0.985242593
        if (JELEM	.eq.	1648	)		coeff = 	1.081880505
        if (JELEM	.eq.	1649	)		coeff = 	0.811207981
        if (JELEM	.eq.	1650	)		coeff = 	0.861062238
        if (JELEM	.eq.	1651	)		coeff = 	1.110565455
        if (JELEM	.eq.	1652	)		coeff = 	0.816573188
        if (JELEM	.eq.	1653	)		coeff = 	1.050430965
        if (JELEM	.eq.	1654	)		coeff = 	1.105786341
        if (JELEM	.eq.	1655	)		coeff = 	0.894895346
        if (JELEM	.eq.	1656	)		coeff = 	0.681424109
        if (JELEM	.eq.	1657	)		coeff = 	0.820344293
        if (JELEM	.eq.	1658	)		coeff = 	0.9106155
        if (JELEM	.eq.	1659	)		coeff = 	1.013577308
        if (JELEM	.eq.	1660	)		coeff = 	1.031581078
        if (JELEM	.eq.	1661	)		coeff = 	0.977542216
        if (JELEM	.eq.	1662	)		coeff = 	1.005994951
        if (JELEM	.eq.	1663	)		coeff = 	0.775180106
        if (JELEM	.eq.	1664	)		coeff = 	1.225591759
        if (JELEM	.eq.	1665	)		coeff = 	0.880648551
        if (JELEM	.eq.	1666	)		coeff = 	1.010194106
        if (JELEM	.eq.	1667	)		coeff = 	1.021713348
        if (JELEM	.eq.	1668	)		coeff = 	0.900844494
        if (JELEM	.eq.	1669	)		coeff = 	0.731672915
        if (JELEM	.eq.	1670	)		coeff = 	0.761510934
        if (JELEM	.eq.	1671	)		coeff = 	0.480143143
        if (JELEM	.eq.	1672	)		coeff = 	1.127884447
        if (JELEM	.eq.	1673	)		coeff = 	1.093421893
        if (JELEM	.eq.	1674	)		coeff = 	0.607715154
        if (JELEM	.eq.	1675	)		coeff = 	1.115951162
        if (JELEM	.eq.	1676	)		coeff = 	0.88509875
        if (JELEM	.eq.	1677	)		coeff = 	0.721808456
        if (JELEM	.eq.	1678	)		coeff = 	1.087286403
        if (JELEM	.eq.	1679	)		coeff = 	0.769417802
        if (JELEM	.eq.	1680	)		coeff = 	1.186358849
        if (JELEM	.eq.	1681	)		coeff = 	1.008604017
        if (JELEM	.eq.	1682	)		coeff = 	1.046107197
        if (JELEM	.eq.	1683	)		coeff = 	0.594276488
        if (JELEM	.eq.	1684	)		coeff = 	0.799256863
        if (JELEM	.eq.	1685	)		coeff = 	0.832190045
        if (JELEM	.eq.	1686	)		coeff = 	0.754527954
        if (JELEM	.eq.	1687	)		coeff = 	0.840190985
        if (JELEM	.eq.	1688	)		coeff = 	0.806231212
        if (JELEM	.eq.	1689	)		coeff = 	1.016069516
        if (JELEM	.eq.	1690	)		coeff = 	0.859345752
        if (JELEM	.eq.	1691	)		coeff = 	0.894189224
        if (JELEM	.eq.	1692	)		coeff = 	1.027325456
        if (JELEM	.eq.	1693	)		coeff = 	1.010432099
        if (JELEM	.eq.	1694	)		coeff = 	0.745892249
        if (JELEM	.eq.	1695	)		coeff = 	1.166751334
        if (JELEM	.eq.	1696	)		coeff = 	0.880688527
        if (JELEM	.eq.	1697	)		coeff = 	1.057711185
        if (JELEM	.eq.	1698	)		coeff = 	0.814344267
        if (JELEM	.eq.	1699	)		coeff = 	1.053922798
        if (JELEM	.eq.	1700	)		coeff = 	1.11009815
        if (JELEM	.eq.	1701	)		coeff = 	0.912059204
        if (JELEM	.eq.	1702	)		coeff = 	1.107610455
        if (JELEM	.eq.	1703	)		coeff = 	0.829579038
        if (JELEM	.eq.	1704	)		coeff = 	0.751859238
        if (JELEM	.eq.	1705	)		coeff = 	1.01906664
        if (JELEM	.eq.	1706	)		coeff = 	1.030054472
        if (JELEM	.eq.	1707	)		coeff = 	1.222645472
        if (JELEM	.eq.	1708	)		coeff = 	0.956640297
        if (JELEM	.eq.	1709	)		coeff = 	0.57277642
        if (JELEM	.eq.	1710	)		coeff = 	0.861667796
        if (JELEM	.eq.	1711	)		coeff = 	1.019695454
        if (JELEM	.eq.	1712	)		coeff = 	0.919309787
        if (JELEM	.eq.	1713	)		coeff = 	0.839334871
        if (JELEM	.eq.	1714	)		coeff = 	1.103165224
        if (JELEM	.eq.	1715	)		coeff = 	0.927123458
        if (JELEM	.eq.	1716	)		coeff = 	0.781290518
        if (JELEM	.eq.	1717	)		coeff = 	0.883544616
        if (JELEM	.eq.	1718	)		coeff = 	1.044329403
        if (JELEM	.eq.	1719	)		coeff = 	1.008591528
        if (JELEM	.eq.	1720	)		coeff = 	0.939341836
        if (JELEM	.eq.	1721	)		coeff = 	1.131907161
        if (JELEM	.eq.	1722	)		coeff = 	0.985608747
        if (JELEM	.eq.	1723	)		coeff = 	1.046411442
        if (JELEM	.eq.	1724	)		coeff = 	0.959007875
        if (JELEM	.eq.	1725	)		coeff = 	0.825473823
        if (JELEM	.eq.	1726	)		coeff = 	0.996618167
        if (JELEM	.eq.	1727	)		coeff = 	0.948994222
        if (JELEM	.eq.	1728	)		coeff = 	1.113193701
        if (JELEM	.eq.	1729	)		coeff = 	0.745400581
        if (JELEM	.eq.	1730	)		coeff = 	1.058172296
        if (JELEM	.eq.	1731	)		coeff = 	0.995190538
        if (JELEM	.eq.	1732	)		coeff = 	0.895917253
        if (JELEM	.eq.	1733	)		coeff = 	0.851584022
        if (JELEM	.eq.	1734	)		coeff = 	1.072767565
        if (JELEM	.eq.	1735	)		coeff = 	0.949202503
        if (JELEM	.eq.	1736	)		coeff = 	0.511324982
        if (JELEM	.eq.	1737	)		coeff = 	0.87571546
        if (JELEM	.eq.	1738	)		coeff = 	1.122396348
        if (JELEM	.eq.	1739	)		coeff = 	1.153735093
        if (JELEM	.eq.	1740	)		coeff = 	0.958259802
        if (JELEM	.eq.	1741	)		coeff = 	0.973416288
        if (JELEM	.eq.	1742	)		coeff = 	0.959752801
        if (JELEM	.eq.	1743	)		coeff = 	1.075972307
        if (JELEM	.eq.	1744	)		coeff = 	1.00245742
        if (JELEM	.eq.	1745	)		coeff = 	0.772647063
        if (JELEM	.eq.	1746	)		coeff = 	0.85856283
        if (JELEM	.eq.	1747	)		coeff = 	0.983477765
        if (JELEM	.eq.	1748	)		coeff = 	0.74725772
        if (JELEM	.eq.	1749	)		coeff = 	1.113375408
        if (JELEM	.eq.	1750	)		coeff = 	1.069162535
        if (JELEM	.eq.	1751	)		coeff = 	0.695726875
        if (JELEM	.eq.	1752	)		coeff = 	1.109372157
        if (JELEM	.eq.	1753	)		coeff = 	0.993170439
        if (JELEM	.eq.	1754	)		coeff = 	1.02731072
        if (JELEM	.eq.	1755	)		coeff = 	1.127734667
        if (JELEM	.eq.	1756	)		coeff = 	1.063526681
        if (JELEM	.eq.	1757	)		coeff = 	0.983237069
        if (JELEM	.eq.	1758	)		coeff = 	1.025988347
        if (JELEM	.eq.	1759	)		coeff = 	0.878148181
        if (JELEM	.eq.	1760	)		coeff = 	1.045569466
        if (JELEM	.eq.	1761	)		coeff = 	0.671543236
        if (JELEM	.eq.	1762	)		coeff = 	1.021544634
        if (JELEM	.eq.	1763	)		coeff = 	1.081013798
        if (JELEM	.eq.	1764	)		coeff = 	0.935751005
        if (JELEM	.eq.	1765	)		coeff = 	0.834432899
        if (JELEM	.eq.	1766	)		coeff = 	0.974336791
        if (JELEM	.eq.	1767	)		coeff = 	0.504402611
        if (JELEM	.eq.	1768	)		coeff = 	0.976796605
        if (JELEM	.eq.	1769	)		coeff = 	1.021904698
        if (JELEM	.eq.	1770	)		coeff = 	1.042988755
        if (JELEM	.eq.	1771	)		coeff = 	0.668662379
        if (JELEM	.eq.	1772	)		coeff = 	1.052087437
        if (JELEM	.eq.	1773	)		coeff = 	0.990613402
        if (JELEM	.eq.	1774	)		coeff = 	1.051442371
        if (JELEM	.eq.	1775	)		coeff = 	1.034254667
        if (JELEM	.eq.	1776	)		coeff = 	0.982914789
        if (JELEM	.eq.	1777	)		coeff = 	0.467112564
        if (JELEM	.eq.	1778	)		coeff = 	0.743301074
        if (JELEM	.eq.	1779	)		coeff = 	0.937567082
        if (JELEM	.eq.	1780	)		coeff = 	0.920833703
        if (JELEM	.eq.	1781	)		coeff = 	1.123385283
        if (JELEM	.eq.	1782	)		coeff = 	0.928109042
        if (JELEM	.eq.	1783	)		coeff = 	0.756610164
        if (JELEM	.eq.	1784	)		coeff = 	0.967756365
        if (JELEM	.eq.	1785	)		coeff = 	0.989672798
        if (JELEM	.eq.	1786	)		coeff = 	1.107230642
        if (JELEM	.eq.	1787	)		coeff = 	0.899182276
        if (JELEM	.eq.	1788	)		coeff = 	0.549681094
        if (JELEM	.eq.	1789	)		coeff = 	0.888533817
        if (JELEM	.eq.	1790	)		coeff = 	0.979531804
        if (JELEM	.eq.	1791	)		coeff = 	0.898145832
        if (JELEM	.eq.	1792	)		coeff = 	0.923655501
        if (JELEM	.eq.	1793	)		coeff = 	0.858219337
        if (JELEM	.eq.	1794	)		coeff = 	0.903724886
        if (JELEM	.eq.	1795	)		coeff = 	0.953951737
        if (JELEM	.eq.	1796	)		coeff = 	0.709152278
        if (JELEM	.eq.	1797	)		coeff = 	0.917452396
        if (JELEM	.eq.	1798	)		coeff = 	0.904494928
        if (JELEM	.eq.	1799	)		coeff = 	0.877918654
        if (JELEM	.eq.	1800	)		coeff = 	0.788696554



             
             
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