!     Subroutines for computing hydrodynamic interaction between two spheres
!     Developed by Ahmad Ababaei (ahmad.ababaei@imgw.pl) under PhD supervision
!     of Professor Bogdan Rosa (bogdan.rosa@imgw.pl)
!     Institute of Meteorology and Water Management — National Research Institute
!     Podleśna 61, 01-673 Warsaw, Poland

      PROGRAM HYDROFORCES
!     -mcmodel=medium -fopenmp export OMP_NUM_THREADS=<n>
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DIMENSION fg(4,4)
      PARAMETER ( acu = 1d-9 )
      LOGICAL opp, ncl

      CALL CPU_TIME(t0)
      OPEN (1, file='output.dat')

! ============ I N P U T S ==============

        alam = 1.00d-0 ! Radius ratio
          a1 = 1d+1    ! Larger drop radius [μm]
         mur = 1d+2    ! Viscosity ratio
         ncl = .TRUE.  ! Non-continuum lubrication ON
!        ncl = .FALSE. ! Non-continuum lubrication OFF
         opp = .TRUE.  ! Orientation: opposing
!        opp = .FALSE. ! Orientation: same

! ========= L O G   D I S T R. ========== 
! Logarithmic distribution of normalized
! gap size ξ = s — 2 in JO84 notation:
      xi_min = 1d-3
      xi_max = 1d+2
      sample = 19d0
      dlt_xi = DLOG ( xi_max / xi_min ) / sample
!     dlt_la = ( 1d0 - 5d-2  ) / sample
           s = 2d0 + xi_min + 3d-16
! ======================================= 
      DO i = 1, INT(sample)+1
      CALL MAP(s,alam,al,be) ! (s,λ) —> (α,β)

! ============ M E T H O D ==============

!     CALL J1915(opp,al,be,T1,T2,acu) ! limit exists

!     CALL SJ26M61EXP(opp,al,be,F1,F2,acu)
!     CALL SJ1926IMP(opp,al,be,F1,F2,acu)

!     CALL H1937vdW(s,alam,1d1,5d-13,F1,F2) ! wrong: needs velocity

!     CALL GCB66T(al,F1,T1,acu) ! limit exists
!     CALL GCB66R(al,F1,T1,acu) ! limit exists for force only

!     CALL ON69T(al,F1,T1,acu) ! no limit
!     CALL ON69R(al,F1,T1,acu) ! no limit

!     CALL FREEROTEQUAL(opp,al,F1,acu)

!     CALL ONM70R(opp,al,be,fg,F1,F2,T1,T2,acu)
!     CALL ONM70T(opp,al,be,fg,F1,F2,T1,T2,acu)

!     CALL FREEROTATION(opp,al,be,F1,F2,acu)

!     CALL WW72(mur,al,F1,acu)

!     CALL HHS73EXP(opp,mur,mur,al,be,F1,F2,acu)
!     CALL HHS73IMP(opp,mur,mur,al,be,F1,F2,acu)

!     CALL RM74(opp,al,be,a1,F1,F2,acu)

      CALL RSD22(opp,mur,mur,al,be,a1,F1,F2,acu)

!     CALL BRI78(mur,al,F1)

!     CALL JO84XA(opp,ncl,s,alam,F1,F2,acu)
!     CALL JO84YA(opp,s,alam,F1,F2,acu)
!     CALL JO84YB(opp,s,alam,F1,F2,T1,T2,acu)
!     CALL JO84YC(opp,s,alam,T1,T2,acu) ! fix CY12
!     CALL JO84XC(opp,s,alam,T1,T2,acu)

!     CALL WAG05ISMX(opp,mur,s,alam,F1,F2)
!     CALL WAG05ISMY(opp,mur,s,alam,F1,F2)
!     CALL ROT(opp,s,alam,F1,F2) ! not done

!     CALL GMS20a(al,be,F1,F2)
!     CALL GMS20b(al,F1) ! wrong ?

! ============ O U T P U T ==============
      WRITE(1,*) s-2d0, F1!, F2, T1, T2
!     WRITE(1,*) F1, F2
      WRITE(*,*) s-2d0, F1, F2!, T1, T2
!     WRITE(1,*) alam, F1, F2
!     WRITE(*,*) alam, F1, F2, T1, T2
! ======================================= 

!     STOP

! ========= L O G   D I S T R. ==========
      s = DLOG ( s - 2d0 ) + dlt_xi
      s = DEXP ( s ) + 2d0
!     alam = alam + dlt_la
      ENDDO
! ======================================= 

      CALL CPU_TIME(t1)
      WRITE(*,*) 'Elapsed time:', t1-t0 ,'s'

      CLOSE ( 1 )
      END PROGRAM HYDROFORCES

! ======================= S U B R O U T I N E S ========================

! ============ M A P P I N G ============
      SUBROUTINE MAP(s,alam,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
!     Zinchenko explicit transform
!     Here we map the regular JO84 (s,λ)
!     geometrical setting in two sets of
!     spherical polar coordinates (JO84 Sec. 2)
!     system to spherical bipolar, known as
!     bispherical, coordinates (ε,k), 
!     to obtain two spheres interfaces, e.g.
!     (α,β) in Stimson & Jeffery (1926)
!     notation or (ξ₁, ξ₂) in Jeffery (1915)
!     notation, using formulae given by Zinchenko
!     for example in Eq. (2.2) of his paper:
!     doi.org/10.1016/0021-8928(78)90051-5.
!     The methodology to derive these, however,
!     is not straightforward, and he kindly
!     sent it to us in a personal communication.

!     Radii ratio:
!     "k" in Zinchenko (1977) notation
!     "λ = 1/k" in JO84 notation

!     Mapping formulae: 
!     (s,λ) —> (ε,k) —> (α,β)
!     They are modified according to our problem setting.
           k = 1d0 / alam                         ! size ratio
         eps = ( s - 2d0 ) * ( alam + 1d0 ) / 2d0 ! normalized clearance: ε.a₁ = gap
      coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
      coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
!     coshal = 1d0 + eps  *   ( 1d0 + k * eps / 2d0 ) / ( 1d0 + k * ( 1d0 + eps ) )
!     coshbe = 1d0 + eps * k**2 * ( 1d0 + eps / 2d0 ) / ( 1d0 + k * ( 1d0 + eps ) )
      al = DACOSH(coshal)
      be =-DACOSH(coshbe)
!     be = asinh(-k*DSINH(al))
!     WRITE(*,*) "ε, k, α, β =",eps,k,al,be

      END SUBROUTINE

! === Jeffery (1915) ===================================================
!     Jeffery, G. B. (1915). On the steady rotation of a solid of revolution in a viscous fluid. Proceedings of the London Mathematical Society, 2(1), 327-338.
      SUBROUTINE J1915(opp,xi1,xi2,G1,G2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      LOGICAL opp

      sum1o = 0d0
      sum2o = 0d0
        rel = 1d3
          i = 0
      DO WHILE ( rel .GT. acu )
          m = DBLE(i)
         IF (opp) THEN
          sum1 = sum1o + DSINH((m+1d0)*xi1-m*xi2)**(-3) + DSINH((m+1d0)*(xi1-xi2))**(-3)
          sum2 = sum2o + DSINH((m+1d0)*xi2-m*xi1)**(-3) + DSINH((m+1d0)*(xi2-xi1))**(-3)
         ELSE
          sum1 = sum1o + DSINH((m+1d0)*xi1-m*xi2)**(-3) - DSINH((m+1d0)*(xi1-xi2))**(-3)
          sum2 = sum2o + DSINH((m+1d0)*xi2-m*xi1)**(-3) - DSINH((m+1d0)*(xi2-xi1))**(-3)
         ENDIF
           rel = MAX(DABS(sum1-sum1o)/DABS(sum1),DABS(sum2-sum2o)/DABS(sum2))
         sum1o = sum1
         sum2o = sum2
             i = i + 1
      ENDDO

!     Equation (12):
      G1 = DSINH(xi1)**3 * sum1 ! normalized by —8πμa₁³Ω₁
      G2 = DSINH(xi2)**3 * sum2 ! normalized by —8πμa₂³Ω₂

      END SUBROUTINE
! ======================================================================

! === Stimson & Jeffery (1926) - Explicit Method =======================
!     Stimson, M., & Jeffery, G. B. (1926). The motion of two spheres in a viscous fluid. Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character, 111(757), 110-116.
!     Maude, A. D. (1961). End effects in a falling-sphere viscometer. British Journal of Applied Physics, 12(6), 293.
      SUBROUTINE SJ26M61EXP(opp,al,be,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      LOGICAL opp

      sum1o = 0d0
      sum2o = 0d0
      sum3o = 0d0
!     sum4o = 0d0
!     sum5o = 0d0
!     sum6o = 0d0
!     sum7o = 0d0
        rel = 1d3
          i = 1
      DO WHILE ( rel .GT. acu )
          n = DBLE(i)
          k = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)

         IF (opp) THEN  ! Maude (1961)
         sum1 = sum1o + ( A2(n,k,al,be) + B2(n,k,al,be) &
                      +   C2(n,k,al,be) + D2(n,k,al,be) &
                        )               / Delta(n,al,be)

         sum2 = sum2o + ( A2(n,k,al,be) - B2(n,k,al,be) &
                      +   C2(n,k,al,be) - D2(n,k,al,be) &
                        )               / Delta(n,al,be)

         sum3 = sum3o + lambda2(n,al) ! (2)

         ELSE ! Stimson & Jeffery (1926)
         sum1 = sum1o + ( A1(n,k,al,be) + B1(n,k,al,be) &
                      +   C1(n,k,al,be) + D1(n,k,al,be) &
                        )               / Delta(n,al,be)

         sum2 = sum2o + ( A1(n,k,al,be) - B1(n,k,al,be) &
                      +   C1(n,k,al,be) - D1(n,k,al,be) &
                        )               / Delta(n,al,be)

!        last page for equal-sized spheres:
         sum3 = sum3o + lambda1(n,al) ! (37)
!        sum4 = sum4o + C_equal(n,k,al)
!        sum5 = sum5o + C(n,k,al,be)/Delta(n,al,be)
!        sum6 = sum6o + C_equal(n,k,al)
!        sum7 = sum7o + C(n,k,al,be)/Delta(n,al,be)
         ENDIF

           rel = MAX(DABS(sum1-sum1o)/DABS(sum1),DABS(sum2-sum2o)/DABS(sum2),DABS(sum3-sum3o)/DABS(sum3))
         sum1o = sum1
         sum2o = sum2
         sum3o = sum3
!        sum4o = sum4
!        sum5o = sum5
!        sum6o = sum6
!        sum7o = sum7
             i = i + 1
      ENDDO

      F1 = 1d0/3d0 * DSINH(al) * DABS(sum1) ! (34) normalized by —6πκa₁V₁
      F2 =-1d0/3d0 * DSINH(be) * DABS(sum2) ! (35) normalized by —6πκa₂V₂
      Fe = 4d0/3d0 * DSINH(al) * DABS(sum3) ! (37) normalized by —6πκaV

      END SUBROUTINE

! ==== F U N C T I O N S :  S & J (1926) ===============================

      FUNCTION A1(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      A1 = 4d0*DEXP(-(n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al-be))
      A1 = A1 + (2d0*n+1d0)**2*DEXP(al-be)*DSINH(al-be)
      A1 = A1 + 2d0*(2d0*n-1d0)*DSINH((n+1d0/2d0)*(al-be))*DCOSH((n+1d0/2d0)*(al+be))
      A1 = A1 - 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(al-be))*DCOSH((n-1d0/2d0)*(al+be))
      A1 = A1 - (2d0*n+1d0)*(2d0*n-1d0)*DSINH(al-be)*DCOSH(al+be)
      A1 = A1 * (2d0*n+3d0) * k
      END

      FUNCTION B1(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      B1 = 2d0*(2d0*n-1d0)*DSINH((n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al+be))
      B1 = B1 - 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(al-be))*DSINH((n-1d0/2d0)*(al+be))
      B1 = B1 + (2d0*n+1d0)*(2d0*n-1d0)*DSINH(al-be)*DSINH(al+be)
      B1 = B1 *-(2d0*n+3d0) * k
      END

      FUNCTION C1(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      C1 = 4d0*DEXP(-(n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al-be))
      C1 = C1 - (2d0*n+1d0)**2*DEXP(be-al)*DSINH(al-be)
      C1 = C1 + 2d0*(2d0*n+1d0)*DSINH((n-1d0/2d0)*(al-be))*DCOSH((n+3d0/2d0)*(al+be))
      C1 = C1 - 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(al-be))*DCOSH((n+1d0/2d0)*(al+be))
      C1 = C1 + (2d0*n+1d0)*(2d0*n+3d0)*DSINH(al-be)*DCOSH(al+be)
      C1 = C1 *-(2d0*n-1d0) * k
      END

      FUNCTION D1(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D1 = 2d0*(2d0*n+1d0)*DSINH((n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al+be))
      D1 = D1 - 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al+be))
      D1 = D1 + (2d0*n+1d0)*(2d0*n+3d0)*DSINH(al-be)*DSINH(al+be)
      D1 = D1 * (2d0*n-1d0) * k
      END

      FUNCTION Delta(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      Delta = 4d0*DSINH((n+1d0/2d0)*(al-be))**2-((2d0*n+1d0)*DSINH(al-be))**2
      END

! ==== F U N C T I O N S :  M (1961) ===================================

      FUNCTION A2(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      A2 = 2d0*(2d0*n-1d0)*DSINH((n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al+be))
      A2 = A2 - 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(al-be))*DSINH((n-1d0/2d0)*(al+be))
      A2 = A2 - (2d0*n+1d0)*(2d0*n-1d0)*DSINH(al-be)*DSINH(al+be)
      A2 = A2 * (2d0*n+3d0) * k
      END

      FUNCTION B2(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
!     B2 =-4d0*DEXP((n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al-be))         !   typo?!
      B2 =-4d0*DEXP(-(n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al-be))
      B2 = B2 - (2d0*n+1d0)**2*DEXP(al-be)*DSINH(al-be)
      B2 = B2 + 2d0*(2d0*n-1d0)*DSINH((n+1d0/2d0)*(al-be))*DCOSH((n+1d0/2d0)*(al+be))
      B2 = B2 - 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(al-be))*DCOSH((n-1d0/2d0)*(al+be))
!     B2 = B2 + (2d0*n+1d0)*(2d0*n-1d0)*DSINH(al-be)*DCOSH(al-be)           !   typo?!
      B2 = B2 + (2d0*n+1d0)*(2d0*n-1d0)*DSINH(al-be)*DCOSH(al+be)
      B2 = B2 *-(2d0*n+3d0) * k
      END

      FUNCTION C2(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
!     C2 = 2d0*(2d0*n+1d0)*DSINH((n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al-be))      !   typo?!
      C2 = 2d0*(2d0*n+1d0)*DSINH((n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al+be))
!     C2 = C2 - 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al-be)) !   typo?!
      C2 = C2 - 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al+be))
      C2 = C2 - (2d0*n+1d0)*(2d0*n+3d0)*DSINH(al-be)*DSINH(al+be)
      C2 = C2 *-(2d0*n-1d0) * k
      END

      FUNCTION D2(n,k,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
!     D2 =-4d0*DEXP(-(n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al+be))        !   typo?!
      D2 =-4d0*DEXP(-(n+1d0/2d0)*(al-be))*DSINH((n+1d0/2d0)*(al-be))
      D2 = D2 + (2d0*n+1d0)**2*DEXP(be-al)*DSINH(al-be)
      D2 = D2 + 2d0*(2d0*n+1d0)*DSINH((n-1d0/2d0)*(al-be))*DCOSH((n+3d0/2d0)*(al+be))
      D2 = D2 - 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(al-be))*DCOSH((n+1d0/2d0)*(al+be))
      D2 = D2 - (2d0*n+1d0)*(2d0*n+3d0)*DSINH(al-be)*DCOSH(al+be)
      D2 = D2 * (2d0*n-1d0) * k
      END

! ==== E Q U A L - S I Z E   S P H E R E S :  S & J (1926) =============

      FUNCTION A_equal(n,k,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      A_equal = 2d0*(1d0-DEXP(-(2d0*n+1d0)*al))+(2d0*n+1d0)*(DEXP(2d0*al)-1d0)
      A_equal = A_equal / ( 2d0*DSINH((2d0*n+1d0)*al)+(2d0*n+1d0)*DSINH(2d0*al) )
      A_equal = A_equal *-(2d0*n+3d0) * k
      END

      FUNCTION C_equal(n,k,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      C_equal = 2d0*(1d0-DEXP(-(2d0*n+1d0)*al))+(2d0*n+1d0)*(1d0-DEXP(-2d0*al))
      C_equal = C_equal / ( 2d0*DSINH((2d0*n+1d0)*al)+(2d0*n+1d0)*DSINH(2d0*al) )
      C_equal = C_equal * (2d0*n-1d0) * k
      END

      FUNCTION lambda1(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      lambda1 = 4d0*DSINH(al*(n+1d0/2d0))**2-((2d0*n+1d0)*DSINH(al))**2
      lambda1 = lambda1 / ( 2d0*DSINH(al*(2d0*n+1d0))+(2d0*n+1d0)*DSINH(al*2d0) )
      lambda1 = 1d0 - lambda1
      lambda1 = lambda1 * n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)
      END

      FUNCTION lambda2(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      lambda2 = 4d0*DCOSH(al*(n+1d0/2d0))**2+((2d0*n+1d0)*DSINH(al))**2
      lambda2 = lambda2 / ( 2d0*DSINH(al*(2d0*n+1d0))-(2d0*n+1d0)*DSINH(al*2d0) )
      lambda2 = 1d0 - lambda2
      lambda2 = lambda2 * n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)
      END

! === Stimson & Jeffery (1926) - Implicit Method =======================
!     Stimson, M., & Jeffery, G. B. (1926). The motion of two spheres in a viscous fluid. Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character, 111(757), 110-116.
!     Zinchenko's suggestion:
!     Instead of using explicit Eqs. (27)—(31) for (34) & (35), we
!     numerically solve the set (26). This gives us the advantage of
!     "choosing" the velocity boundary condition such that we can have
!     the forces when the spheres are moving not only in the same direction
!     but also an "opposing" direction, that is, Maude's solution.
      SUBROUTINE SJ1926IMP(opp,al,be,F1,F2,acu)
      USE OMP_LIB
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, DIMENSION(4,5) :: T
      INTEGER nd
      LOGICAL opp

!!$OMP PARALLEL PRIVATE(id)
!     nd = OMP_GET_NUM_THREADS()
!     id = OMP_GET_THREAD_NUM()
!     IF ( id .eq. 0 ) WRITE(*,*) "Open-MP version with threads = ", nd
!!$OMP END PARALLEL

      sum1o= 0d0
      sum2o= 0d0
       rel = 1d3
         i = 1
      DO WHILE ( rel .GT. acu )
!!$OMP PARALLEL DEFAULT(PRIVATE) FIRSTPRIVATE(i) SHARED(al,be,opp,acu) REDUCTION(+:sum1,sum2)
!       nd = OMP_GET_NUM_THREADS()
!       id = OMP_GET_THREAD_NUM()
!        n = DBLE((i-1)*nd+id+1)
         n = DBLE(i)
!        PRINT*,"i,n,id+1",i,n,id+1
         k = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)
         T(1,1) = DCOSH((n-1d0/2d0)*al)
         T(1,2) = DSINH((n-1d0/2d0)*al)
         T(1,3) = DCOSH((n+3d0/2d0)*al)
         T(1,4) = DSINH((n+3d0/2d0)*al)
         T(1,5) =-k*( (2d0*n+3d0)*DEXP(-(n-1d0/2d0)*al) - (2d0*n-1d0)*DEXP(-(n+3d0/2d0)*al) )

         T(2,1) = DCOSH((n-1d0/2d0)*be)
         T(2,2) = DSINH((n-1d0/2d0)*be)
         T(2,3) = DCOSH((n+3d0/2d0)*be)
         T(2,4) = DSINH((n+3d0/2d0)*be)
         IF (opp) THEN
         T(2,5) =+k*( (2d0*n+3d0)*DEXP((n-1d0/2d0)*be) - (2d0*n-1d0)*DEXP((n+3d0/2d0)*be) ) ! opposing direction
         ELSE
         T(2,5) =-k*( (2d0*n+3d0)*DEXP((n-1d0/2d0)*be) - (2d0*n-1d0)*DEXP((n+3d0/2d0)*be) ) ! same direction
         ENDIF

         T(3,1) = (2d0*n-1d0)*DSINH((n-1d0/2d0)*al)
         T(3,2) = (2d0*n-1d0)*DCOSH((n-1d0/2d0)*al)
         T(3,3) = (2d0*n+3d0)*DSINH((n+3d0/2d0)*al)
         T(3,4) = (2d0*n+3d0)*DCOSH((n+3d0/2d0)*al)
         T(3,5) = (2d0*n-1d0)*(2d0*n+3d0)*k*( DEXP(-(n-1d0/2d0)*al) - DEXP(-(n+3d0/2d0)*al) )

         T(4,1) = (2d0*n-1d0)*DSINH((n-1d0/2d0)*be)
         T(4,2) = (2d0*n-1d0)*DCOSH((n-1d0/2d0)*be)
         T(4,3) = (2d0*n+3d0)*DSINH((n+3d0/2d0)*be)
         T(4,4) = (2d0*n+3d0)*DCOSH((n+3d0/2d0)*be)
         IF (opp) THEN
         T(4,5) =+(2d0*n-1d0)*(2d0*n+3d0)*k*( DEXP((n-1d0/2d0)*be) - DEXP((n+3d0/2d0)*be) ) ! opposing direction
         ELSE
         T(4,5) =-(2d0*n-1d0)*(2d0*n+3d0)*k*( DEXP((n-1d0/2d0)*be) - DEXP((n+3d0/2d0)*be) ) ! same direction
         ENDIF

         CALL GAUSS(4,5,T)

          sum1 = sum1o + SUM(T(1:4,5))
          sum2 = sum2o + T(1,5) - T(2,5) + T(3,5) - T(4,5)
!!$OMP END PARALLEL
           rel = MAX(DABS((sum1-sum1o)/DABS(sum1)),DABS((sum2-sum2o)/DABS(sum2)))
         sum1o = sum1
         sum2o = sum2
!        WRITE(*,*)"n_max,sum1,sum2,sum1o,sum2o,rel",i,sum1,sum2,sum1o,sum2o,rel
!        pause
             i = i + 1
      ENDDO

      F1 =-DSINH(al)/3d0 *      sum1  ! (34) normalized by —6πκa₁V₁
      F2 =-DSINH(be)/3d0 * DABS(sum2) ! (35) normalized by —6πκa₂V₂

      END SUBROUTINE
! ======================================================================

! === Hamaker (1937) - van der Waals ===================================
!     Hamaker, H. C. (1937). The London—van der Waals attraction between spherical particles. physica, 4(10), 1058-1072.
      SUBROUTINE H1937vdW(s,alam,aa,A,F1,F2)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

!     A  = 5d-13    ! Hamaker constant g.cm2/s2 
      a1 = aa / 1d4 ! radius in cm
      a2 = alam * a1
      mu = 1.7d-4   ! dyn visc of air
      pi = 4d0*DATAN(1d0)

      !     Dimensional formulation: doi.org/10.1016/S0031-8914(37)80203-7
!     r  = (a1+a2)/2d0*s
!     D1 = r**2-(a1+a2)**2
!     D2 = r**2-(a1-a2)**2
!     F  = A/6d0*(4d0*a1*a2)**3*r/D1**2/D2**2 ! d(8) / dC
!     F1 = F/(6d0*pi*mu*a1)
!     F2 = F/(6d0*pi*mu*a2)
      
!     Dimensionless formulation: doi.org/10.1017/S002211208400286X
      D1 = (s**2-4d0)*(1d0+alam)**2
      D2 = s**2*(1d0+alam)**2-4d0*(1d0-alam)**2
      F  = A/6d0*(16d0*alam)**3*s*(1d0+alam)**2/D1**2/D2**2 ! d(2.2a) / dr
      F1 = F/(6d0*pi*mu*a1)*2d0/(a1+a2) ! Note1: dPhi / dr = dPhi / ds * ds / dr
      F2 = F/(6d0*pi*mu*a2)*2d0/(a1+a2) ! Note2:   ds / dr = 2 / ( a1 + a2 )

      END SUBROUTINE
! ======================================================================

! === Goldman, Cox, Brenner (1966) - Translation =======================
!     Goldman, A. J., Cox, R. G., & Brenner, H. (1966). The slow motion of two identical arbitrarily oriented spheres through a viscous fluid. Chemical Engineering Science, 21(12), 1151-1170.

!     "Cross effects" due to translation-rotation
!     coupling: for any alpha the ratio of F1 from
!     rotating spheres subroutine to T1 of this
!     (translating spheres) subroutine is 4/3, given
!     in (4.40) in the paper

      SUBROUTINE GCB66T(al,F1,T1,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: At
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

      F1o= 0d0
      T1o= 0d0
      iN = 50
1     ALLOCATE ( T(iN,4), At(-1:iN+1) )
      At = 0d0
      DO i = 1, iN
         n = DBLE(i)

         T(i,1) = (n-1d0)*(gam(n-1d0,al)-1d0)-(n-1d0)*(2d0*n-3d0)/(2d0*n-1d0)*(gam(n,al)-1d0)
         T(i,2) = (2d0*n+1d0)-5d0*gam(n,al)-n*(2d0*n-1d0)/(2d0*n+1d0)*(gam(n-1d0,al)+1d0)   &
                + (n+1d0)*(2d0*n+3d0)/(2d0*n+1d0)*(gam(n+1d0,al)-1d0)
         T(i,3) = (n+2d0)*(2d0*n+5d0)/(2d0*n+3d0)*(gam(n,al)+1d0)-(n+2d0)*(gam(n+1d0,al)+1d0)
         T(i,4) = DSQRT(2d0)*DEXP(-(n+5d-1)*al)*( DEXP(al)/DCOSH((n-5d-1)*al)               &
                - 2d0/DCOSH((n+5d-1)*al)+DEXP(-al)/DCOSH((n+15d-1)*al) )
      ENDDO

      CALL TDMA(iN,T)

      At(1:iN) = T(1:iN,4)

      sumF = 0d0
      sumT = 0d0
      DO i = 0, iN
         n = DBLE(i)
         Bt_n = 2d0*(n-1d0)/(2d0*n-1d0)*(gam(n,al)-1d0)*At(i-1)     &
              - 2d0*gam(n,al)*At(i) + 2d0*(n+2d0)/(2d0*n+3d0)       &
              * (gam(n,al)+1d0)*At(i+1)

         Dt_n = DSQRT(8d0)*DEXP(-(n+5d-1)*al)/DCOSH((n+5d-1)*al)    &
              - n*(n-1d0)/(2d0*n-1d0)*(gam(n,al)-1d0)      *At(i-1) &
              + (n+1d0)*(n+2d0)/(2d0*n+3d0)*(gam(n,al)+1d0)*At(i+1)
         sumF = sumF + Dt_n + n*(n+1d0) * Bt_n
         sumT = sumT + 2d0*n* (n+1d0)*(2d0+DEXP(-(2d0*n+1d0)*al))   &
                            * At(i) + (2d0-DEXP(-(2d0*n+1d0)*al))   &
              * ( n*(n+1d0) * Bt_n / DTANH(al) -(2d0*n+1d0-1d0      &
              /  DTANH(al)) * Dt_n )
      ENDDO
      DEALLOCATE ( T, At )

      F1 = DSQRT(2d0)/6d0 * DSINH(al) * sumF ! (3.57) normalized by —6πμaU
      T1 =-DSINH(al)**2/DSQRT(288d0)  * sumT ! (3.58) normalized by —8πμa²U

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(T1-T1o)/DABS(T1))
      IF ( rel .GT. acu ) THEN
         iN = INT(2.0 * FLOAT(iN)) ! 100% increase
         F1o= F1
         T1o= T1
!        WRITE(*,*) "n_max, F1, T1 = ", iN,F1,T1
         GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === Goldman, Cox, Brenner (1966) - Rotation ==========================
!     Goldman, A. J., Cox, R. G., & Brenner, H. (1966). The slow motion of two identical arbitrarily oriented spheres through a viscous fluid. Chemical Engineering Science, 21(12), 1151-1170.

!     "Cross effects" due to translation-rotation
!     coupling: for any alpha the ratio of F1 from this
!     (rotating spheres) subroutine to T1 of the
!     translating spheres subroutine is 4/3, given
!     in (4.40) in the paper

      SUBROUTINE GCB66R(al,F1,T1,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: Ar
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

      F1o= 0d0
      T1o= 0d0
      iN = 50
1     ALLOCATE ( T(iN,4), Ar(-1:iN+1) )
      Ar = 0d0
      DO i = 1, iN
         n = DBLE(i)

         T(i,1) = (n-1d0)*(gam(n-1d0,al)-1d0)-(n-1d0)*(2d0*n-3d0)/(2d0*n-1d0)*(gam(n,al)-1d0)
         T(i,2) = (2d0*n+1d0)-5d0*gam(n,al)-n*(2d0*n-1d0)/(2d0*n+1d0)*(gam(n-1d0,al)+1d0)   &
                + (n+1d0)*(2d0*n+3d0)/(2d0*n+1d0)*(gam(n+1d0,al)-1d0)
         T(i,3) = (n+2d0)*(2d0*n+5d0)/(2d0*n+3d0)*(gam(n,al)+1d0)-(n+2d0)*(gam(n+1d0,al)+1d0)
         T(i,4) = DSQRT(2d0) * DEXP(-(n+5d-1)*al) / DSINH(al)   *                           &
                ( (2d0*n+1d0)/ DCOSH((n+5d-1)*al) *                                         &
                ( DEXP(al)/(2d0*n-1d0) + DEXP(-al)/(2d0*n+3d0)  )                           &
                - (2d0*n-1d0)/(2d0*n+1d0) / DCOSH((n- 5d-1)*al)                             &
                - (2d0*n+3d0)/(2d0*n+1d0) / DCOSH((n+15d-1)*al) )
      ENDDO

      CALL TDMA(iN,T)

      Ar(1:iN) = T(1:iN,4)

      sumF = 0d0
      sumT = 0d0
      DO i = 0, iN
         n = DBLE(i)
         Br_n = 4d0*tau(n,al)/DSINH(al)/DCOSH((n+5d-1)*al)          &
              + 2d0*(n-1d0)/(2d0*n-1d0)*(gam(n,al)-1d0)*Ar(i-1)     &
              - 2d0*gam(n,al)*Ar(i) + 2d0*(n+2d0)/(2d0*n+3d0)       &
              * (gam(n,al)+1d0)*Ar(i+1)

         Dr_n = DSQRT(8d0)/DSINH(al) / DCOSH((n+05d-1)*al)          &
              * (   n**2  /(2d0*n-1d0)*DEXP(-(n-05d-1)*al)          &
              - (n+1d0)**2/(2d0*n+3d0)*DEXP(-(n+15d-1)*al)  )       &
              -  n*(n-1d0)/(2d0*n-1d0)*(gam(n,al)-1d0)     *Ar(i-1) &
              + (n+1d0)*(n+2d0)/(2d0*n+3d0)*(gam(n,al)+1d0)*Ar(i+1)

         sumF = sumF + Dr_n + n*(n+1d0) * Br_n
         sumT = sumT + 2d0*n*(n+1d0)*(2d0+DEXP(-(2d0*n+1d0)*al))    &
                            * Ar(i) +(2d0-DEXP(-(2d0*n+1d0)*al))    &
              * ( n*(n+1d0) * Br_n  / DTANH(al)-(2d0*n+1d0-1d0      &
              /  DTANH(al)) * Dr_n  )
      ENDDO
      DEALLOCATE ( T, Ar )

      F1 = DSQRT(2d0)/6d0*DSINH(al)**2       * sumF ! (4.24) normalized by —6πμa²Ω
      T1 = 1d0/3d0-DSINH(al)**3/DSQRT(288d0) * sumT ! (4.25) normalized by —8πμa³Ω

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(T1-T1o)/DABS(T1))
      IF ( rel .GT. acu ) THEN
         iN = INT(2.0 * FLOAT(iN)) ! 100% increase
         F1o= F1
         T1o= T1
!        WRITE(*,*) "n_max, F1, T1 = ", iN,F1,T1
         GOTO 1
      ENDIF

      END SUBROUTINE

! === F U N C T I O N S :  G C B (1966) ================================

      FUNCTION gam(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      gam = DTANH((n+5d-1)*al)/DTANH(al)
      END

      FUNCTION tau(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DOUBLE PRECISION :: n,al
      tau = (DEXP(-(n+15d-1)*al)/(2d0*n+3d0)-DEXP(-(n-5d-1)*al)/(2d0*n-1d0))/DSQRT(2d0)
      END
! ======================================================================

! === O'Neill (1969) - Rotation ========================================
!     O'Neill, M. E. (1969). Exact solutions of the equations of slow viscous flow generated by the asymmetrical motion of two equal spheres. Applied Scientific Research, 21(1), 452-466.
      SUBROUTINE ON69R(al,f1,g1,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: An
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

      f1o= 0d0
      g1o= 0d0
      iN = 50
1     ALLOCATE ( T(iN,4), An(-1:iN+1) )
      An = 0d0
      DO i = 1, iN
         n = DBLE(i)
          fctr1 = (2d0*n-3d0)*gma(n,al)-(2d0*n-1d0)*gma(n-1d0,al)+2d0
          fctr2 = (2d0*n+3d0)*gma(n+1d0,al)-(2d0*n+5d0)*gma(n,al)-2d0
         T(i,1) = fctr1 * (n-1d0)/(2d0*n-1d0)
         T(i,2) = fctr1 * -n/(2d0*n+1d0) + fctr2 * -(n+1d0)/(2d0*n+1d0)
         T(i,3) = fctr2 * (n+2d0)/(2d0*n+3d0)
         T(i,4) = DSQRT(2d0)/DCOSH(al) * &
                ( (gma(n-1d0,al)-1d0/DTANH(al))*(2d0*n-1d0)/(2d0*n+1d0)*DEXP(-al) &
                - (gma(n,al)-1d0/DTANH(al))*((2d0*n+1d0)/(2d0*n-1d0)*DEXP(al)+(2d0*n+1d0)/(2d0*n+3d0)*DEXP(-al)) &
                + (gma(n+1d0,al)-1d0/DTANH(al))*(2d0*n+3d0)/(2d0*n+1d0)*DEXP(al) )
      ENDDO

      CALL TDMA(iN,T)

      An(1:iN) = T(1:iN,4)

      sumf = 0d0
      sumg = 0d0
      DO i = 0, iN
         n = DBLE(i)
         D_n = (n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1)*(gma(n,al)+1d0) &
             -  n     *(n-1d0)/(2d0*n-1d0)*An(i-1)*(gma(n,al)-1d0) &
             + DSQRT(8d0) / DCOSH(al) * (gma(n,al)-1d0/DTANH(al))  &
             * ( n**2*DEXP(al)/(2d0*n-1d0)-(n+1d0)**2*DEXP(-al)/(2d0*n+3d0) )
        sumf = sumf + D_n
        sumg = sumg + D_n * (2d0*n+1d0-1d0/DTANH(al))
      ENDDO
      DEALLOCATE ( T, An )

      f1 =-DSQRT(2d0)/3d0 * DSINH(al)**2 * sumf ! (5.1) normalized by —6πμa²Ω
      g1 = DSQRT(2d0)/4d0 * DSINH(al)**3 * sumg ! (5.2) normalized by —8πμa³Ω

      rel = MAX(DABS(f1-f1o)/DABS(f1),DABS(g1-g1o)/DABS(g1))
      IF ( rel .GT. acu ) THEN
         iN = INT(2.0 * FLOAT(iN)) ! 100% increase
        f1o = f1
        g1o = g1
!       WRITE(*,*) "n_max, f1, g1 = ", iN,f1,g1
        GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === O'Neill (1969) - Translation =====================================
!     O'Neill, M. E. (1969). Exact solutions of the equations of slow viscous flow generated by the asymmetrical motion of two equal spheres. Applied Scientific Research, 21(1), 452-466.
      SUBROUTINE ON69T(al,f2,g2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: An
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

      f2o= 0d0
      g2o= 0d0
      iN = 50
1     ALLOCATE ( T(iN,4), An(-1:iN+1) )
      An = 0d0
      DO i = 1, iN
         n = DBLE(i)
          fctr1 = (2d0*n-3d0)*gma(n,al)-(2d0*n-1d0)*gma(n-1d0,al)+2d0
          fctr2 = (2d0*n+3d0)*gma(n+1d0,al)-(2d0*n+5d0)*gma(n,al)-2d0
         T(i,1) = fctr1 * (n-1d0)/(2d0*n-1d0)
         T(i,2) = fctr1 * -n/(2d0*n+1d0) + fctr2 * -(n+1d0)/(2d0*n+1d0)
         T(i,3) = fctr2 * (n+2d0)/(2d0*n+3d0)
         T(i,4) = DSQRT(2d0)*DTANH(al)*(2d0*gma(n,al)-gma(n-1d0,al)-gma(n+1d0,al))
      ENDDO

      CALL TDMA(iN,T)

      An(1:iN) = T(1:iN,4)

      sumf = 0d0
      sumg = 0d0
      DO i = 0, iN
         n = DBLE(i)
         D_n = (n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1)*(gma(n,al)+1d0) &
             -  n     *(n-1d0)/(2d0*n-1d0)*An(i-1)*(gma(n,al)-1d0) &
             + DSQRT(8d0)*(DTANH((n+5d-1)*al)**(-1)-1d0)
        sumf = sumf + D_n
        sumg = sumg + D_n * (2d0*n+1d0-1d0/DTANH(al))
      ENDDO
      DEALLOCATE ( T, An )

      f2 = DSQRT(2d0)/3d0 * DSINH(al)    * sumf ! (5.3) normalized by —6πμaU
      g2 =-DSQRT(2d0)/4d0 * DSINH(al)**2 * sumg ! (5.4) normalized by —8πμa²U

      rel = MAX(DABS(f2-f2o)/DABS(f2),DABS(g2-g2o)/DABS(g2))
      IF ( rel .GT. acu ) THEN
         iN = INT(2.0 * FLOAT(iN)) ! 100% increase
        f2o = f2
        g2o = g2
!       WRITE(*,*) "n_max, f2, g2 = ", iN,f2,g2,rel
        GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === Forces acting on a pair of freely rotating equal-sized spheres ===

      SUBROUTINE FREEROTEQUAL(opp,al,F,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      LOGICAL opp

      IF (opp) THEN
          CALL ON69T(al,F1,T1,acu)
          F_tra = F1
          T_tra = T1

          CALL ON69R(al,F1,T1,acu)
          F_rot = F1
          T_rot = T1
      ELSE
          CALL GCB66T(al,F1,T1,acu)
          F_tra = F1
          T_tra = T1

          CALL GCB66R(al,F1,T1,acu)
          F_rot = F1
          T_rot = T1
      ENDIF

      F = F_tra - F_rot * T_tra / T_rot

      END SUBROUTINE
! ======================================================================

! === O'Neill & Majumdar (1970) - Rotation =============================
!     O'Neill, M. E., & Majumdar, R. (1970). Asymmetrical slow viscous fluid motions caused by the translation or rotation of two spheres. Part I: The determination of exact solutions for any values of the ratio of radii and separation parameters. Zeitschrift für angewandte Mathematik und Physik ZAMP, 21(2), 164-179.
      SUBROUTINE ONM70R(opp,al,be,fg,F1,F2,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: An, Bn
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: T
      DOUBLE PRECISION :: fg(4,4)
      LOGICAL opp

      F1o= 0d0
      F2o= 0d0
      T1o= 0d0
      T2o= 0d0
      iN = 25
1     ALLOCATE ( T(2*iN,8), An(-1:iN+1), Bn(-1:iN+1) )
       j = 1 ! first droplet
2      T = 0d0
      DO i = 1, iN
         n = DBLE(i)
         fctr1 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)+2d0*kappa(al,be)
         fctr2 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)+2d0*kappa(al,be)
         fctr3 = (2d0*n-1d0)*(gamm(n-1d0,al,be)-dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)-dlta(n,al,be))-4d0
         fctr4 = (2d0*n+5d0)*(gamm(n,al,be)-dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)-dlta(n+1d0,al,be))+4d0
         IF (i.GT.1) THEN
         T(2*i-1,2) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (39)
         T(2*i-1,3) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (39)
         ENDIF
         T(2*i-1,4) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (39)
         T(2*i-1,5) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (39)
         IF (i.LT.iN) THEN
         T(2*i-1,6) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (39)
         T(2*i-1,7) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (39)
         ENDIF
         T(2*i-1,8) = D_1R(n,al,be)

         fctr1 = (2d0*n-1d0)*(gamm(n-1d0,al,be)+dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)+dlta(n,al,be))+4d0
         fctr2 = (2d0*n+5d0)*(gamm(n,al,be)+dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)+dlta(n+1d0,al,be))-4d0
         fctr3 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)-2d0*kappa(al,be)
         fctr4 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)-2d0*kappa(al,be)

         IF (i.GT.1) THEN
         T(2*i,1) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (40)
         T(2*i,2) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (40)
         ENDIF
         T(2*i,3) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (40)
         T(2*i,4) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (40)
         IF (i.LT.iN) THEN
         T(2*i,5) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (40)
         T(2*i,6) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (40)
         ENDIF
         T(2*i,8) = D_2R(n,al,be)
      ENDDO

      CALL THOMAS(2*iN,3,3,T)

      An = 0d0
      Bn = 0d0
      DO i = 1, iN
         An(i) = T(2*i  ,8)
         Bn(i) = T(2*i-1,8)
      ENDDO

       f11 = 0d0
       g11 = 0d0
       f12 = 0d0
       g12 = 0d0
      DO i = 0, iN
         n = DBLE(i)
        En =-(kappa(al,be)+alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(gamm(n,al,be) &
           - dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*Bn(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1)+(lambda(n,al)/DSINH(al) &
           - DSQRT(2d0)*(2d0*n+1d0)*DEXP(-(n+5d-1)*al))*DSINH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

        Fn = (gamm(n,al,be)+dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(kappa(al,be) &
           - alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*An(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1)-(lambda(n,al)/DSINH(al) &
           - DSQRT(2d0)*(2d0*n+1d0)*DEXP(-(n+5d-1)*al))*DCOSH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

       f11 = f11 + ( En + Fn )
       f12 = f12 + ( En - Fn )
       g11 = g11 + ( En + Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH( al) )
       g12 = g12 + ( En - Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH(-be) )
      ENDDO

      f11 = DSQRT(2d0)/3d0*DSINH( al)**2*f11 ! (5.7)
      f12 = DSQRT(2d0)/3d0*DSINH(-be)**2*f12 ! (5.9)
      g11 = DSQRT(2d0)/4d0*DSINH( al)**3*g11 ! (5.8)
      g12 = DSQRT(2d0)/4d0*DSINH(-be)**3*g12 ! (5.10) typo !

!     IF ( j .EQ. 1 ) WRITE(*,*)"f11,f12,g11,g12=",f11,f12,g11,g12
!     IF ( j .EQ. 2 ) WRITE(*,*)"f12,f11,g12,g11=",f12,f11,g12,g11

      IF ( j .EQ. 1 ) THEN
       fg(1,1) = f11
       fg(1,2) = f12
       fg(1,3) = g11
       fg(1,4) = g12
           tmp =-al
            al =-be
            be = tmp
             j = 2 ! second droplet
            GOTO 2
      ENDIF
      DEALLOCATE ( T, An, Bn )
      fg(2,1) = f12
      fg(2,2) = f11
      fg(2,3) = g12
      fg(2,4) = g11
           tmp=-al
           al =-be
           be = tmp

!     Normalization: —6πμa²Ω  &  —8πμa³Ω
      IF (opp) THEN
         F1 =-fg(1,1) - fg(2,1)
         F2 = fg(1,2) + fg(2,2)
         T1 = fg(1,3) + fg(2,3)
         T2 = fg(1,4) + fg(2,4)
      ELSE
         F1 =-fg(1,1) + fg(2,1)
         F2 =-fg(1,2) + fg(2,2)
         T1 = fg(1,3) - fg(2,3)
         T2 =-fg(1,4) + fg(2,4)
      ENDIF

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2),DABS(T1-T1o)/DABS(T1),DABS(T2-T2o)/DABS(T2))
      IF ( rel .GT. acu ) THEN
         iN = INT(1.5 * FLOAT(iN)) ! 50% increase
        F1o = F1
        F2o = F2
        T1o = T1
        T2o = T2
!       WRITE(*,*) "n_max, F1, F2, T1, T2 = ", iN,F1,F2,T1,T2,rel
        GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === O'Neill & Majumdar (1970) - Translation ==========================
!     O'Neill, M. E., & Majumdar, R. (1970). Asymmetrical slow viscous fluid motions caused by the translation or rotation of two spheres. Part I: The determination of exact solutions for any values of the ratio of radii and separation parameters. Zeitschrift für angewandte Mathematik und Physik ZAMP, 21(2), 164-179.
      SUBROUTINE ONM70T(opp,al,be,fg,F1,F2,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: An, Bn
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: T
      DOUBLE PRECISION :: fg(4,4)
      LOGICAL opp

      F1o= 0d0
      F2o= 0d0
      T1o= 0d0
      T2o= 0d0
      iN = 25
1     ALLOCATE ( T(2*iN,8), An(-1:iN+1), Bn(-1:iN+1) )
       j = 1 ! first droplet
2      T = 0d0

      DO i = 1, iN
         n = DBLE(i)
         fctr1 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)+2d0*kappa(al,be)
         fctr2 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)+2d0*kappa(al,be)
         fctr3 = (2d0*n-1d0)*(gamm(n-1d0,al,be)-dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)-dlta(n,al,be))-4d0
         fctr4 = (2d0*n+5d0)*(gamm(n,al,be)-dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)-dlta(n+1d0,al,be))+4d0
         IF (i.GT.1) THEN
         T(2*i-1,2) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (39)
         T(2*i-1,3) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (39)
         ENDIF
         T(2*i-1,4) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (39)
         T(2*i-1,5) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (39)
         IF (i.LT.iN) THEN
         T(2*i-1,6) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (39)
         T(2*i-1,7) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (39)
         ENDIF
         T(2*i-1,8) = D_1T(n,al,be)

         fctr1 = (2d0*n-1d0)*(gamm(n-1d0,al,be)+dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)+dlta(n,al,be))+4d0
         fctr2 = (2d0*n+5d0)*(gamm(n,al,be)+dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)+dlta(n+1d0,al,be))-4d0
         fctr3 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)-2d0*kappa(al,be)
         fctr4 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)-2d0*kappa(al,be)
         IF (i.GT.1) THEN
         T(2*i,1) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (40)
         T(2*i,2) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (40)
         ENDIF
         T(2*i,3) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (40)
         T(2*i,4) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (40)
         IF (i.LT.iN) THEN
         T(2*i,5) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (40)
         T(2*i,6) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (40)
         ENDIF
         T(2*i,8) = D_2T(n,al,be)
      ENDDO

      CALL THOMAS(2*iN,3,3,T)

      An = 0d0
      Bn = 0d0
      DO i = 1, iN
         An(i) = T(2*i  ,8)
         Bn(i) = T(2*i-1,8)
      ENDDO

       f21 = 0d0
       g21 = 0d0
       f22 = 0d0
       g22 = 0d0
      DO i = 0, iN
         n = DBLE(i)
        En =-(kappa(al,be)+alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(gamm(n,al,be) &
           - dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*Bn(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1)-DSQRT(8d0)*DEXP(-(n+5d-1)*al) &
           * DSINH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

        Fn = (gamm(n,al,be)+dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(kappa(al,be) &
           - alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*An(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1)+DSQRT(8d0)*DEXP(-(n+5d-1)*al) &
           * DCOSH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

       f21 = f21 + ( En + Fn )
       f22 = f22 + ( En - Fn )
       g21 = g21 + ( En + Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH( al) )
       g22 = g22 + ( En - Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH(-be) )
      ENDDO

      f21 = DSQRT(2d0)/3d0*DSINH( al)  * f21 ! (5.11)
      f22 = DSQRT(2d0)/3d0*DSINH(-be)  * f22 ! (5.11)
      g21 = DSQRT(2d0)/4d0*DSINH( al)**2*g21 ! (5.12)
      g22 = DSQRT(2d0)/4d0*DSINH(-be)**2*g22 ! (5.13)

!     IF ( j .EQ. 1 ) WRITE(*,*)"f21,f22,g21,g22=",f21,f22,g21,g22
!     IF ( j .EQ. 2 ) WRITE(*,*)"f22,f21,g22,g21=",f22,f21,g22,g21

      IF ( j .EQ. 1 ) THEN ! first droplet
         fg(3,1) = f21
         fg(3,2) = f22
         fg(3,3) = g21
         fg(3,4) = g22
         IF (opp) THEN
         F1 = f21
         F2 =-f22
         T1 =-g21
         T2 =-g22
         ELSE
         F1 = f21
         F2 = f22
         T1 =-g21
         T2 = g22
         ENDIF
        tmp =-al
         al =-be
         be = tmp
          j = 2 ! second droplet
         GOTO 2
      ENDIF
      DEALLOCATE ( T, An, Bn )
      tmp=-al
      al =-be
      be = tmp
      fg(4,1) = f22
      fg(4,2) = f21
      fg(4,3) = g22
      fg(4,4) = g21

!     Normalization: —6πμaV  &  —8πμa²V
      IF (opp) THEN
         F1 = F1 - f22
         F2 = F2 + f21
         T1 = T1 + g22
         T2 = T2 + g21
      ELSE
         F1 = F1 + f22
         F2 = F2 + f21
         T1 = T1 - g22
         T2 = T2 + g21
      ENDIF

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2),DABS(T1-T1o)/DABS(T1),DABS(T2-T2o)/DABS(T2))
      IF ( rel .GT. acu ) THEN
         iN = INT(1.3 * FLOAT(iN)) ! 30% increase
        F1o = F1
        F2o = F2
        T1o = T1
        T2o = T2
!       WRITE(*,*) "n_max, F1, F2, T1, T2 = ", iN,F1,F2,T1,T2,rel
        GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === F U N C T I O N S :  O M (1970) ==================================

      FUNCTION alpha(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      alpha = DSINH(al-be)/DSINH(al)/DSINH(be)/DSINH((n+5d-1)*(al-be))
      alpha = alpha * DSINH((n+5d-1)*(al+be))
      END

      FUNCTION gamm(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      gamm = DSINH(al-be)/DSINH(al)/DSINH(be)/DSINH((n+5d-1)*(al-be))
      gamm = gamm * DCOSH((n+5d-1)*(al+be))
      END

      FUNCTION dlta(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      dlta = DSINH(al-be)/DSINH(al)/DSINH(be)/DSINH((n+5d-1)*(al-be))
      dlta = dlta * DCOSH((n+5d-1)*(al-be))
      END

      FUNCTION kappa(al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      kappa = DSINH(al+be)/DSINH(al)/DSINH(be)
      END

      FUNCTION lambda(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      lambda = DEXP(-(n+5d-1)*al)/DSQRT(2d0)
      lambda = lambda * ( DEXP(-al)/(2d0*n+3d0) - DEXP(al)/(2d0*n-1d0) )
      END

      FUNCTION D_1R(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_1R =-(2d0*n+1d0)**2*DSINH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))
      D_1R = D_1R * ( DEXP(al)/(2d0*n-1d0) + DEXP(-al)/(2d0*n+3d0) )
      D_1R = D_1R + (2d0*n-1d0)*DSINH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_1R = D_1R + (2d0*n+3d0)*DSINH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_1R = D_1R * DSQRT(8d0) *DEXP(-(n+05d-1)*al)/(2d0*n+1d0)/DSINH(al)
      END

      FUNCTION D_2R(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_2R =-(2d0*n+1d0)**2*DCOSH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))
      D_2R = D_2R * ( DEXP(al)/(2d0*n-1d0) + DEXP(-al)/(2d0*n+3d0) )
      D_2R = D_2R + (2d0*n-1d0)*DCOSH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_2R = D_2R + (2d0*n+3d0)*DCOSH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_2R = D_2R * DSQRT(8d0) *DEXP(-(n+05d-1)*al)/(2d0*n+1d0)/DSINH(al)
      END

      FUNCTION D_1T(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_1T =  2d0 * DEXP(-(n+05d-1)*al)*DSINH((n+05d-1)*be)/DSINH((n+05d-1)*(al-be))
      D_1T = D_1T - DEXP(-(n-05d-1)*al)*DSINH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_1T = D_1T - DEXP(-(n+15d-1)*al)*DSINH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_1T = D_1T * DSQRT(8d0)
      END

      FUNCTION D_2T(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_2T =  2d0 * DEXP(-(n+05d-1)*al)*DCOSH((n+05d-1)*be)/DSINH((n+05d-1)*(al-be))
      D_2T = D_2T - DEXP(-(n-05d-1)*al)*DCOSH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_2T = D_2T - DEXP(-(n+15d-1)*al)*DCOSH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_2T = D_2T * DSQRT(8d0)
      END

! ======================================================================

! === Forces acting on a pair of freely rotating unequal-sized spheres =
!     In O'Neill & Majumdar (1970), by setting the torques equal to zero,
!     i.e. (5.16) = (5.17) = 0, we allow the particles to freely rotate.
!     This leads to a system of equations that describes the angular
!     velocities as functions of translational velocities, which then
!     could be used in (5.14) and (5.15) to compute torque-free forces
!     acting on two spheres translating normal to their line of centers.

!     y: out normal to this plane                 ↺ Ω₁
!     ⨀ → x                                    a₁ ◯ → V₁
!     ↓ 
!     z                                           ↺ Ω₂
!                                              a₂ ◯ → V₂

!                        T R A N S L A T I O N                 R O T A T I O N
!                     ---------------------------     ---------------------------------
!     F₁-x = —6πμa₁ { f21(α,β) V₁ + f22(-β,-α) V₂  -  f11(α,β) a₁ Ω₁ + f12(-β,-α) a₁ Ω₂  } → free rotation
!                     > 0           < 0               < 0              > 0
!     F₂-x = —6πμa₂ { f22(α,β) V₁ + f21(-β,-α) V₂  -  f12(α,β) a₂ Ω₁ + f11(-β,-α) a₂ Ω₂  } → free rotation
!                     < 0           > 0               > 0              < 0
!     T₁-y = —8πμa₁²{-g21(α,β) V₁ - g22(-β,-α) V₂  +  g11(α,β) a₁ Ω₁ - g12(-β,-α) a₁ Ω₂  } = 0
!                     < 0           > 0               > 0              < 0
!     T₂-y = —8πμa₂²{ g22(α,β) V₁ + g21(-β,-α) V₂  -  g12(α,β) a₂ Ω₁ + g11(-β,-α) a₂ Ω₂  } = 0
!                     > 0           < 0               < 0              > 0
 
!     To verify the signs used above we consider two particles with
!     translational velocities along x+ and rotational velocities
!     around y+ and the resistances are assumed to have the same
!     fixed signs (there are exceptions) as in O'Neill & Majumdar.
!     Then the rest of the signs can be inferred from the direction
!     of force/torque induced by each disturbance component V₁,V₂,Ω₁,Ω₂.

      SUBROUTINE FREEROTATION(opp,al,be,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION :: fg(4,4), Om(2,3)
      LOGICAL opp

      k = DSINH(al) / DSINH(-be)
      CALL ONM70R(opp,al,be,fg,F1,F2,T1,T2,acu)
      CALL ONM70T(opp,al,be,fg,F1,F2,T1,T2,acu)

      Om(1,1) = fg(1,3)
      Om(2,1) =-fg(1,4)*k
      Om(1,2) =-fg(2,3)/k
      Om(2,2) = fg(2,4)

      IF (opp) THEN
         Om(1,3) = fg(3,3) - fg(4,3)
         Om(2,3) =-fg(3,4) + fg(4,4)
      ELSE
         Om(1,3) = fg(3,3) + fg(4,3)
         Om(2,3) =-fg(3,4) - fg(4,4)
      ENDIF

      CALL GAUSS(2,3,Om)
!     WRITE(*,*) 'Free Rotation a x Ω', Om(:,3)

      IF (opp) THEN
         F1 = fg(3,1) - fg(4,1) - fg(1,1) * Om(1,3) + fg(2,1) * Om(2,3) / k
         F2 = fg(3,2) - fg(4,2) + fg(2,2) * Om(2,3) - fg(1,2) * Om(1,3) * k
      ELSE
         F1 = fg(3,1) + fg(4,1) - fg(1,1) * Om(1,3) + fg(2,1) * Om(2,3) / k
         F2 = fg(3,2) + fg(4,2) + fg(2,2) * Om(2,3) - fg(1,2) * Om(1,3) * k
      ENDIF

      F2 = DABS(F2)

      END SUBROUTINE
! ======================================================================

! === Wacholder & Weihs (1972) =========================================
!     Wacholder, E., & Weihs, D. (1972). Slow motion of a fluid sphere in the vicinity of another sphere or a plane boundary. Chemical Engineering Science, 27(10), 1817-1828.
      SUBROUTINE WW72(si,al,F,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

      sum1o = 0d0
        rel = 1d3
          i = 1
      DO WHILE ( rel .GT. acu )
         n  = DBLE(i)
       sum1 = sum1o + lambdaWW72(n,al,si) ! (2.46)
        rel = DABS(sum1-sum1o)/DABS(sum1)
!       WRITE(*,*) "n_max, sum1 = ", i,sum1
      sum1o = sum1
          i = i + 1
      ENDDO

      lambdaWW = 4d0*(si+1d0)/(3d0*si+2d0)*DSINH(al)*sum1 ! Table 1
      F = 4d0/3d0*DSINH(al)*sum1    ! (2.45) normalized ONLY by —6πμaU
      F = F * (si+1d0)/(si+2d0/3d0) ! —4πμaV(1.5μᵣ+1)/(μᵣ+1)

      END SUBROUTINE

! === F U N C T I O N S :  W W (1972) ================================

      FUNCTION lambdaWW72(n,al,si)
      IMPLICIT DOUBLE PRECISION (A-Z)
      lambdaWW72 = 2d0*DSINH((2d0*n+1d0)*al)+(2d0*n+1d0)*DSINH(2d0*al)
      lambdaWW72 = lambdaWW72 - 4d0*DSINH((n+1d0/2d0)*al)**2
      lambdaWW72 = lambdaWW72 + ( (2d0*n+1d0)*DSINH(al) )**2
      lambdaWW72 = lambdaWW72 * si + (2d0*n+3d0)*DEXP(2d0*al)/2d0
      lambdaWW72 = lambdaWW72 + 2d0*DEXP(-(2d0*n+1d0)*al)-(2d0*n-1d0)*DEXP(-2d0*al)/2d0
      lambdaWW72 = lambdaWW72 * n*(n+1d0)/(2d0*n+3d0)/(2d0*n-1d0)
!          denom = 2d0*DSINH((2d0*n+1d0)*al)+(2d0*n-1d0)*DSINH(2d0*al)  ! typo?!
           denom = 2d0*DSINH((2d0*n+1d0)*al)+(2d0*n+1d0)*DSINH(2d0*al)  ! guess!
           denom = denom * si + 2d0*(DCOSH((2d0*n+1d0)*al)+DCOSH(2d0*al))
      lambdaWW72 = lambdaWW72 / denom
      END

! ======================================================================

! === Haber, Hetsroni, Solan (1973) - Explicit Method ==================
!     Haber, S., Hetsroni, G., & Solan, A. (1973). On the low Reynolds number motion of two droplets. International Journal of Multiphase Flow, 1(1), 57-71.
      SUBROUTINE HHS73EXP(opp,lama,lamb,al,be,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      LOGICAL opp

          j = 1 ! first droplet
1     sumao = 0d0
      sumbo = 0d0
        rel = 1d3
          i = 1
      DO WHILE ( rel .GT. acu )
         n  = DBLE(i)
         k1 = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)

       suma = sumao + ( delta0(n,k1,al,be)                           &
            + lama  *   delta1(n,k1,al,be)                           &
            + lamb  *   delta2(n,k1,al,be)                           &
            + lama*lamb*delta3(n,k1,al,be) ) / Delt(n,al,be,lama,lamb)

       sumb = sumbo + ( deltabar0(n,k1,al,be)                        &
            + lama  *   deltabar1(n,k1,al,be)                        &
            + lamb  *   deltabar2(n,k1,al,be)                        &
            + lama*lamb*deltabar3(n,k1,al,be))/Delt(n,al,be,lama,lamb)

         rel = MAX(DABS(suma-sumao)/DABS(suma),DABS(sumb-sumbo)/DABS(sumb))
!        WRITE(*,*) "n_max, suma, sumb = ", i,suma,sumb
         sumao = suma
         sumbo = sumb
             i = i + 1
      ENDDO

!     Lambda resistances are computed from [32] & [33]
!     Typo: "c" parameter must be "c**2" which cancels
!     out with the same factor in K1 function.
!     Similarly, factors sqrt(2) and (2n+1)**2 in delta
!     functions are ignored, not to be computed twice.

      IF ( j .EQ. 1 ) THEN
      Lambda11 = DSINH(al)/3d0 * suma
      Lambda12 = DSINH(al)/3d0 * sumb
           tmp =-al
            al =-be
            be = tmp
             j = 2 ! second droplet
            GOTO 1
      ENDIF
      
      Lambda22 = DSINH(al)/3d0 * suma
      Lambda21 = DSINH(al)/3d0 * sumb

      IF (opp) THEN
         F1 = Lambda11 - Lambda12
         F2 = Lambda22 - Lambda21
      ELSE
         F1 = Lambda11 + Lambda12
         F2 = Lambda22 + Lambda21
      ENDIF

!     Hadamard–Rybczyński normalization:
      F1 = F1*(lama+1d0)/(lama+2d0/3d0) ! —4πμa₁V₁(1.5μᵣ+1)/(μᵣ+1)
      F2 = F2*(lamb+1d0)/(lamb+2d0/3d0) ! —4πμa₂V₂(1.5μᵣ+1)/(μᵣ+1)

      END SUBROUTINE

! === F U N C T I O N S :  H H S (1973)  A P P E N D I X  C ============

      FUNCTION delta0(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      delta0 = (2d0*n+3d0)*DSINH((n+3d0/2d0)*(al-be))*DEXP(-(n-1d0/2d0)*(al+be))
      delta0 = delta0-(2d0*n-1d0)*DEXP(-(n+3d0/2d0)*(al+be))*DSINH((n-1d0/2d0)*(al-be))
      delta0 = delta0*4d0*k1
      END

      FUNCTION deltabar0(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      deltabar0 =-(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al-be))
      deltabar0 = deltabar0+(2d0*n-1d0)*DEXP(-(n+3d0/2d0)*(al-be))*DSINH((n-1d0/2d0)*(al-be))
      deltabar0 = deltabar0*4d0*k1
      END

      FUNCTION delta1(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      delta1 = 2d0*(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(2d0*n+1d0)*be)
      delta1 = delta1-(2d0*n+3d0)**2*DEXP(-(n-1d0/2d0)*(al+be))*DCOSH((n+3d0/2d0)*(al-be))
      delta1 = delta1-(2d0*n-1d0)**2*DEXP(-(n+3d0/2d0)*(al+be))*DCOSH((n-1d0/2d0)*(al-be))
      delta1 = delta1-(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*(al+be))*DSINH((n+3d0/2d0)*(al-be))
      delta1 = delta1-(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n+3d0/2d0)*(al+be))*DSINH((n-1d0/2d0)*(al-be)) ! note: typo
      delta1 =-delta1*k1
      END

      FUNCTION deltabar1(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      deltabar1 =-2d0*(2d0*n-1d0)*(2d0*n+3d0)*DCOSH(2d0*be)
      deltabar1 = deltabar1+(2d0*n+3d0)**2*DEXP(-(n-1d0/2d0)*(al-be))*DCOSH((n+3d0/2d0)*(al-be))
      deltabar1 = deltabar1+(2d0*n-1d0)**2*DEXP(-(n+3d0/2d0)*(al-be))*DCOSH((n-1d0/2d0)*(al-be))
      deltabar1 = deltabar1+(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al-be))
      deltabar1 = deltabar1+(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n+3d0/2d0)*(al-be))*DSINH((n-1d0/2d0)*(al-be))
      deltabar1 =-deltabar1*k1
      END

      FUNCTION delta2(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      delta2 = 2d0*(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(2d0*n+1d0)*al)
      delta2 = delta2-(2d0*n+3d0)**2*DEXP(-(n-1d0/2d0)*(al+be))*DCOSH((n+3d0/2d0)*(al-be))
      delta2 = delta2-(2d0*n-1d0)**2*DEXP(-(n+3d0/2d0)*(al+be))*DCOSH((n-1d0/2d0)*(al-be))
      delta2 = delta2+(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*(al+be))*DSINH((n+3d0/2d0)*(al-be))
      delta2 = delta2+(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n+3d0/2d0)*(al+be))*DSINH((n-1d0/2d0)*(al-be))
      delta2 =-delta2*k1
      END

      FUNCTION deltabar2(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      deltabar2 =-2d0*(2d0*n-1d0)*(2d0*n+3d0)*DCOSH(2d0*al)
      deltabar2 = deltabar2+(2d0*n+3d0)**2*DEXP(-(n-1d0/2d0)*(al-be))*DCOSH((n+3d0/2d0)*(al-be))
      deltabar2 = deltabar2+(2d0*n-1d0)**2*DEXP(-(n+3d0/2d0)*(al-be))*DCOSH((n-1d0/2d0)*(al-be))
      deltabar2 = deltabar2+(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al-be))
      deltabar2 = deltabar2+(2d0*n-1d0)*(2d0*n+3d0)*DEXP(-(n+3d0/2d0)*(al-be))*DSINH((n-1d0/2d0)*(al-be))
      deltabar2 =-deltabar2*k1
      END

      FUNCTION delta3(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      delta3 = (2d0*n-1d0)*(2d0*n+3d0)*(DEXP(-(2d0*n+1d0)*be)-DEXP(-(2d0*n+1d0)*al))
      delta3 = delta3-(2d0*n+1d0)*(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*(al+be))*DSINH((n+3d0/2d0)*(al-be))
      delta3 = delta3-(2d0*n+1d0)*(2d0*n-1d0)*DEXP(-(n+3d0/2d0)*(al+be))*DSINH((n-1d0/2d0)*(al-be))
      delta3 =-delta3*2d0*k1 ! typo: negative sign mentioned by Zinchenko
      END

      FUNCTION deltabar3(n,k1,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      deltabar3 =-2d0*(2d0*n-1d0)*(2d0*n+1d0)*(2d0*n+3d0)*DSINH(al-be)*DCOSH(al+be)
      deltabar3 = deltabar3+(2d0*n+1d0)**3*DSINH(2d0*(al-be))
      deltabar3 = deltabar3+2d0*(2d0*n+1d0)**2*DCOSH(2d0*(al-be))
      deltabar3 = deltabar3-8d0*DEXP(-(2d0*n+1d0)*(al-be))-2d0*(2d0*n-1d0)*(2d0*n+3d0)
      deltabar3 =-deltabar3*k1
      END

      FUNCTION Delt(n,al,be,lama,lamb)
      IMPLICIT DOUBLE PRECISION (A-Z)
      Delt = 4d0*DSINH((n-1d0/2d0)*(al-be))*DSINH((n+3d0/2d0)*(al-be))
      Delt = Delt+(lama+lamb)*(2d0*DSINH((2d0*n+1d0)*(al-be))-(2d0*n+1d0)*DSINH(2d0*(al-be)))
      Delt = Delt+(lama*lamb)*(4d0*DSINH((n+1d0/2d0)*(al-be))**2-((2d0*n+1d0)*DSINH(al-be))**2)
      END      
! ======================================================================

! === Haber, Hetsroni, Solan (1973) - Implicit Method ==================
!     Haber, S., Hetsroni, G., & Solan, A. (1973). On the low Reynolds number motion of two droplets. International Journal of Multiphase Flow, 1(1), 57-71.
      SUBROUTINE HHS73IMP(opp,lama,lamb,al,be,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, DIMENSION(8,9) :: T
      LOGICAL opp

      sum1o = 0d0
      sum2o = 0d0
        rel = 1d3
          i = 1
      DO WHILE ( rel .GT. acu )
         n = DBLE(i)
        k1 = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)
        k2 = n*(n+1d0)/4d0

!        Eqs. [B-4] to [B-11] of HHS73:

         T(1,1) = DCOSH((n-1d0/2d0)*al)
         T(1,2) = DSINH((n-1d0/2d0)*al)
         T(1,3) = DCOSH((n+3d0/2d0)*al)
         T(1,4) = DSINH((n+3d0/2d0)*al)
         T(1,5:8) = 0d0
         T(1,9) = k1*f(n,al)

         T(2,1) = DCOSH((n-1d0/2d0)*be)
         T(2,2) = DSINH((n-1d0/2d0)*be)
         T(2,3) = DCOSH((n+3d0/2d0)*be)
         T(2,4) = DSINH((n+3d0/2d0)*be)
         T(2,5:8) = 0d0
         IF (opp) THEN
         T(2,9) =-k1*f(n,-be) ! α and β moving in an opposing direction
         ELSE
         T(2,9) =+k1*f(n,-be) ! α and β moving in the same direction
         ENDIF

         T(3,1:4) = 0d0
         T(3,5) = DEXP(-(n-1d0/2d0)*al)
         T(3,6) = DEXP(-(n+3d0/2d0)*al)
         T(3,7:8) = 0d0
         T(3,9) = k1*f(n,al)

         T(4,1:6) = 0d0
         T(4,7) = DEXP((n-1d0/2d0)*be)
         T(4,8) = DEXP((n+3d0/2d0)*be)
         IF (opp) THEN
         T(4,9) =-k1*f(n,-be) ! α and β moving in an opposing direction
         ELSE
         T(4,9) =+k1*f(n,-be) ! α and β moving in the same direction
         ENDIF

         T(5,1) = (n-1d0/2d0)*DSINH((n-1d0/2d0)*al)
         T(5,2) = (n-1d0/2d0)*DCOSH((n-1d0/2d0)*al)
         T(5,3) = (n+3d0/2d0)*DSINH((n+3d0/2d0)*al)
         T(5,4) = (n+3d0/2d0)*DCOSH((n+3d0/2d0)*al)
         T(5,5) = (n-1d0/2d0)*DEXP(-(n-1d0/2d0)*al)
!        T(5,6) =-(n+3d0/2d0)*DEXP(+(n+3d0/2d0)*al) !HHS73 typo?!
         T(5,6) = (n+3d0/2d0)*DEXP(-(n+3d0/2d0)*al) !WW72
         T(5,7:9) = 0d0

         T(6,1) = (n-1d0/2d0)*DSINH((n-1d0/2d0)*be)
         T(6,2) = (n-1d0/2d0)*DCOSH((n-1d0/2d0)*be)
         T(6,3) = (n+3d0/2d0)*DSINH((n+3d0/2d0)*be)
         T(6,4) = (n+3d0/2d0)*DCOSH((n+3d0/2d0)*be)
         T(6,5:6) = 0d0
         T(6,7) =-(n-1d0/2d0)*DEXP(+(n-1d0/2d0)*be)
         T(6,8) =-(n+3d0/2d0)*DEXP(+(n+3d0/2d0)*be)
         T(6,9) = 0d0

         T(7,1) = (n-1d0/2d0)**2*DCOSH((n-1d0/2d0)*al)
         T(7,2) = (n-1d0/2d0)**2*DSINH((n-1d0/2d0)*al)
         T(7,3) = (n+3d0/2d0)**2*DCOSH((n+3d0/2d0)*al)
         T(7,4) = (n+3d0/2d0)**2*DSINH((n+3d0/2d0)*al)
         T(7,5) =-(n-1d0/2d0)**2*DEXP(-(n-1d0/2d0)*al)*lama
         T(7,6) =-(n+3d0/2d0)**2*DEXP(-(n+3d0/2d0)*al)*lama
         T(7,7:8) = 0d0
         T(7,9) = (1d0-lama)*k2*g(n,al)

         T(8,1) = (n-1d0/2d0)**2*DCOSH((n-1d0/2d0)*be)
         T(8,2) = (n-1d0/2d0)**2*DSINH((n-1d0/2d0)*be)
         T(8,3) = (n+3d0/2d0)**2*DCOSH((n+3d0/2d0)*be)
         T(8,4) = (n+3d0/2d0)**2*DSINH((n+3d0/2d0)*be)
         T(8,5:6) = 0d0
         T(8,7) =-(n-1d0/2d0)**2*DEXP(+(n-1d0/2d0)*be)*lamb
         T(8,8) =-(n+3d0/2d0)**2*DEXP(+(n+3d0/2d0)*be)*lamb
         IF (opp) THEN
         T(8,9) =-(1d0-lamb)*k2*g(n,-be) ! α and β moving in an opposing direction
         ELSE
         T(8,9) =+(1d0-lamb)*k2*g(n,-be) ! α and β moving in the same direction
         ENDIF

         CALL GAUSS(8,9,T)

         sum1 = sum1o + SUM(T(1:4,9))                ! Summation of Eq. [29] in HHS73
         sum2 = sum2o + T(1,9)-T(2,9)+T(3,9)-T(4,9)  ! Summation of Eq. [30] in HHS73
!        WRITE(*,*) "n_max, sum1, sum2 = ", i,sum1,sum2
         rel  = MAX(DABS((sum1-sum1o))/DABS(sum1),DABS((sum2-sum2o))/DABS(sum2))
         sum1o= sum1
         sum2o= sum2
            i = i + 1
      ENDDO

!     Rigid sphere Stokes drag normalization:      
      F1 = -DSINH(al)/3d0 *      sum1  ! normalized by —6πμa₁V₁
      F2 = -DSINH(be)/3d0 * DABS(sum2) ! normalized by —6πμa₂V₂

!     Hadamard–Rybczyński normalization:
      F1 = F1*(lama+1d0)/(lama+2d0/3d0) ! —4πμa₁V₁(1.5μᵣ+1)/(μᵣ+1)
      F2 = F2*(lamb+1d0)/(lamb+2d0/3d0) ! —4πμa₂V₂(1.5μᵣ+1)/(μᵣ+1)

      END SUBROUTINE

! === F U N C T I O N S :  H H S (1973)  A P P E N D I X  B ============
      FUNCTION f(n,x)
      IMPLICIT DOUBLE PRECISION (A-Z)
      f = (2d0*n-1d0)*DEXP(-(n+3d0/2d0)*x)-(2d0*n+3d0)*DEXP(-(n-1d0/2d0)*x)
      END

      FUNCTION g(n,x)
      IMPLICIT DOUBLE PRECISION (A-Z)
      g = (2d0*n+3d0)*DEXP(-(n+3d0/2d0)*x)-(2d0*n-1d0)*DEXP(-(n-1d0/2d0)*x)
      END
! ======================================================================

! === Reed & Morrison (1974) ===========================================
!     Reed, L. D., & Morrison Jr, F. A. (1974). Particle interactions in viscous flow at small values of Knudsen number. Journal of Aerosol Science, 5(2), 175-189.
      SUBROUTINE RM74(opp,xi1,xi2,a1,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T
      LOGICAL opp

      Cm = 1.0d0 ! Momentum accommodation (assumed)
      mfp= 0.1d0 ! Air mean free path
       c = a1 * DSINH(xi1)
      a2 = a1 * DSINH(xi1) /-DSINH(xi2)
      Clc= Cm * mfp / c
      C1 = Cm * mfp / a1
      C2 = Cm * mfp / a2

      F1o= 0d0
      F2o= 0d0
      iN = 50
1     ALLOCATE ( T(4*iN,14) )
       T = 0d0

      DO i = 1, iN
         n = DBLE(i)
         k = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)

         ! (19): ξ₁
         T(4*i-3,08) = DCOSH((n-05d-1)*xi1)
         T(4*i-3,09) = DSINH((n-05d-1)*xi1)
         T(4*i-3,10) = DCOSH((n+15d-1)*xi1)
         T(4*i-3,11) = DSINH((n+15d-1)*xi1)
         T(4*i-3,14)=-k*( (2d0*n+3d0)*DEXP(-(n-05d-1)*xi1) - (2d0*n-1d0)*DEXP(-(n+15d-1)*xi1) )
         ! the factors "Uc^2/sq(2)" are ignored due to (22) & (23)

         ! (19): ξ₂
         T(4*i-2,07) = DCOSH((n-05d-1)*xi2)
         T(4*i-2,08) = DSINH((n-05d-1)*xi2)
         T(4*i-2,09) = DCOSH((n+15d-1)*xi2)
         T(4*i-2,10) = DSINH((n+15d-1)*xi2)
         T(4*i-2,14) =-k*( (2d0*n+3d0)*DEXP((n-05d-1)*xi2) - (2d0*n-1d0)*DEXP((n+15d-1)*xi2) )
         IF (opp) T(4*i-2,14) = -T(4*i-2,14)
         
         ! (20): ξ₁
         IF (i.GT.1) THEN
         T(4*i-1,02) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi1)
         T(4*i-1,03) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi1)
         T(4*i-1,04) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         T(4*i-1,05) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         ENDIF
         T(4*i-1,06) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi1)
         T(4*i-1,07) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi1)
         T(4*i-1,08) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi1)
         T(4*i-1,09) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi1)
         IF (i.LT.iN) THEN
         T(4*i-1,10) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         T(4*i-1,11) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         T(4*i-1,12) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi1)
         T(4*i-1,13) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi1)
         ENDIF
         T(4*i-1,14) =-Clc*n*(n+1d0)*( DCOSH(xi1)* &
         ((2d0*n-1d0)*DEXP(-(n-05d-1)*xi1)-(2d0*n+3d0)*DEXP(-(n+15d-1)*xi1)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP(-(n-15d-1)*xi1)-(2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)-(2d0*n+5d0)*DEXP(-(n+25d-1)*xi1)))&
         +2d0*n*(n+1d0)      *            (DEXP(-(n-05d-1)*xi1)      -      DEXP(-(n+15d-1)*xi1) )

         ! (20): ξ₂
         IF (i.GT.1) THEN
         T(4*i,01) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi2)
         T(4*i,02) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi2)
         T(4*i,03) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         T(4*i,04) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         ENDIF
         T(4*i,05) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi2)
         T(4*i,06) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi2)
         T(4*i,07) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi2)
         T(4*i,08) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi2)
         IF (i.LT.iN) THEN
         T(4*i,09) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         T(4*i,10) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         T(4*i,11) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi2)
         T(4*i,12) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi2)
         ENDIF
         IF (opp) THEN
         T(4*i,14) =-Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         +2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ELSE
         T(4*i,14) =+Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         -2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ENDIF
      ENDDO

      CALL THOMAS(4*iN,7,5,T)

      F1 = SUM(T(:,14))
      F2 = 0d0
      DO i = 1, iN
         an = T(4*i-3,14)
         bn = T(4*i-2,14)
         cn = T(4*i-1,14)
         dn = T(4*i  ,14)
         F2 = F2 + an - bn + cn - dn
      ENDDO
      DEALLOCATE ( T )

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( rel .GT. acu ) THEN
          iN = INT(1.5 * FLOAT(iN)) ! 50% increase
         F1o = F1
         F2o = F2
!        WRITE(*,*) 'n_max, F1, F2, rel = ', iN,F1,F2,rel
         GOTO 1
      ENDIF

!     Normalized by (1): —6πμaᵢVᵢ(1+2Cl/aᵢ)/(1+3Cl/aᵢ)
      F1 = DSINH(xi1)/3d0 * (1d0+3d0*C1)/(1d0+2d0*C1) * DABS(F1)  ! (22)
      F2 =-DSINH(xi2)/3d0 * (1d0+3d0*C2)/(1d0+2d0*C2) * DABS(F2)  ! (23)

      END SUBROUTINE
! ======================================================================                 

! === Rother, Stark, Davis (2022) ======================================
!     Rother, M. A., Stark, J. K., & Davis, R. H. (2022). Gravitational collision efficiencies of small viscous drops at finite Stokes numbers and low Reynolds numbers. International Journal of Multiphase Flow, 146, 103876.
      SUBROUTINE RSD22(opp,mu1,mu2,eta1,eta2,a1,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,L-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T
      LOGICAL opp

      IF ( mu1 .EQ. 0d0 ) mu1 = 1d-6
      IF ( mu2 .EQ. 0d0 ) mu2 = 1d-6
      CS = 1.0d0 ! Momentum accommodation (assumed)
      mfp= 0.1d0 ! Air mean free path [μm]
       c = a1 * DSINH(eta1)
      a2 = a1 * DSINH(eta1) /-DSINH(eta2)
      Clc= CS * mfp / c
      C1 = CS * mfp / a1
      C2 = CS * mfp / a2
      C1 = (1d0+2d0/3d0/mu1+2d0*C1)/(1d0+1d0/mu1+3d0*C1)
      C2 = (1d0+2d0/3d0/mu2+2d0*C2)/(1d0+1d0/mu2+3d0*C2)

      F1o= 0d0
      F2o= 0d0
      iN = 50
1     ALLOCATE ( T(4*iN,14) )
       T = 0d0

      DO i = 1, iN
         n = DBLE(i)

         ! (A.11)
         T(4*i-3,08) = 1d0
         T(4*i-3,09) = DEXP((n-05d-1)*(eta2-eta1))
         T(4*i-3,10) = 1d0
         T(4*i-3,11) = DEXP((n+15d-1)*(eta2-eta1))
         T(4*i-3,14) = DEXP((n-05d-1)*-eta1)/(2d0*n-1d0)-DEXP((n+15d-1)*-eta1)/(2d0*n+3d0) ! typo

         ! (A.12)
         T(4*i-2,07) = DEXP((n-05d-1)*(eta2-eta1))
         T(4*i-2,08) = 1d0
         T(4*i-2,09) = DEXP((n+15d-1)*(eta2-eta1))
         T(4*i-2,10) = 1d0
         T(4*i-2,14) = DEXP((n-05d-1)*eta2)/(2d0*n-1d0)-DEXP((n+15d-1)*eta2)/(2d0*n+3d0) ! typo
         IF (opp) T(4*i-2,14) = -T(4*i-2,14)

         ! (A.13)
         fct1 =-Clc*(n-1d0)/(2d0*n-1d0)
         IF (i.GT.1) THEN
         T(4*i-1,02) = fct1 * (n-15d-1)**2
         T(4*i-1,03) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*(eta2-eta1))
         T(4*i-1,04) = fct1 * (n+05d-1)**2
         T(4*i-1,05) = fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         ENDIF

         fct2 = Clc*DCOSH(eta1) + (mu1*(2d0*n+1d0))**-1
         T(4*i-1,06) = fct2 * (n-05d-1)**2 + (n-05d-1)
         T(4*i-1,07) = fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*(eta2-eta1)) - (n-05d-1) * DEXP((n-05d-1)*(eta2-eta1))
         T(4*i-1,08) = fct2 * (n+15d-1)**2 + (n+15d-1)
         T(4*i-1,09) = fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*(eta2-eta1)) - (n+15d-1) * DEXP((n+15d-1)*(eta2-eta1))

         fct3 =-Clc*(n+2d0)/(2d0*n+3d0)
         IF (i.LT.iN) THEN
         T(4*i-1,10) = fct3 * (n+05d-1)**2
         T(4*i-1,11) = fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         T(4*i-1,12) = fct3 * (n+25d-1)**2
         T(4*i-1,13) = fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*(eta2-eta1))
         ENDIF

         T(4*i-1,14) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*-eta1)/(2d0*n-3d0) &
                     - fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*-eta1)/(2d0*n+1d0) &
                     + fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*-eta1)/(2d0*n-1d0) &
                     - fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*-eta1)/(2d0*n+3d0) &
                     + fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*-eta1)/(2d0*n+1d0) &
                     - fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*-eta1)/(2d0*n+5d0) &
                     +        (n+15d-1)    * DEXP((n+15d-1)*-eta1)/(2d0*n+3d0) &
                     -        (n-05d-1)    * DEXP((n-05d-1)*-eta1)/(2d0*n-1d0)

         ! (A.14)
         fct1 =-Clc*(n-1d0)/(2d0*n-1d0)
         IF (i.GT.1) THEN
!        T(4*i,01) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*eta2)
         T(4*i,01) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*(eta2-eta1)) ! typo
         T(4*i,02) = fct1 * (n-15d-1)**2
         T(4*i,03) = fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         T(4*i,04) = fct1 * (n+05d-1)**2
         ENDIF

         fct2 = Clc*DCOSH(eta2) + (mu2*(2d0*n+1d0))**-1
!        T(4*i,05) = fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*eta2) - (n-05d-1) * DEXP((n-05d-1)*(eta2-eta1))
         T(4*i,05) = fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*(eta2-eta1)) - (n-05d-1) * DEXP((n-05d-1)*(eta2-eta1)) ! typo
         T(4*i,06) = fct2 * (n-05d-1)**2 + (n-05d-1)
         T(4*i,07) = fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*(eta2-eta1)) - (n+15d-1) * DEXP((n+15d-1)*(eta2-eta1))
         T(4*i,08) = fct2 * (n+15d-1)**2 + (n+15d-1)

         fct3 =-Clc*(n+2d0)/(2d0*n+3d0)
         IF (i.LT.iN) THEN
         T(4*i,09) = fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         T(4*i,10) = fct3 * (n+05d-1)**2
         T(4*i,11) = fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*(eta2-eta1))
         T(4*i,12) = fct3 * (n+25d-1)**2
         ENDIF

         T(4*i,14) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*eta2)/(2d0*n-3d0) &
                   - fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*eta2)/(2d0*n+1d0) &
                   + fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*eta2)/(2d0*n-1d0) &
                   - fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*eta2)/(2d0*n+3d0) &
                   + fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*eta2)/(2d0*n+1d0) & ! typo
                   - fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*eta2)/(2d0*n+5d0) &
                   +        (n+15d-1)    * DEXP((n+15d-1)*eta2)/(2d0*n+3d0) &
                   -        (n-05d-1)    * DEXP((n-05d-1)*eta2)/(2d0*n-1d0)

         IF (opp) T(4*i,14) = -T(4*i,14)
      ENDDO

      CALL THOMAS(4*iN,7,5,T)

      F1 = 0d0
      F2 = 0d0
      DO i = 1, iN
         n = DBLE(i)
        En = T(4*i-3,14)
        Fn = T(4*i-2,14)
        Gn = T(4*i-1,14)
        Hn = T(4*i  ,14)
        F1 = F1 + n*(n+1d0) * ( En * DEXP(-(n-05d-1)*eta1) + Gn * DEXP(-(n+15d-1)*eta1) )
        F2 = F2 + n*(n+1d0) * ( Fn * DEXP( (n-05d-1)*eta2) + Hn * DEXP( (n+15d-1)*eta2) )
      ENDDO
      DEALLOCATE ( T )

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( rel .GT. acu ) THEN
          iN = INT(1.5 * FLOAT(iN)) ! 50% increase
         F1o = F1
         F2o = F2
!        WRITE(*,*) 'n_max, F1, F2, rel = ', iN,F1,F2,rel
         GOTO 1
      ENDIF

      F1 = 2d0/3d0 * DSINH(eta1) / C1 * DABS(F1)
      F2 =-2d0/3d0 * DSINH(eta2) / C2 * DABS(F2)

      END SUBROUTINE
! ======================================================================                 

! === Beshkov, Radoev, Ivanov (1978) ===================================
!     Beshkov, V. N., Radoev, B. P., & Ivanov, I. B. (1978). Slow motion of two droplets and a droplet towards a fluid or solid interface. International Journal of Multiphase Flow, 4(5-6), 563-570.
!     Kim, S., & Karrila, S. J. (2013). Microhydrodynamics: principles and selected applications. Courier Corporation.      
      SUBROUTINE BRI78(mur,al,F1)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

      sumo = 0d0
      DO i = 1, 100000
         n = DBLE(i)
        Kn = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0) ! two typos in Beshkov et al.

        N0 = 2d0*(2d0*n+1d0)*DSINH(2*al)+4d0*DCOSH(2*al)-4d0*DEXP(-(2d0*n+1d0)*al)
        N1 = (2d0*n+1d0)**2*DCOSH(2*al)+2d0*(2d0*n+1d0)*DSINH(2*al)-(2d0*n-1d0)*(2d0*n+3d0)+4d0*DEXP(-(2d0*n+1d0)*al) ! typo in (9.39) of "Microhydrodynamics" by Kim & Karrila
        D0 = 4d0*DSINH((n-1d0/2d0)*al)*DSINH((n+3d0/2d0)*al)
        D1 = 2d0*DSINH((2d0*n+1d0)*al)-(2d0*n+1d0)*DSINH(2*al) ! next typo in Beshkov et al.

      suma = sumo + Kn*( N0 + mur * N1 ) / ( D0 + mur * D1 )

         criterion = DABS((suma-sumo))/DABS(suma)
         IF ( criterion .LT. 1d-10 ) THEN
!           WRITE(*,*) "n_max, sum = ", i,suma
            GOTO 1
         ENDIF
         sumo = suma
      ENDDO

1     F1 = 2d0/3d0*DSINH(al) * suma  ! normalized by —6πμaV

      END SUBROUTINE
! ======================================================================

! === Jeffrey & Onishi (1984) - XA =====================================
!     Jeffrey, D. J., & Onishi, Y. (1984). Calculation of the resistance and mobility functions for two unequal rigid spheres in low-Reynolds-number flow. Journal of Fluid Mechanics, 139, 261-290.
!     Townsend, A. K. (2018). Generating, from scratch, the near-field asymptotic forms of scalar resistance functions for two unequal rigid spheres in low-Reynolds-number flow. arXiv preprint arXiv:1802.08226.      
      SUBROUTINE JO84XA(opp,ncl,s,alam,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fm1,fm2
      DOUBLE PRECISION Kn, lmbd0
      LOGICAL opp, ncl
      
!     IF ( acu .LT. 1d-6 ) WRITE(*,*) "Consider decreasing accuracy!"

      rlam = 1d0/alam
       F1o = 0d0
       F2o = 0d0
        n0 = 100
1     ALLOCATE ( fm1(0:n0), fm2(0:n0) )
      CALL coeffXA(n0,alam,rlam,fm1,fm2)

      g11=2.d0*alam*alam/(1.d0+alam)**3
      g21=0.2d0*alam*(1.d0+7.d0*alam+alam*alam)/(1.d0+alam)**3
      g31=1.d0/42.d0*(1.d0+18.d0*alam-29.d0*alam**2+18.d0*alam**3+alam**4)/(1.d0+alam)**3
      g12=2.d0*rlam*rlam/(1.d0+rlam)**3
      g22=0.2d0*rlam*(1.d0+7.d0*rlam+rlam*rlam)/(1.d0+rlam)**3
      g32=1.d0/42.d0*(1.d0+18.d0*rlam-29.d0*rlam**2+18.d0*rlam**3+rlam**4)/(1.d0+rlam)**3

      XA11=0d0
      XA12=0d0
      XA21=0d0
      XA22=0d0
      AX11=0d0
      AX12=0d0
      AX21=0d0
      AX22=0d0

      sc = 1.d0-4.d0/s/s
      sd = (s+2.d0)/(s-2.d0)
      xi = s-2.d0

      DO m = 1, n0
        m1 = m - 2
        IF(m.EQ.2) m1 = -2 

        dm = DBLE(m)
       dm1 = DBLE(m1)
       dm2 = DBLE(m+2)

       IF (MOD(m,2).EQ.0) THEN
          ccc = fm1(m)/(2.d0+2.d0*alam)**m-g11-g21*2.d0/dm+4.d0*g31/dm/dm1
         AX11 = AX11 + ccc
         XA11 = XA11 + ccc*(2/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m-g12-g22*2.d0/dm+4.d0*g32/dm/dm1
         AX22 = AX22 + ccc
         XA22 = XA22 + ccc*(2.d0/s)**m

       ELSE
          ccc = fm1(m)/(2.d0+2.d0*alam)**m-g11-g21*2.d0/dm+4.d0*g31/dm/dm1 ! corr: dm1->JO84 | dm2->T19
         AX12 = AX12 + ccc
         XA12 = XA12 + ccc*(2.d0/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m-g12-g22*2.d0/dm+4.d0*g32/dm/dm1 ! corr: dm1->JO84 | dm2->T19
         AX21 = AX21 + ccc
         XA21 = XA21 + ccc*(2.d0/s)**m
       ENDIF
      ENDDO
      DEALLOCATE ( fm1, fm2 )

      XA11 = XA11 + g11/sc - g21*DLOG(sc)-g31*sc*DLOG(sc)+1.d0-g11
      XA22 = XA22 + g12/sc - g22*DLOG(sc)-g32*sc*DLOG(sc)+1.d0-g12
      XA12 = XA12 + 2.d0/s*g11/sc + g21*DLOG(sd)+ g31*sc*DLOG(sd) + 4.d0*g31/s
      XA21 = XA21 + 2.d0/s*g12/sc + g22*DLOG(sd)+ g32*sc*DLOG(sd) + 4.d0*g32/s
      XA12 = XA12 *-2.d0/(1.d0+alam)
      XA21 = XA21 *-2.d0/(1.d0+rlam)

      IF ( ncl ) THEN
        lmbd0 = 1d-5 ! cm (0.1 micron)
           a1 = 10d0 / 1d4
           a2 = alam * a1

           Kn = 2d0 * lmbd0 / ( a1 + a2 )
!          Kn = lmbd0 * ( a1 + a2 ) / ( 2d0 * a1 * a2 )
           Kn = 1d-2

         dlt0 = ( s - 2d0 ) / Kn
         XA11c= XA11
         XA11 = XA11 + g11 * ( f_nc(dlt0)/Kn - 1d0/(s-2d0) )
         XA22 = XA22 + g12 * ( f_nc(dlt0)/Kn - 1d0/(s-2d0) )
         XA12 = XA12 + 2d0 * g11/(1d0+alam) * ( 1d0/(s-2d0) - f_nc(dlt0)/Kn )
         XA21 = XA21 + 2d0 * g12/(1d0+rlam) * ( 1d0/(s-2d0) - f_nc(dlt0)/Kn )
      ENDIF
      WRITE(*,*)"XAxx",XA11,XA12,XA21,XA22

      IF (opp) THEN
         F1 = XA11 - (1d0+alam)/2d0 * XA12
         F2 = XA22 - (1d0+rlam)/2d0 * XA21
      ELSE
         F1 = XA11 + (1d0+alam)/2d0 * XA12
         F2 = XA22 + (1d0+rlam)/2d0 * XA21
      ENDIF

!     rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
!     IF ( rel .GT. acu ) THEN
!        n0 = INT(1.1 * FLOAT(n0)) ! 10% increase
!        F1o= F1
!        F2o= F2
!        WRITE(*,*) "n_max, F1, F2 = ", n0, F1, F2, rel
!        GOTO 1
!     ENDIF

!     WRITE (2,*) s-2d0, XA11c, XA11

      AX11 = AX11 + 1d0 - g11/4d0
      AX22 = AX22 + 1d0 - g12/4d0
      AX12 = AX12 + g11/4d0 + 2d0*g21*DLOG(2d0) + 2d0*g31    ! typo in T19
      AX21 = AX21 + g12/4d0 + 2d0*g22*DLOG(2d0) + 2d0*g32    ! typo in T19
      AX12 =-AX12*2.d0/(1.d0+alam)
      AX21 =-AX21*2.d0/(1.d0+rlam)
!     WRITE (*,*) AX11,AX12,AX21,AX22

!     Nearly touching spheres (3.17) & (3.18)
      XA11=g11/xi+g21*DLOG(1d0/xi)+g31*xi*DLOG(1d0/xi)+AX11
      XA22=g12/xi+g22*DLOG(1d0/xi)+g32*xi*DLOG(1d0/xi)+AX22
      XA12=g11/xi+g21*DLOG(1d0/xi)+g31*xi*DLOG(1d0/xi)-(1.d0+alam)/2d0*AX12
      XA21=g12/xi+g22*DLOG(1d0/xi)+g32*xi*DLOG(1d0/xi)-(1.d0+rlam)/2d0*AX21
!     XA12=-XA12                                     ! T19
      XA12=-XA12*2.d0/(1.d0+alam)                    ! JO84
!     XA21=-XA21                                     ! T19
      XA21=-XA21*2.d0/(1.d0+rlam)                    ! JO84
!     WRITE (*,*) XA11,XA12
!     WRITE (*,*) XA11-(1d0+alam)/2d0*XA12,XA11+(1d0+alam)/2d0*XA12 ! JO84
!     WRITE (*,*) XA11-XA12,XA11+XA12                ! T19

      END SUBROUTINE
! ======================================================================

      SUBROUTINE coeffXA(n0,alam,rlam,fm1,fm2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     DOUBLE PRECISION P(0:2307277),V(0:2307277)
!     DOUBLE PRECISION FAC1(151,151)
!     DOUBLE PRECISION P(0:17160960),V(0:17160960)
!     DOUBLE PRECISION FAC1(301,301)
      DOUBLE PRECISION P(0:75840222),V(0:75840222)
      DOUBLE PRECISION FAC1(501,501)
      DOUBLE PRECISION XN,XS
      DOUBLE PRECISION FCTR1,FCTR2,FCTR3,FCTR4,FCTRQ,FCTRV,FACTOR
      DOUBLE PRECISION PSUM,VSUM,QSUM,FAC
      DOUBLE PRECISION fm1(0:n0),fm2(0:n0)
      INTEGER MMAX,MAXS,NM,N,IS,M,NQ,IQ,INDX2,JS,IP,KS,J

      MAXS=n0*2
      MMAX=MAXS+1
! First tabulate FAC1

      DO 10  N=1,(MMAX+1)/2
        XN=DBLE(N)
        FAC1(N,1)=XN+1D0
        DO 10 IS=2,(MMAX+1)/2
          XS=DBLE(IS)
 10       FAC1(N,IS)=FAC1(N,IS-1)*(XN+XS)/XS

! We first set the initial conditions for translation.
      P(INDX2(1,0,0))=1D0
      V(INDX2(1,0,0))=1D0

! We start at N=1, Q=0, P=M and proceed keeping M fixed.
! NM is the solution of the equation M=N+P+Q-1=NM+NM-1+0-1, except
! the last time through when only N=1 and 2 are calculated.
      DO 500 M=1,MMAX
       NM=M/2+1
       IF(M.EQ.MMAX)NM=2
       DO 400 N=1,NM

! The equation for NQ is M=N+P+Q-1=N+N-1+NQ-1
        NQ=M-2*N+2
        XN=DBLE(N)
        FCTR1=(XN+.5D0)/(XN+1D0)*XN
        FCTR2=-XN*(XN-.5D0)/(XN+1D0)
        FCTR3=-XN*(2D0*XN*XN-.5D0)/(XN+1D0)
        FCTRV=-2D0*XN/((XN+1D0)*(2D0*XN+3D0))

        DO 300 IQ=0,NQ
! Now that M, N, Q are specified, I know P (denoted IP) as well.
! We obtain JS from the equation JS-1=IQ-JS.
! We also need KS defined by KS-1=IQ-KS-1.
         IP=M-IQ-N+1
         JS=(IQ+1)/2
         KS=IQ/2

         PSUM=0D0
         VSUM=0D0
! If JS is less than 1, we skip the loop on IS and set the elements to 0
         IF (JS.LT.1) GO TO 250
         DO 200 IS=1,JS

          XS=DBLE(IS)
          FAC=FAC1(N,IS)

! We search the summations in Jeffrey & Onishi (1984), finding first
! all INDX2es that start IS,IQ-IS after which we test the final argument
! which can be IP-N+1, IP-N or IP-N-1. The range of N is chosen to keep IP-N+1
! non-negative, and the range of IS is chosen to keep the combination IS,IQ-IS
! legal. Thus the other possibilities require if tests and jumps.
          FACTOR=(2D0*XN*XS-XN-XS+2D0)/(XN+XS)

          PSUM=PSUM+FAC*FCTR1*FACTOR*P(INDX2(IS,IQ-IS,IP-N+1))/(2D0*XS-1D0)

          IF(IP.GT.N) THEN
            PSUM=PSUM+FAC*FCTR2*P(INDX2(IS,IQ-IS,IP-N-1))
            VSUM=VSUM+FAC*FCTRV*P(INDX2(IS,IQ-IS,IP-N-1))
          ENDIF
          IF(IS.LT.JS) THEN
            PSUM=PSUM+FAC*FCTR3*V(INDX2(IS,IQ-IS-2,IP-N+1))/(2D0*XS+1D0)
          ENDIF
  200    CONTINUE
  250    P(INDX2(N,IP,IQ))=PSUM
  300    V(INDX2(N,IP,IQ))=PSUM+VSUM
  400  CONTINUE
  500 CONTINUE

      fm1(0) = 1d0
      fm2(0) = 1d0
      fm1(1) = 3d0 * alam
      fm2(1) = 3d0 * rlam

      DO k = 2, n0
         fm1(k) = 0d0
         fm2(k) = 0d0
         DO  iq = 1, k
             fm1(k) = fm1(k) + P(INDX2(1,k-iq,iq)) * alam**iq
             fm2(k) = fm2(k) + P(INDX2(1,k-iq,iq)) * rlam**iq
         ENDDO
         fm1(k) = fm1(k) * 2d0**k
         fm2(k) = fm2(k) * 2d0**k
!        WRITE(*,*)"k:  f(k) = ",k, fm1(k), fm2(k)
      ENDDO

      END SUBROUTINE

! ==== F U N C T I O N :  N O N - C O N T I N U U M   L U B ============
!     Sundararajakumar, R. R., & Koch, D. L. (1996). Non-continuum lubrication flows between particles colliding in a gas. Journal of Fluid Mechanics, 313, 283-308.
!     Dhanasekaran, J., Roy, A., & Koch, D. L. (2021). Collision rate of bidisperse spheres settling in a compressional non-continuum gas flow. Journal of Fluid Mechanics, 910.
      FUNCTION f_nc ( dlt0 )
      IMPLICIT DOUBLE PRECISION (A-Z)

      pi = 4d0 * DATAN(1d0)

      IF ( dlt0 .LE. 0.26d0 ) THEN
      t0   = DLOG ( 1d0 / dlt0 ) + 0.4513d0
      f_nc = pi/6d0 * ( DLOG(t0) - 1d0/t0 - 1d0/t0**2 - 2d0/t0**3 ) + 2.587d0*dlt0**2 + 1.419d0*dlt0 + 0.3847d0

      ELSEIF ( dlt0 .GT. 0.26d0 .AND. dlt0 .LE. 5.08d0  ) THEN
      f_nc = 5.607d-4*dlt0**4 - 9.275d-3*dlt0**3 + 6.067d-2*dlt0**2 - 0.2082d0*dlt0 + 0.4654d0 + 0.05488d0 / dlt0

      ELSEIF ( dlt0 .GT. 5.08d0 .AND. dlt0 .LE. 10.55d0 ) THEN
      f_nc =-1.182d-4*dlt0**3 + 3.929d-3*dlt0**2 - 5.017d-2*dlt0 + 0.3102d0

      ELSEIF ( dlt0 .GT. 10.55d0 ) THEN
      f_nc = 0.0452d0 * ( ( 6.649d0 + dlt0 ) * DLOG ( 1d0+6.649d0/dlt0 ) - 6.649d0 )

      ENDIF
      END
! ================      
!     How, M. L. S., Koch, D. L., & Collins, L. R. (2021). Non-continuum tangential lubrication gas flow between two spheres. Journal of Fluid Mechanics, 920.
      FUNCTION f_nc2 ( dlt0 ) ! f_nc works better
      IMPLICIT DOUBLE PRECISION (A-Z)

      pi = 4d0 * DATAN(1d0)

      IF ( dlt0 .LE. 0.35d0 ) THEN
      t0    = DLOG ( 1d0 / dlt0 ) + 0.8d0
      f_nc2 = pi/6d0 * ( DLOG(t0) - 1d0/t0 - 1d0/t0**2 - 2d0/t0**3 ) &
            + 1.5475d0*dlt0**2 + 0.7896d0*dlt0 + 0.4094d0

      ELSEIF ( dlt0 .GT. 0.35d0 .AND. dlt0 .LE. 9d-1  ) THEN
      f_nc2 = 6d-5*dlt0**-3 - 172d-5*dlt0**-2 + 169d-4*dlt0**-1 &
            + 1769d-4*DLOG(dlt0**-1) + 3744d-4

      ELSEIF ( dlt0 .GT. 9d-1 .AND. dlt0 .LE. 4.4d0 ) THEN
      f_nc2 = 5d-4*dlt0**4 - 835d-5*dlt0**3 + 5605d-5*dlt0**2 - 2196d-4*dlt0 + 5661d-4 &
            - 55d-4*dlt0**-1 - 1965d-5*dlt0**-2 + 231d-4*dlt0**-3

      ELSEIF ( dlt0 .GT. 4.4d0 ) THEN
      k1    = 1.016d0
      f_nc2 = 1d0/(18d0*k1**2) * ( ( 6d0*k1 + dlt0 ) * DLOG ( 1d0+6d0*k1/dlt0 ) - 6d0*k1 )

      ENDIF

      END
! ======================================================================

! === Jeffrey & Onishi (1984) - YA =====================================
!     Jeffrey, D. J., & Onishi, Y. (1984). Calculation of the resistance and mobility functions for two unequal rigid spheres in low-Reynolds-number flow. Journal of Fluid Mechanics, 139, 261-290.
      SUBROUTINE JO84YA(opp,s,alam,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fm1,fm2
      LOGICAL opp

      IF ( acu .LT. 1d-6 ) WRITE(*,*) "Consider decreasing accuracy!"

      rlam = 1d0/alam
       F1o = 0d0
       F2o = 0d0
        n0 = 100
1     ALLOCATE ( fm1(0:n0), fm2(0:n0) )
      CALL coeffYA(n0,alam,rlam,fm1,fm2)

      g21=4.d0/15.d0*alam*(2.d0+alam+2.d0*alam*alam)/(1.d0+alam)**3
      g31=2.d0/375.d0*(16.d0-45.d0*alam+58.d0*alam**2-45.d0*alam**3+16.d0*alam**4)/(1.d0+alam)**3
      g22=4.d0/15.d0*rlam*(2.d0+rlam+2.d0*rlam*rlam)/(1.d0+rlam)**3
      g32=2.d0/375.d0*(16.d0-45.d0*rlam+58.d0*rlam**2-45.d0*rlam**3+16.d0*rlam**4)/(1.d0+rlam)**3

      sc = 1.d0-4.d0/s/s
      sd = (s+2.d0)/(s-2.d0)

      YA11=0d0
      YA12=0d0
      YA21=0d0
      YA22=0d0
      AY11=0d0
      AY12=0d0
      AY21=0d0
      AY22=0d0

      DO m=1,n0  
       
         m1 = m - 2 
         IF(m.EQ.2) m1 = -2 

         dm1=DBLE(m1)
         dm=DBLE(m)

         IF (MOD(m,2).EQ.0) THEN
          ccc = fm1(m)/(2.d0+2.d0*alam)**m - g21*2.d0/dm + 4.d0*g31/dm/dm1
         AY11 = AY11 + ccc
         YA11 = YA11 + ccc*(2.d0/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m - g22*2.d0/dm + 4.d0*g32/dm/dm1
         AY22 = AY22 + ccc
         YA22 = YA22 + ccc*(2.d0/s)**m

         ELSE
          ccc = fm1(m)/(2.d0+2.d0*alam)**m - g21*2.d0/dm + 4.d0*g31/dm/dm1
         AY12 = AY12 + ccc
         YA12 = YA12 + ccc*(2.d0/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m - g22*2.d0/dm + 4.d0*g32/dm/dm1
         AY21 = AY21 + ccc
         YA21 = YA21 + ccc*(2.d0/s)**m
         ENDIF
      ENDDO
      DEALLOCATE ( fm1, fm2 )

      YA11 = YA11 - g21*DLOG(sc) - g31*sc*DLOG(sc) + 1.d0
      YA22 = YA22 - g22*DLOG(sc) - g32*sc*DLOG(sc) + 1.d0 
      YA12 = YA12 + g21*DLOG(sd) + g31*sc*DLOG(sd) + 4.d0*g31/s
      YA21 = YA21 + g22*DLOG(sd) + g32*sc*DLOG(sd) + 4.d0*g32/s
      YA12 =-YA12 * 2d0/(1d0+alam)
      YA21 =-YA21 * 2d0/(1d0+rlam)

      IF (opp) THEN
         F1 = YA11 - (1d0+alam)/2d0 * YA12
         F2 = YA22 - (1d0+rlam)/2d0 * YA21
      ELSE
         F1 = YA11 + (1d0+alam)/2d0 * YA12
         F2 = YA22 + (1d0+rlam)/2d0 * YA21
      ENDIF

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
!     IF ( rel .GT. acu ) THEN
         n0 = INT(1.1 * FLOAT(n0)) ! 10% increase
         F1o= F1
         F2o= F2
!        WRITE(*,*) "n_max, F1, F2 = ", n0, F1, F2, rel
!        GOTO 1
!     ENDIF

      AY11 = AY11 + 1d0
      AY22 = AY22 + 1d0
      AY12 = AY12 + 2d0*g21*DLOG(2d0) + 2d0*g31
      AY21 = AY21 + 2d0*g22*DLOG(2d0) + 2d0*g32
      AY12 =-AY12*2.d0/(1.d0+alam)
      AY21 =-AY21*2.d0/(1.d0+rlam)
!     WRITE (*,*) AY11,AY12,AY21,AY22

      END SUBROUTINE
          
! ======================================================================
      SUBROUTINE coeffYA(n0,alam,rlam,fm1,fm2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION fm1(0:n0),fm2(0:n0)
!     DOUBLE PRECISION P(0:2307277)
!     DOUBLE PRECISION V(0:2307277)
!     DOUBLE PRECISION Q(0:2307277),FAC1(151,151)
      DOUBLE PRECISION P(0:16409500)
      DOUBLE PRECISION V(0:16409500)
      DOUBLE PRECISION Q(0:16409500),FAC1(301,301)
      DOUBLE PRECISION XN,XS
      DOUBLE PRECISION FCTR1,FCTR2,FCTR3,FCTR4,FCTRQ,FCTRV,FACTOR
      DOUBLE PRECISION PSUM,VSUM,QSUM,FAC
      INTEGER MMAX,MAXS,NM,N,IS,M,NQ,IQ,INDX2,JS,IP,KS,J
            
      MAXS=n0*2
      MMAX=MAXS+1
! First tabulate FAC1
      DO 10 N=1,(MMAX+1)/2
        XN=DBLE(N)
        FAC1(N,1)=1D0
        DO 10 IS=2,(MMAX+1)/2
          XS=DBLE(IS)
10    FAC1(N,IS)=FAC1(N,IS-1)*(XN+XS)/(XS-1D0)
! We first set the initial conditions for translation.
      P(INDX2(1,0,0))=1D0
      V(INDX2(1,0,0))=1D0
      Q(INDX2(1,0,0))=0D0
! We start at N=1, Q=0, P=M and proceed keeping M fixed.
!  NM is the solution of the equation M=N+P+Q-1=NM+NM-1+0-1, except
!  the last time through when only N=1 and 2 are calculated.
      DO 500 M=1,MMAX
       NM=M/2+1
       IF(M.EQ.MMAX)NM=2
       DO 400 N=1,NM
! The equation for NQ is M=N+P+Q-1=N+N-1+NQ-1
        NQ=M-2*N+2
        XN=DBLE(N)
        FCTR1=(XN+.5D0)/(XN+1D0)
        FCTR2=XN*(XN-.5D0)/(XN+1D0)
        FCTR3=XN*(2D0*XN*XN-.5D0)/(XN+1D0)
        FCTR4=(8D0*XN*XN-2D0)/(3D0*XN+3D0)
        FCTRV=2D0*XN/((XN+1D0)*(2D0*XN+3D0))
        FCTRQ=3D0/(2D0*XN*(XN+1D0))
        DO 300 IQ=0,NQ
! Now that M, N, Q are specified, I know P (denoted IP) as well.
!  We obtain JS from the equation JS-1=IQ-JS.
!  We also need KS defined by KS-1=IQ-KS-1.
         IP=M-IQ-N+1
         JS=(IQ+1)/2
         KS=IQ/2
         PSUM=0D0
         VSUM=0D0
         QSUM=0D0
! If JS is less than 1, we skip the loop on IS and set the elements to 0
         IF (JS.LT.1) GO TO 250
         DO 200 IS=1,JS
          XS=DBLE(IS)
          FAC=FAC1(N,IS)
! We search the summations in Jeffrey & Onishi (1984), finding first
!  all INDX2es that start IS,IQ-IS after which we test the final argument
!  which can be IP-N+1, IP-N or IP-N-1. The range of N is chosen to keep IP-N+1
!  non-negative, and the range of IS is chosen to keep the combination IS,IQ-IS
!  legal. Thus the other possibilities require if tests and jumps.

          FACTOR=2D0*(XN*XS+1D0)**2/(XN+XS)
          PSUM=PSUM+FAC*FCTR1*(4D0+XN*XS-FACTOR)*P(INDX2(IS,IQ-IS,IP-N+1))/(XS*(2D0*XS-1D0))
          IF(IP.GE.N) THEN
            QSUM=QSUM-FAC*FCTRQ*P(INDX2(IS,IQ-IS,IP-N))/XS
          ENDIF
          IF(IP.GT.N) THEN
            PSUM=PSUM+FAC*FCTR2*P(INDX2(IS,IQ-IS,IP-N-1))
! A typographical error in equn (4.9) of Jeffrey & Onishi (1984)
!    has been corrected in the following statement.
            VSUM=VSUM+FAC*FCTRV*P(INDX2(IS,IQ-IS,IP-N-1))
          ENDIF
          IF(IS.LE.KS) THEN
            PSUM=PSUM-FAC*FCTR4*Q(INDX2(IS,IQ-IS-1,IP-N+1))
            IF(IP.GE.N) THEN
              QSUM=QSUM+FAC*XS*Q(INDX2(IS,IQ-IS-1,IP-N))/(XN+1D0)
             ENDIF
          ENDIF
          IF(IS.LT.JS) THEN
            PSUM=PSUM+FAC*FCTR3*V(INDX2(IS,IQ-IS-2,IP-N+1))/(2D0*XS+1D0)
          ENDIF
  200    CONTINUE
  250    P(INDX2(N,IP,IQ))=PSUM
         V(INDX2(N,IP,IQ))=PSUM+VSUM
  300    Q(INDX2(N,IP,IQ))=QSUM
  400  CONTINUE
  500 CONTINUE

      fm1(0) = 1d0
      fm2(0) = 1d0
      fm1(1) = 3d0/2d0 * alam
      fm2(1) = 3d0/2d0 * rlam

      DO k = 2, n0
         fm1(k) = 0d0
         fm2(k) = 0d0
         DO iq = 0, k
            fm1(k) = fm1(k) + P(INDX2(1,k-iq,iq)) * alam**iq
            fm2(k) = fm2(k) + P(INDX2(1,k-iq,iq)) * rlam**iq
         ENDDO
         fm1(k) = fm1(k) * 2d0**k
         fm2(k) = fm2(k) * 2d0**k
!        WRITE(*,*)"k:  f(k) = ",k, fm1(k)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === Jeffrey & Onishi (1984) - YB =====================================
!     Jeffrey, D. J., & Onishi, Y. (1984). Calculation of the resistance and mobility functions for two unequal rigid spheres in low-Reynolds-number flow. Journal of Fluid Mechanics, 139, 261-290.
!     Townsend, A. K. (2018). Generating, from scratch, the near-field asymptotic forms of scalar resistance functions for two unequal rigid spheres in low-Reynolds-number flow. arXiv preprint arXiv:1802.08226.
      SUBROUTINE JO84YB(opp,s,alam,F1,F2,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fm1,fm2
      LOGICAL opp

      IF ( acu .LT. 1d-6 ) WRITE(*,*) "Consider decreasing accuracy!"

      rlam = 1d0/alam
       F1o = 0d0
       F2o = 0d0
        n0 = 100
1     ALLOCATE ( fm1(0:n0), fm2(0:n0) )
      CALL coeffYB(n0,alam,rlam,fm1,fm2)

      g21=-1.d0/5.d0*alam*(4.d0+alam)/(1.d0+alam)**2
      g31=-1.d0/250.d0*(32.d0-33.d0*alam+83.d0*alam**2+43.d0*alam**3)/(1.d0+alam)**2
      g22=-1.d0/5.d0*rlam*(4.d0+rlam)/(1.d0+rlam)**2
      g32=-1.d0/250.d0*(32.d0-33.d0*rlam+83.d0*rlam**2+43.d0*rlam**3)/(1.d0+rlam)**2

      sc = 1.d0-4.d0/s/s
      sd = (s+2.d0)/(s-2.d0)

      YB11=0d0
      YB12=0d0
      YB21=0d0
      YB22=0d0
      BY11=0d0
      BY12=0d0
      BY21=0d0
      BY22=0d0

      DO m=1,n0  
      
         m1 = m - 2 
         IF(m.EQ.2) m1 = -2 

         dm1=DBLE(m1)
         dm=DBLE(m)

         IF(MOD(m,2).EQ.1) THEN
          ccc = fm1(m)/(2.d0+2.d0*alam)**m - g21*2.d0/dm + 4.d0*g31/dm/dm1
         BY11 = BY11 + ccc
         YB11 = YB11 + ccc*(2.d0/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m - g22*2.d0/dm + 4.d0*g32/dm/dm1
         BY22 = BY22 + ccc
         YB22 = YB22 + ccc*(2.d0/s)**m

         ELSE
          ccc = fm1(m)/(2.d0+2.d0*alam)**m - g21*2.d0/dm + 4.d0*g31/dm/dm1
         BY12 = BY12 + ccc
         YB12 = YB12 + ccc*(2.d0/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m - g22*2.d0/dm + 4.d0*g32/dm/dm1
         BY21 = BY21 + ccc
         YB21 = YB21 + ccc*(2.d0/s)**m
         ENDIF
      ENDDO
      DEALLOCATE ( fm1, fm2 )

      YB11 = YB11 + g21*DLOG(sd) + g31*sc*DLOG(sd) + 4.d0*g31/s
      YB22 = YB22 + g22*DLOG(sd) + g32*sc*DLOG(sd) + 4.d0*g32/s
      YB12 = YB12 - g21*DLOG(sc) - g31*sc*DLOG(sc)
      YB21 = YB21 - g22*DLOG(sc) - g32*sc*DLOG(sc)
      YB12 = YB12 *-4d0/(1d0+alam)**2
      YB21 = YB21 *-4d0/(1d0+rlam)**2

      IF (opp) THEN
         F1 = 2d0/3d0 * (-YB11 - (1d0+alam)**2/4d0 * YB21 )
         F2 = 2d0/3d0 * ( YB22 + (1d0+rlam)**2/4d0 * YB12 )
         T1 = 1d0/2d0 * (-YB11 + (1d0+alam)**2/4d0 * YB12 )
         T2 = 1d0/2d0 * ( YB22 - (1d0+rlam)**2/4d0 * YB21 )
      ELSE
         F1 = 2d0/3d0 * (-YB11 + (1d0+alam)**2/4d0 * YB21 )
         F2 = 2d0/3d0 * ( YB22 - (1d0+rlam)**2/4d0 * YB12 )
         T1 = 1d0/2d0 * (-YB11 - (1d0+alam)**2/4d0 * YB12 )
         T2 = 1d0/2d0 * ( YB22 + (1d0+rlam)**2/4d0 * YB21 )
      ENDIF

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( rel .GT. acu ) THEN
         n0 = INT(1.1 * FLOAT(n0)) ! 10% increase
         F1o= F1
         F2o= F2
!        WRITE(*,*) "n_max, F1, F2, T1, T2 = ", n0, F1, F2, T1, T2, rel
         GOTO 1
      ENDIF

!     JO84: (seems wrong)      
!     BY11 = BY11
!     BY22 = BY22
!     BY11 = BY11 + 2d0*g21*DLOG(2d0) + 2d0*g31
!     BY22 = BY22 + 2d0*g22*DLOG(2d0) + 2d0*g32

!     T19: (seems wrong)      
!     BY11 = BY11 + 2d0*g21*DLOG(2d0) - 2d0*g31
!     BY22 = BY22 + 2d0*g22*DLOG(2d0) - 2d0*g32
!     BY11 = BY11 - g31
!     BY22 = BY22 - g32

!     correct?      
      BY11 = BY11 + 2d0*g21*DLOG(2d0) + 2d0*g31
      BY22 = BY22 + 2d0*g22*DLOG(2d0) + 2d0*g32
      BY12 =-BY12*4.d0/(1.d0+alam)**2
      BY21 =+BY21*4.d0/(1.d0+rlam)**2
      BY22 =-BY22
!     WRITE (*,*) BY11,BY12,BY21,BY22,"BY not fixed"

      END SUBROUTINE

! ======================================================================
      SUBROUTINE coeffYB(n0,alam,rlam,fm1,fm2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION fm1(0:n0),fm2(0:n0)
!     DOUBLE PRECISION P(0:2307277)
!     DOUBLE PRECISION V(0:2307277)
!     DOUBLE PRECISION Q(0:2307277),FAC1(151,151)
      DOUBLE PRECISION P(0:16409500)
      DOUBLE PRECISION V(0:16409500)
      DOUBLE PRECISION Q(0:16409500),FAC1(301,301)
      DOUBLE PRECISION XN,XS
      DOUBLE PRECISION FCTR1,FCTR2,FCTR3,FCTR4,FCTRQ,FCTRV,FACTOR
      DOUBLE PRECISION PSUM,VSUM,QSUM,FAC
      INTEGER MMAX,MAXS,NM,N,IS,M,NQ,IQ,INDX2,JS,IP,KS,J
            
      MAXS=n0*2
      MMAX=MAXS+1
! First tabulate FAC1
      DO 10 N=1,(MMAX+1)/2
        XN=DBLE(N)
        FAC1(N,1)=1D0
        DO 10 IS=2,(MMAX+1)/2
          XS=DBLE(IS)
10    FAC1(N,IS)=FAC1(N,IS-1)*(XN+XS)/(XS-1D0)
! We first set the initial conditions for translation.
      P(INDX2(1,0,0))=1D0
      V(INDX2(1,0,0))=1D0
      Q(INDX2(1,0,0))=0D0
! We start at N=1, Q=0, P=M and proceed keeping M fixed.
!  NM is the solution of the equation M=N+P+Q-1=NM+NM-1+0-1, except
!  the last time through when only N=1 and 2 are calculated.
      DO 500 M=1,MMAX
       NM=M/2+1
       IF(M.EQ.MMAX)NM=2
       DO 400 N=1,NM
! The equation for NQ is M=N+P+Q-1=N+N-1+NQ-1
        NQ=M-2*N+2
        XN=DBLE(N)
        FCTR1=(XN+.5D0)/(XN+1D0)
        FCTR2=XN*(XN-.5D0)/(XN+1D0)
        FCTR3=XN*(2D0*XN*XN-.5D0)/(XN+1D0)
        FCTR4=(8D0*XN*XN-2D0)/(3D0*XN+3D0)
        FCTRV=2D0*XN/((XN+1D0)*(2D0*XN+3D0))
        FCTRQ=3D0/(2D0*XN*(XN+1D0))
        DO 300 IQ=0,NQ
! Now that M, N, Q are specified, I know P (denoted IP) as well.
!  We obtain JS from the equation JS-1=IQ-JS.
!  We also need KS defined by KS-1=IQ-KS-1.
         IP=M-IQ-N+1
         JS=(IQ+1)/2
         KS=IQ/2
         PSUM=0D0
         VSUM=0D0
         QSUM=0D0
! If JS is less than 1, we skip the loop on IS and set the elements to 0
         IF (JS.LT.1) GO TO 250
         DO 200 IS=1,JS
          XS=DBLE(IS)
          FAC=FAC1(N,IS)
! We search the summations in Jeffrey & Onishi (1984), finding first
!  all INDX2es that start IS,IQ-IS after which we test the final argument
!  which can be IP-N+1, IP-N or IP-N-1. The range of N is chosen to keep IP-N+1
!  non-negative, and the range of IS is chosen to keep the combination IS,IQ-IS
!  legal. Thus the other possibilities require if tests and jumps.

          FACTOR=2D0*(XN*XS+1D0)**2/(XN+XS)
          PSUM=PSUM+FAC*FCTR1*(4D0+XN*XS-FACTOR)*P(INDX2(IS,IQ-IS,IP-N+1))/(XS*(2D0*XS-1D0))
          IF(IP.GE.N) THEN
            QSUM=QSUM-FAC*FCTRQ*P(INDX2(IS,IQ-IS,IP-N))/XS
          ENDIF
          IF(IP.GT.N) THEN
            PSUM=PSUM+FAC*FCTR2*P(INDX2(IS,IQ-IS,IP-N-1))
! A typographical error in equn (4.9) of Jeffrey & Onishi (1984)
!    has been corrected in the following statement.
            VSUM=VSUM+FAC*FCTRV*P(INDX2(IS,IQ-IS,IP-N-1))
          ENDIF
          IF(IS.LE.KS) THEN
            PSUM=PSUM-FAC*FCTR4*Q(INDX2(IS,IQ-IS-1,IP-N+1))
            IF(IP.GE.N) THEN
              QSUM=QSUM+FAC*XS*Q(INDX2(IS,IQ-IS-1,IP-N))/(XN+1D0)
             ENDIF
          ENDIF
          IF(IS.LT.JS) THEN
            PSUM=PSUM+FAC*FCTR3*V(INDX2(IS,IQ-IS-2,IP-N+1))/(2D0*XS+1D0)
          ENDIF
  200    CONTINUE
  250    P(INDX2(N,IP,IQ))=PSUM
         V(INDX2(N,IP,IQ))=PSUM+VSUM
  300    Q(INDX2(N,IP,IQ))=QSUM
  400  CONTINUE
  500 CONTINUE

      fm1(0) = 0d0
      fm2(0) = 0d0
      fm1(1) = 0d0
      fm2(1) = 0d0

      DO k = 2, n0
         fm1(k) = 0d0
         fm2(k) = 0d0
         DO iq = 0, k
            fm1(k) = fm1(k) + Q(INDX2(1,k-iq,iq)) * alam**iq
            fm2(k) = fm2(k) + Q(INDX2(1,k-iq,iq)) * rlam**iq
         ENDDO
         fm1(k) = fm1(k) * 2d0**k
         fm2(k) = fm2(k) * 2d0**k
      ENDDO
      DO k = 0, n0
         fm1(k) = 2d0 * fm1(k)
         fm2(k) = 2d0 * fm2(k)
!        WRITE(*,*)"k:  f(k) = ",k, fm1(k)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === Jeffrey & Onishi (1984) - YC =====================================
!     Jeffrey, D. J., & Onishi, Y. (1984). Calculation of the resistance and mobility functions for two unequal rigid spheres in low-Reynolds-number flow. Journal of Fluid Mechanics, 139, 261-290.
!     Townsend, A. K. (2018). Generating, from scratch, the near-field asymptotic forms of scalar resistance functions for two unequal rigid spheres in low-Reynolds-number flow. arXiv preprint arXiv:1802.08226.
!     http://ryuon.sourceforge.net/twobody/errata.html
      SUBROUTINE JO84YC(opp,s,alam,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fm1,fm2
      LOGICAL opp

      IF ( acu .LT. 1d-6 ) WRITE(*,*) "Consider decreasing accuracy!"

      rlam = 1d0/alam
       T1o = 0d0
       T2o = 0d0
        n0 = 150
1     ALLOCATE ( fm1(0:n0), fm2(0:n0) )
      CALL coeffYC(n0,alam,rlam,fm1,fm2)

      g21=2.d0/5.0*alam/(1.d0+alam)
      g31=1.d0/125.d0*(8.d0+6.d0*alam+33.d0*alam*alam)/(1.d0+alam)
!     g41=4.d0/5.d0*alam*alam/(1.d0+alam)**4                                  ! JO84 looks wrong
!     g51=4.d0/125.d0*alam*(43.d0-24.d0*alam+43.d0*alam*alam)/(1.d0+alam)**4  ! JO84
!     g51=2.d0/125.d0*alam*(43.d0-24.d0*alam+43.d0*alam*alam)/(1.d0+alam)**4  ! T19 & Ladd
      g51=1.d0/125.d0*alam*(43.d0-24.d0*alam+43.d0*alam*alam)/(1.d0+alam)**4  ! Ichiki
      g41=1.d0/10.d0*alam*alam/(1.d0+alam)                                    ! KK looks correct
!     g51=1.d0/250.d0*alam*(43.d0-24.d0*alam+43.d0*alam*alam)/(1.d0+alam)     ! KK

      g22=2.d0/5.0*rlam/(1.d0+rlam)
      g32=1.d0/125.d0*(8.d0+6.d0*rlam+33.d0*rlam*rlam)/(1.d0+rlam)
!     g42=4.d0/5.d0*rlam*rlam/(1.d0+rlam)**4                                  ! JO84 looks wrong
!     g52=4.d0/125.d0*rlam*(43.d0-24.d0*rlam+43.d0*rlam*rlam)/(1.d0+rlam)**4  ! JO84
!     g52=2.d0/125.d0*rlam*(43.d0-24.d0*rlam+43.d0*rlam*rlam)/(1.d0+rlam)**4  ! T19 & Ladd
      g52=1.d0/125.d0*rlam*(43.d0-24.d0*rlam+43.d0*rlam*rlam)/(1.d0+rlam)**4  ! Ichiki
      g42=1.d0/10.d0*rlam*rlam/(1.d0+rlam)                                    ! KK looks correct
!     g52=1.d0/250.d0*rlam*(43.d0-24.d0*rlam+43.d0*rlam*rlam)/(1.d0+rlam)     ! KK

      sc = 1.d0-4.d0/s/s
      sd = (s+2.d0)/(s-2.d0)

      YC11=0d0
      YC12=0d0
      YC21=0d0
      YC22=0d0
      CY11=0d0
      CY12=0d0
      CY21=0d0
      CY22=0d0

      DO m=1,n0
         m1 = m - 2
         IF(m.EQ.2) m1 = -2 

         dm1=DBLE(m1)
         dm=DBLE(m)

         IF(MOD(m,2).EQ.0) THEN
          ccc = fm1(m)/(2.d0+2.d0*alam)**m - g21*2.d0/dm + 4.d0*g31/dm/dm1
         CY11 = CY11 + ccc
         YC11 = YC11 + ccc*(2/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m - g22*2.d0/dm + 4.d0*g32/dm/dm1
         CY22 = CY22 + ccc
         YC22 = YC22 + ccc*(2.d0/s)**m

         ELSE
          ccc = fm1(m)/(2.d0+2.d0*alam)**m - g41*2.d0/dm + 4.d0*g51/dm/dm1
         CY12 = CY12 + ccc
         YC12 = YC12 + ccc*(2.d0/s)**m
          ccc = fm2(m)/(2.d0+2.d0*rlam)**m - g42*2.d0/dm + 4.d0*g52/dm/dm1
         CY21 = CY21 + ccc
         YC21 = YC21 + ccc*(2.d0/s)**m
         ENDIF
      ENDDO
      DEALLOCATE ( fm1, fm2 )

      YC11 = YC11 - g21*DLOG(sc) - g31*sc*DLOG(sc) + 1.d0
      YC22 = YC22 - g22*DLOG(sc) - g32*sc*DLOG(sc) + 1.d0
      YC12 = YC12 + g41*DLOG(sd) + g51*sc*DLOG(sd) + 4.d0*g51/s
      YC21 = YC21 + g42*DLOG(sd) + g52*sc*DLOG(sd) + 4.d0*g52/s
      YC12 = YC12 * 8d0/(1d0+alam)**3
      YC21 = YC21 * 8d0/(1d0+rlam)**3

!     WRITE(*,*)"YC11,YC12,YC21,YC22",YC11,YC12,YC21,YC22
!     WRITE(*,*)"YC11, (1d0+alam)**3/8d0 * YC12",YC11, (1d0+alam)**3/8d0 * YC12
!     WRITE(*,*)"YC22, (1d0+rlam)**3/8d0 * YC21",YC22, (1d0+rlam)**3/8d0 * YC21

      IF (opp) THEN
         T1 = YC11 - (1d0+alam)**3/8d0 * YC12
         T2 = YC22 - (1d0+rlam)**3/8d0 * YC21
      ELSE
         T1 = YC11 + (1d0+alam)**3/8d0 * YC12
         T2 = YC22 + (1d0+rlam)**3/8d0 * YC21
      ENDIF

      rel = MAX(DABS(T1-T1o)/DABS(T1),DABS(T2-T2o)/DABS(T2))
      IF ( rel .GT. acu ) THEN
         n0 = INT(1.1 * FLOAT(n0)) ! 10% increase
         T1o= T1
         T2o= T2
!        WRITE(*,*) "n_max, T1, T2 = ", n0, T1, T2, rel
         GOTO 1
      ENDIF

      CY11 = CY11 + 1d0
      CY22 = CY22 + 1d0
      CY12 = CY12 + 2d0*g41*DLOG(2d0) - 2d0*g51 ! Doesn't correspond with Table 6 JO84
      CY21 = CY21 + 2d0*g42*DLOG(2d0) - 2d0*g52
!     WRITE (*,*) CY11,CY12,CY21,CY22

      END SUBROUTINE

! ======================================================================
      SUBROUTINE coeffYC(n0,alam,rlam,fm1,fm2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION fm1(0:n0),fm2(0:n0)
!     DOUBLE PRECISION P(0:2307277)
!     DOUBLE PRECISION V(0:2307277)
!     DOUBLE PRECISION Q(0:2307277),FAC1(151,151)
      DOUBLE PRECISION P(0:16409500)
      DOUBLE PRECISION V(0:16409500)
      DOUBLE PRECISION Q(0:16409500),FAC1(301,301)
      DOUBLE PRECISION XN,XS
      DOUBLE PRECISION FCTR1,FCTR2,FCTR3,FCTR4,FCTRQ,FCTRV,FACTOR
      DOUBLE PRECISION PSUM,VSUM,QSUM,FAC
      INTEGER MMAX,MAXS,NM,N,IS,M,NQ,IQ,INDX2,JS,IP,KS,J

      MAXS=n0*2
      MMAX=MAXS+1
! First tabulate FAC1
      DO 10 N=1,(MMAX+1)/2
        XN=DBLE(N)
        FAC1(N,1)=1D0
        DO 10 IS=2,(MMAX+1)/2
          XS=DBLE(IS)
10    FAC1(N,IS)=FAC1(N,IS-1)*(XN+XS)/(XS-1D0)
! We first set the initial conditions for translation.
      P(INDX2(1,0,0))=0D0
      V(INDX2(1,0,0))=0D0
      Q(INDX2(1,0,0))=1D0
! We start at N=1, Q=0, P=M and proceed keeping M fixed.
!  NM is the solution of the equation M=N+P+Q-1=NM+NM-1+0-1, except
!  the last time through when only N=1 and 2 are calculated.
      DO 500 M=1,MMAX
       NM=M/2+1
       IF(M.EQ.MMAX)NM=2
       DO 400 N=1,NM
! The equation for NQ is M=N+P+Q-1=N+N-1+NQ-1
        NQ=M-2*N+2
        XN=DBLE(N)
        FCTR1=(XN+.5D0)/(XN+1D0)
        FCTR2=XN*(XN-.5D0)/(XN+1D0)
        FCTR3=XN*(2D0*XN*XN-.5D0)/(XN+1D0)
        FCTR4=(8D0*XN*XN-2D0)/(3D0*XN+3D0)
        FCTRV=2D0*XN/((XN+1D0)*(2D0*XN+3D0))
        FCTRQ=3D0/(2D0*XN*(XN+1D0))
        DO 300 IQ=0,NQ
! Now that M, N, Q are specified, I know P (denoted IP) as well.
!  We obtain JS from the equation JS-1=IQ-JS.
!  We also need KS defined by KS-1=IQ-KS-1.
         IP=M-IQ-N+1
         JS=(IQ+1)/2
         KS=IQ/2
         PSUM=0D0
         VSUM=0D0
         QSUM=0D0
! If JS is less than 1, we skip the loop on IS and set the elements to 0
         IF (JS.LT.1) GO TO 250
         DO 200 IS=1,JS
          XS=DBLE(IS)
          FAC=FAC1(N,IS)
! We search the summations in Jeffrey & Onishi (1984), finding first
!  all INDX2es that start IS,IQ-IS after which we test the final argument
!  which can be IP-N+1, IP-N or IP-N-1. The range of N is chosen to keep IP-N+1
!  non-negative, and the range of IS is chosen to keep the combination IS,IQ-IS
!  legal. Thus the other possibilities require if tests and jumps.

          FACTOR=2D0*(XN*XS+1D0)**2/(XN+XS)
          PSUM=PSUM+FAC*FCTR1*(4D0+XN*XS-FACTOR)*P(INDX2(IS,IQ-IS,IP-N+1))/(XS*(2D0*XS-1D0))
          IF(IP.GE.N) THEN
            QSUM=QSUM-FAC*FCTRQ*P(INDX2(IS,IQ-IS,IP-N))/XS
          ENDIF
          IF(IP.GT.N) THEN
            PSUM=PSUM+FAC*FCTR2*P(INDX2(IS,IQ-IS,IP-N-1))
! A typographical error in equn (4.9) of Jeffrey & Onishi (1984)
!    has been corrected in the following statement.
            VSUM=VSUM+FAC*FCTRV*P(INDX2(IS,IQ-IS,IP-N-1))
          ENDIF
          IF(IS.LE.KS) THEN
            PSUM=PSUM-FAC*FCTR4*Q(INDX2(IS,IQ-IS-1,IP-N+1))
            IF(IP.GE.N) THEN
              QSUM=QSUM+FAC*XS*Q(INDX2(IS,IQ-IS-1,IP-N))/(XN+1D0)
             ENDIF
          ENDIF
          IF(IS.LT.JS) THEN
            PSUM=PSUM+FAC*FCTR3*V(INDX2(IS,IQ-IS-2,IP-N+1))/(2D0*XS+1D0)
          ENDIF
  200    CONTINUE
  250    P(INDX2(N,IP,IQ))=PSUM
         V(INDX2(N,IP,IQ))=PSUM+VSUM
  300    Q(INDX2(N,IP,IQ))=QSUM
  400  CONTINUE
  500 CONTINUE

      fm1(0) = 1d0
      fm2(0) = 1d0
      fm1(1) = 0d0
      fm2(1) = 0d0

      DO k = 2, n0
         fm1(k) = 0d0
         fm2(k) = 0d0
         DO iq = 1, k
            fm1(k) = fm1(k) + Q(INDX2(1,k-iq,iq)) * alam**(iq+MOD(k,2))
            fm2(k) = fm2(k) + Q(INDX2(1,k-iq,iq)) * rlam**(iq+MOD(k,2))
         ENDDO
         fm1(k) = fm1(k) * 2d0**k
         fm2(k) = fm2(k) * 2d0**k
!        WRITE(*,*)"k:  f(k) = ",k, fm1(k)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === Jeffrey & Onishi (1984) - XC =====================================
!     Jeffrey, D. J., & Onishi, Y. (1984). Calculation of the resistance and mobility functions for two unequal rigid spheres in low-Reynolds-number flow. Journal of Fluid Mechanics, 139, 261-290.
!     http://ryuon.sourceforge.net/twobody/errata.html
!     Rosa, B., Wang, L. P., Maxey, M. R., & Grabowski, W. W. (2011). An accurate and efficient method for treating aerodynamic interactions of cloud droplets. Journal of Computational Physics, 230(22), 8109-8133.
      SUBROUTINE JO84XC(opp,s,alam,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fm1,fm2
!     REAL*16, ALLOCATABLE, DIMENSION(:) :: fm1,fm2
      LOGICAL opp

      IF ( acu .LT. 1d-6 ) WRITE(*,*) "Consider decreasing accuracy!"

      rlam = 1d0/alam
       T1o = 0d0
       T2o = 0d0
        n0 = 50
1     ALLOCATE ( fm1(0:2*n0+1), fm2(0:2*n0+1) )
      CALL coeffXC(n0,alam,rlam,fm1,fm2)

      sc = 1.d0-4.d0/s/s
      sd = (s+2.d0)/(s-2.d0)

      XC11=0d0
      XC12=0d0
      XC21=0d0
      XC22=0d0

      DO k=1,n0

         dk=DBLE(k)

         ccc = fm1(2*k)/(1d0+alam)**(2*k) - 2d0**(2*k+1)/dk/(2d0*dk-1d0) * alam**2/4d0/(1d0+alam)
        XC11 = XC11 + ccc * s**(-2*k)
         ccc = fm2(2*k)/(1d0+rlam)**(2*k) - 2d0**(2*k+1)/dk/(2d0*dk-1d0) * rlam**2/4d0/(1d0+rlam)
        XC22 = XC22 + ccc * s**(-2*k)

!        ccc = fm1(2*k+1)/(1d0+alam)**(2*k+1) - 2d0**(2*k+2)/dk/(2d0*dk+1d0) * alam**2    /(1d0+alam)        ! JO84
         ccc = fm1(2*k+1)/(1d0+alam)**(2*k+1) - 2d0**(2*k+2)/dk/(2d0*dk+1d0) * alam**2/4d0/(1d0+alam)        ! Ichiki: missing 4
!        ccc = fm1(2*k+1)/(1d0+alam)**(2*k+1) - 2d0**(2*k-1)/dk*(3d0*dk-2d0)/(2d0*dk-1d0)*alam**2/(1d0+alam) ! RWMG11: (44)
        XC12 = XC12 + ccc * s**(-2*k-1)
         ccc = fm2(2*k+1)/(1d0+rlam)**(2*k+1) - 2d0**(2*k+2)/dk/(2d0*dk+1d0) * rlam**2/4d0/(1d0+rlam)
        XC21 = XC21 + ccc * s**(-2*k-1)
      ENDDO
      DEALLOCATE ( fm1, fm2 )

      XC11 = XC11 + alam**2/(2d0+2d0*alam)*DLOG(sc) + alam**2/(1d0+alam)/s*DLOG(sd) + 1d0
      XC22 = XC22 + rlam**2/(2d0+2d0*rlam)*DLOG(sc) + rlam**2/(1d0+rlam)/s*DLOG(sd) + 1d0
!     XC12 =-8d0/(1d0+alam)**3*XC12+4d0*alam**2/(1d0+alam)**4*DLOG(sd)+8d0*alam**2/(1d0+alam)**4/s*DLOG(sc)     ! JO84
      XC12 =-8d0/(1d0+alam)**3*XC12+4d0*alam**2/(1d0+alam)**4*DLOG(sd)+8d0*alam**2/(1d0+alam)**4/s*DLOG(sc) &   ! Ichiki: missing last factor
                                                                     -16d0*alam**2/(1d0+alam)**4/s              ! Ichiki: missing last factor
!     XC12 =-8d0/(1d0+alam)**3*XC12+4d0*(alam/s)**2/(1d0+alam)**4*DLOG(sd)+8d0*alam**2/(1d0+alam)**4/s*DLOG(sc) ! RWMG11: (44)
      XC21 =-8d0/(1d0+rlam)**3*XC21+4d0*rlam**2/(1d0+rlam)**4*DLOG(sd)+8d0*rlam**2/(1d0+rlam)**4/s*DLOG(sc) &
                                                                     -16d0*rlam**2/(1d0+rlam)**4/s

      IF (opp) THEN
         T1 = XC11 - (1d0+alam)**3/8d0 * XC12
         T2 = XC22 - (1d0+rlam)**3/8d0 * XC21
      ELSE
         T1 = XC11 + (1d0+alam)**3/8d0 * XC12
         T2 = XC22 + (1d0+rlam)**3/8d0 * XC21
      ENDIF

      rel = MAX(DABS(T1-T1o)/DABS(T1),DABS(T2-T2o)/DABS(T2))
!     IF ( rel .GT. acu ) THEN
         n0 = INT(1.1 * FLOAT(n0)) ! 10% increase
         T1o= T1
         T2o= T2
!        WRITE(*,*) "n_max, T1, T2 = ", n0, T1, T2, rel
!        GOTO 1
!     ENDIF

      END SUBROUTINE

! ======================================================================
      SUBROUTINE coeffXC(n0,alam,rlam,fm1,fm2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     DOUBLE PRECISION Q(0:2284175),FAC0(151,151),QSUM,XN,XS
      DOUBLE PRECISION Q(0:16409500),FAC0(301,301),QSUM,XN,XS
      DOUBLE PRECISION fm1(0:2*n0+1),fm2(0:2*n0+1)
!     REAL*16          fm1(0:2*n0+1),fm2(0:2*n0+1)
      INTEGER MAXS,MMAX,NM,IS,M,NQ,INDX2,KS,N,IP,IQ

      MAXS=n0*2
      MMAX=MAXS+1
! First tabulate FAC0
      DO 10 N=1,(MMAX+1)/2
        XN=DBLE(N)
        FAC0(N,1)=1D0+XN
        DO 10 IS=2,(MMAX+1)/2
          XS=DBLE(IS)
   10 FAC0(N,IS)=FAC0(N,IS-1)*(XN+XS)/XS
! We set the initial conditions for rotation.
      Q(INDX2(1,0,0))=1D0
      DO 500 M=1,MMAX
        NM=M/2+1
        IF(M.EQ.MMAX)NM=1
        DO 400 N=1,NM
          NQ=M-2*N+2
          XN=DBLE(N)
          DO 300 IQ=0,NQ
            IP=M-IQ-N+1
            KS=IQ/2
            QSUM=0D0
            IF(KS.LT.1.OR.IP.LT.N) GOTO 300
            DO 200 IS=1,KS
              XS=DBLE(IS)
  200       QSUM=QSUM+FAC0(N,IS)*XS*Q(INDX2(IS,IQ-IS-1,IP-N))
  300     Q(INDX2(N,IP,IQ))=QSUM/(XN+1D0)
  400   CONTINUE
  500 CONTINUE

      fm1(0) = 1d0
      fm2(0) = 1d0
      fm1(1) = 0d0
      fm2(1) = 0d0

      DO k = 2, 2*n0+1
         fm1(k) = 0d0
         fm2(k) = 0d0
         DO iq = 1, k
            fm1(k) = fm1(k) + Q(INDX2(1,k-iq,iq)) * alam**(iq+MOD(k,2))
            fm2(k) = fm2(k) + Q(INDX2(1,k-iq,iq)) * rlam**(iq+MOD(k,2))
         ENDDO
         fm1(k) = fm1(k) * 2d0**k
         fm2(k) = fm2(k) * 2d0**k
!        WRITE(*,*)"k:  f(k) = ",k, fm1(k), fm2(k)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === F U N C T I O N S :  J O (1984) ==================================
      INTEGER FUNCTION INDX2(N,IP,IQ)
!  This function maps 3-dimensional arrays onto 1-dimensional arrays.
!  For each array element P(n,p,q) we define M=n+p+q-1, which is the primary
!  variable. For each value of M, we find that the array elements are
!  non-zero only for 1=< n =<NM = M/2+1. Further, for each fixed M and n,
!  the array elements that are non-zero are:
!                     P(n,m-n+1,0) -----> P(n,n-1,m-2n+2).
!  INDX2 uses these facts to set up a triangular scheme for each M to convert
!   it to a linear array. The array is filled as follows:
!   n,p,q=   1,M,0    1,M-1,1  .......................  1,1,M-1  1,0,M
!            *****    2,M-1,0  2,M-2,1 ...............  2,1,M-2
!            *****    *******     ....................
!            *****    *******     n,M-n+1,0  ...  n,n-1,M-2n+2
!  Note that the integer arithmetic relies on correct divisibility
! ---------------------------- Arguments -------------------------
      INTEGER N,IP,IQ
! -------------------------- Local Variables ---------------------
      INTEGER M,M2,LM
      M=IP+IQ+N-1
      M2=M/2
      LM=( M2*(M2+1) )/2
      LM=LM*(M2+1)+ ( LM*(M2+2) )/3
! LM is calculated as above to avoid overflow.
! If MMAX=200, THEN LM = 681750
      LM=LM+(M-2*M2)*(M2+1)*(M2+1)
      INDX2=LM+(N-1)*(M-N+3)+IQ
      RETURN
      END
! ======================================================================

! === Wang, Ayala, Grabowski (2005) - ISM, X-type setting: FIG3 (a) ====
!     Wang, L. P., Ayala, O., & Grabowski, W. W. (2005). Improved formulations of the superposition method. Journal of the atmospheric sciences, 62(4), 1255-1266.
!     Consider two spheres following (same orientation+) or
!     approaching/retreating from (opposing orientation-) each other
!     along their line of centers with the same velocity Vx.
!     The system of equations (19) & (20) for u = (u1x, u1y, u2x, u2y)
!     in x-y plane will be:
!     u1x = (L2+B2) * (V2x-u2x)  : V2x = +/- V1x
!     u1y = -B2 * u2y            : V2y = 0
!     u2x = (L1+B1) * (V1x-u1x)
!     u2y = -B1 * u1y            : V1y = 0
!     where L and B are the first and second factors in Eq. (1),
!     respectively. In what follows velocities are normalized by V1x.

      SUBROUTINE WAG05ISMX(opp,VR,s,alam,F1,F2)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION, DIMENSION(4,5) :: T
      LOGICAL opp

      A11 = (2d0+3d0*VR)/2d0/(1d0+VR)
      B11 = VR/4d0/(1d0+VR)

      AR = 2d0/(1d0+alam)/s           ! a1/r
!     L1 = 3d0/4d0*AR*(1d0-AR**2)     ! rigid
!     B1 = 1d0/4d0*AR*(3d0+1d0*AR**2) ! rigid
      L1 = A11/2d0*AR-3d0*B11*AR**3   ! fluid (includes rigid)
      B1 = A11/2d0*AR+B11*AR**3       ! fluid (includes rigid)

      AR = 2d0/(1d0+alam)/s*alam      ! a2/r
!     L2 = 3d0/4d0*AR*(1d0-AR**2)     ! rigid
!     B2 = 1d0/4d0*AR*(3d0+1d0*AR**2) ! rigid
      L2 = A11/2d0*AR-3d0*B11*AR**3   ! fluid (includes rigid)
      B2 = A11/2d0*AR+B11*AR**3       ! fluid (includes rigid)

      T(:,:) = 0d0
      T(1,1) = 1d0
      T(1,3) = L2+B2
      IF (opp) THEN
         T(1,5) =-(L2+B2)
      ELSE
         T(1,5) =+(L2+B2)
      ENDIF
      T(2,2) = 1d0
      T(2,4) = B2
      T(3,1) = L1+B1
      T(3,3) = 1d0
      T(3,5) = L1+B1
      T(4,2) = B1
      T(4,4) = 1d0

      CALL GAUSS(4,5,T)

!     Hadamard–Rybczyński normalization: —4πμaᵢVᵢ(1.5μᵣ+1)/(μᵣ+1)
      F1 = 1d0 - T(1,5)
      IF (opp) THEN
         F2 = 1d0 + T(3,5)
      ELSE
         F2 = 1d0 - T(3,5)
      ENDIF

      END SUBROUTINE
! ======================================================================

! === Wang, Ayala, Grabowski (2005) - ISM, Y-type setting: FIG3 (b) ====
!     Wang, L. P., Ayala, O., & Grabowski, W. W. (2005). Improved formulations of the superposition method. Journal of the atmospheric sciences, 62(4), 1255-1266.
!     Consider two spheres moving side by side perpENDicular to their 
!     line of centers with the same velocity Vy. The the second sphere
!     has the same or an opposing orientation V2y = +/- V1y to the first
!     one. The system of equations (19) & (20) for u = (u1x, u1y, u2x, u2y)
!     in x-y plane will be:
!     u1x = -(L2+B2) * u2x  : V2x = 0
!     u1y = B2 * (V2y-u2y)  : V2y = +/- V1y
!     u2x = -(L1+B1) * u1x  : V1x = 0
!     u2y = B1 * (V1y-u1y)
!     where L and B are the first and second factors in Eq. (1),
!     respectively. In what follows velocities are normalized by V1y.

      SUBROUTINE WAG05ISMY(opp,VR,s,alam,F1,F2)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      DOUBLE PRECISION, DIMENSION(4,5) :: T
      LOGICAL opp

      A11 = (2d0+3d0*VR)/2d0/(1d0+VR)
      B11 = VR/4d0/(1d0+VR)

      AR = 2d0/(1d0+alam)/s           ! a1/r
!     L1 = 3d0/4d0*AR*(1d0-AR**2)     ! rigid
!     B1 = 1d0/4d0*AR*(3d0+1d0*AR**2) ! rigid
      L1 = A11/2d0*AR-3d0*B11*AR**3   ! fluid (includes rigid)
      B1 = A11/2d0*AR+B11*AR**3       ! fluid (includes rigid)

      AR = 2d0/(1d0+alam)/s*alam      ! a2/r
!     L2 = 3d0/4d0*AR*(1d0-AR**2)     ! rigid
!     B2 = 1d0/4d0*AR*(3d0+1d0*AR**2) ! rigid
      L2 = A11/2d0*AR-3d0*B11*AR**3   ! fluid (includes rigid)
      B2 = A11/2d0*AR+B11*AR**3       ! fluid (includes rigid)

      T(:,:) = 0d0
      T(1,1) = 1d0
      T(1,3) = L2+B2
      T(2,2) = 1d0
      T(2,4) = B2
      IF (opp) THEN
         T(2,5) =-B2
      ELSE
         T(2,5) =+B2
      ENDIF
      T(3,1) = L1+B1
      T(3,3) = 1d0
      T(4,2) = B1
      T(4,4) = 1d0
      T(4,5) = B1

      CALL GAUSS(4,5,T)

!     Hadamard–Rybczyński normalization: —4πμaᵢVᵢ(1.5μᵣ+1)/(μᵣ+1)
      F1 = 1d0 - T(2,5)
      IF (opp) THEN
         F2 = 1d0 + T(4,5)
      ELSE
         F2 = 1d0 - T(4,5)
      ENDIF

      END SUBROUTINE
! ======================================================================

! === Rotational Superposition =========================================

      SUBROUTINE ROT(opp,s,alam,F1,F2)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      LOGICAL opp

      AR = 2d0/(1d0+alam)/s        ! a1/r
      A1 = AR**3
      AR = 2d0/(1d0+alam)/s*alam   ! a2/r
      A2 = AR**3

      IF (opp) THEN
         F1 =-A2
      ELSE
         F1 =+A2
      ENDIF

      F2 =-A1

      END SUBROUTINE
! ======================================================================

! === Goddard, Mills-Williams, Sun (2020) - Normal Interaction =========
!     Goddard, B. D., Mills-Williams, R. D., & Sun, J. (2020). The singular hydrodynamic interactions between two spheres in Stokes flow. Physics of Fluids, 32(6), 062001.
!     This is the corrected form of the solution provided by Maude (1961): typos found and mentioned earlier
      SUBROUTINE GMS20a(eta1,eta2,F1,F2)
      IMPLICIT DOUBLE PRECISION (a-h,k-z)
      sum1o = 0d0
      sum2o = 0d0
      DO i  = 1, 10000
         n  = DBLE(i)
         sum1 = sum1o + ( an(n,eta1,eta2) + bn(n,eta1,eta2)   &
                      +   cn(n,eta1,eta2) + dn(n,eta1,eta2)   &
                        )                 / Deltan(n,eta1,eta2)

         sum2 = sum2o + ( an(n,eta1,eta2) - bn(n,eta1,eta2)   &
                      +   cn(n,eta1,eta2) - dn(n,eta1,eta2)   &
                        )                 / Deltan(n,eta1,eta2)

         criterion = MAX(DABS((sum1-sum1o)/sum1), DABS((sum2-sum2o)/sum2))
         IF ( criterion .lt. 1d-10 ) THEN
!           WRITE(*,*) "n_max, sum1, sum2 = ", i,sum1,sum2
            GOTO 3
         ENDIF
         sum1o = sum1
         sum2o = sum2
      ENDDO

3     F1 = -1d0/3d0 * DSINH(eta1) * sum1 ! (38)
      F2 = -1d0/3d0 * DSINH(eta2) * sum2 ! (39)

      END SUBROUTINE

! === Goddard, Mills-Williams, Sun (2020) - Tangential Interaction =====
!     Goddard, B. D., Mills-Williams, R. D., & Sun, J. (2020). The singular hydrodynamic interactions between two spheres in Stokes flow. Physics of Fluids, 32(6), 062001.
      SUBROUTINE GMS20b(al,F1)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: At
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

      F1o= 0d0
      iN = 10
1     ALLOCATE ( T(iN,4), At(-1:iN+1) )
      At = 0d0
      DO i = 1, iN
         n = DBLE(i)

         T(i,1) = (n-1d0)*(gma(n-1d0,al)-1d0)-(n-1d0)*(2d0*n-3d0)/(2d0*n-1d0)*(gma(n,al)-1d0)
         T(i,2) = (2d0*n+1d0)-5d0*gma(n,al)-n*(2d0*n-1d0)/(2d0*n+1d0)*(gma(n-1d0,al)+1d0)   &
                + (n+1d0)*(2d0*n+3d0)/(2d0*n+1d0)*(gma(n+1d0,al)-1d0)
         T(i,3) = (n+2d0)*(2d0*n+5d0)/(2d0*n+3d0)*(gma(n,al)+1d0)-(n+2d0)*(gma(n+1d0,al)+1d0)
         T(i,4) =-DSQRT(2d0)*DEXP(-(n+5d-1)*al)*( DEXP(al)/DSINH((n-5d-1)*al)   &
                - 2d0/DSINH((n+5d-1)*al)+DEXP(al)/DSINH((n+15d-1)*al) )
      ENDDO

      CALL TDMA(iN,T)

      At(1:iN) = T(1:iN,4)

      sumF = 0d0
      DO i = 0, iN
         n = DBLE(i)
         C_n = 2d0*(n-1d0)/(2d0*n-1d0)*(gma(n,al)-1d0)*At(i-1)     &
             - 2d0*gma(n,al)*At(i) + 2d0*(n+2d0)/(2d0*n+3d0)       &  ! typo
             * (gma(n,al)+1d0)*At(i+1)

         E_n = DSQRT(8d0)*DEXP(-(n+5d-1)*al)/DSINH((n+5d-1)*al)    &
             - n*(n-1d0)/(2d0*n-1d0)*(gma(n,al)-1d0)      *At(i-1) &
             + (n+1d0)*(n+2d0)/(2d0*n+3d0)*(gma(n,al)+1d0)*At(i+1)
        sumF = sumF  + E_n + n*(n+1d0) * C_n
      ENDDO
      DEALLOCATE ( T, At )

      F1 = DSQRT(2d0)/6d0 * DSINH(al) * sumF ! (3.57)
      criterion = DABS(F1-F1o)/DABS(F1)
      IF ( criterion .GT. 1d-10 ) THEN
         iN = INT(1.2 * FLOAT(iN)) ! 20% increase
         F1o= F1
!        WRITE(*,*) "n_max, F1 = ", iN,F1
         GOTO 1
      ENDIF

      END SUBROUTINE

! === F U N C T I O N S :  G M S (2020) ================================

      FUNCTION an(n,eta1,eta2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      an = (2d0*n+1d0)*(n-1d0/2d0)*(DCOSH(2*eta1)-DCOSH(2*eta2))
      an = an - 2d0*(2d0*n-1d0)*DSINH((n+1d0/2d0)*(eta1-eta2))*DSINH((n+1d0/2d0)*(eta1+eta2))
      an = an + 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(eta1-eta2))*DSINH((n-1d0/2d0)*(eta1+eta2))
      an = an * (2d0*n+3d0)
      END

      FUNCTION bn(n,eta1,eta2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      bn = (2d0*n+1d0)*(n-1d0/2d0)*(DSINH(2*eta2)-DSINH(2*eta1))
      bn = bn - 2d0*(2d0*n-1d0)*DSINH((n+1d0/2d0)*(eta1-eta2))*DCOSH((n+1d0/2d0)*(eta1+eta2))
      bn = bn + 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(eta1-eta2))*DCOSH((n-1d0/2d0)*(eta1+eta2))
      bn = bn + 4d0*DEXP((eta2-eta1)*(n+1d0/2d0))*DSINH((n+1d0/2d0)*(eta1-eta2))
      bn = bn + (2d0*n+1d0)**2*DEXP(eta1-eta2)*DSINH(eta1-eta2)
      bn = bn *-(2d0*n+3d0)
      END

      FUNCTION cn(n,eta1,eta2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      cn = (2d0*n+1d0)*(n+3d0/2d0)*(DCOSH(2*eta1)-DCOSH(2*eta2))
      cn = cn + 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(eta1-eta2))*DSINH((n+1d0/2d0)*(eta1+eta2))
      cn = cn + 2d0*(2d0*n+1d0)*DSINH((n+3d0/2d0)*(eta1+eta2))*DSINH((n-1d0/2d0)*(eta2-eta1))
      cn = cn *-(2d0*n-1d0)
      END

      FUNCTION dn(n,eta1,eta2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      dn = (2d0*n+1d0)*(n+3d0/2d0)*(DSINH(2*eta1)-DSINH(2*eta2))
      dn = dn + 2d0*(2d0*n+3d0)*DSINH((n+1d0/2d0)*(eta1-eta2))*DCOSH((n+1d0/2d0)*(eta1+eta2))
      dn = dn + 2d0*(2d0*n+1d0)*DSINH((n-1d0/2d0)*(eta2-eta1))*DCOSH((n+3d0/2d0)*(eta1+eta2))
      dn = dn + 4d0*DEXP((eta2-eta1)*(n+1d0/2d0))*DSINH((n+1d0/2d0)*(eta1-eta2))
      dn = dn - (2d0*n+1d0)**2*DEXP(eta2-eta1)*DSINH(eta1-eta2)
      dn = dn * (2d0*n-1d0)
      END

      FUNCTION Deltan(n,eta1,eta2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      Deltan = 4d0*DSINH((n+1d0/2d0)*(eta1-eta2))**2
      Deltan = Deltan - ( (2d0*n+1d0) * DSINH(eta1-eta2) )**2
      Deltan = Deltan * (2d0*n-1d0)*(2d0*n+3d0) / ( n * (n+1d0) )
      END

      FUNCTION gma(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      gma = ( DTANH(al)*DTANH((n+5d-1)*al) )**(-1)
      END
! ======================================================================

! === G A U S S   E L I M I N A T I O N ================================
!     Written by Alexander Zinchenko (University of Colorado Boulder)

!     For every n, write the system in a matrix form AX=b,
!     where X=(A_n, B_n, C_n, D_n) and b is the RHS vector,
!     and solve this system by Gauss elimination with pivoting.
!     For solving a system of N linear equations for N unknowns with
!     M-N (M minus N) right-side vectors b (so, you can have solutions
!     simutaneously for matrix A and several RHS vectors b). The 
!     parameter (NMAX=100, MMAX=200) is just an example to make it
!     suitable for solving with <=100 unknowns with <= 200-100 RHS 
!     vectors. In your specific case (4 eqns with one RHS vector),
!     NMAX=4, MMAX=5 would suffice, but you do not have to make that
!     change, just make sure, one way or the other, that dimensions of
!     T are consistent in GAUSS and the calling routine.

!     To use GAUSS for AX=b, put A into the (N,N) part of matrix T,
!     (i.e. T_{iJ}=A_{ij} for i,j <=N) and set 
!     T(1,N+1)=b_1, T(2,N+1)=b_2, ... T(N,N+1)=b_N (if you have one
!     RHS vector b). After call GAUSS(N,N+1,T) you will get the 
!     solution X_i= T(i,N+1) for i<=N. In your case, N=4.

!     If you want the solutions for several different b-vectors with 
!     the same matrix A, Then these vectors must be placed into the
!     columns (N+1), (N+2), ... M of matrix T before calling GAUSS,
!     and the solutions will then be found in those columns after GAUSS.
!     Of course, all calcs must be in DOUBLE PRECISION.

      SUBROUTINE GAUSS(N,M,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     PARAMETER (NMAX=100, MMAX=200)
!     DIMENSION T(NMAX,MMAX)
      DIMENSION T(N,M) ! My modification
      INTEGER S
      
      DO 1 I=1,N
         S=I
         R=T(I,I)
         DO 2 K=I+1,N
            IF(DABS(T(K,I)).GT.DABS(R)) THEN
              S=K
              R=T(K,I)
            ENDIF
2        CONTINUE
         T(S,I)=T(I,I)
         DO 3 J=I+1,M
            U=T(S,J)/R
            T(S,J)=T(I,J)
            T(I,J)=U
            DO 4 K=I+1,N
4              T(K,J)=T(K,J)-T(K,I)*U
3        CONTINUE
1     CONTINUE
      DO 5 I=N-1,1,-1
         DO 5 J=N+1,M
            DO 5 K=I+1,N
               R=T(I,K)
               T(I,J)=T(I,J)-R*T(K,J)
5     CONTINUE
      RETURN
      END SUBROUTINE
! ======================================================================

! === T R I D I A G O N A L   M A T R I X   A L G O R I T H M ==========
!     Four vectors a, b, c and d are put into a single array T(a|b|c|d)
!     with four columns: N*4. The solution is returned in the 4th
!     column, overwritten on d.
      SUBROUTINE TDMA(N,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,4)

      DO I = 2, N
         T(I,2) = T(I,2) - T(I,1) / T(I-1,2) * T(I-1,3)
         T(I,4) = T(I,4) - T(I,1) / T(I-1,2) * T(I-1,4)
      ENDDO

      T(N,4) = T(N,4) / T(N,2)
      DO I = N-1, 1, -1
         T(I,4) = ( T(I,4) - T(I,3) * T(I+1,4) ) / T(I,2)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === G A U S S — S E I D E L ==========================================
!     The system of N equations, Ax = b, is given in a single array T(A|b|x)
!     where A is T(1:N,1:N), b is T(1:N,N+1), and the initial condition
!     until reaching final solution, with residual R, is the last column
!     T(1:N,N+2).
      SUBROUTINE GAUSS_SEIDEL(N,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,N+2)

      RF = 5d-1 ! Relaxation Factor
1     DO I = 1, N
         IF (ABS(T(I,I)).LT.1D-8) STOP "0 on the main diagonal!"
         S = 0D0
         DO J = 1, N
            IF (I.NE.J) S = S + T(I,J) * T(J,N+2)
         ENDDO
         T(I,N+2) = ( 1d0 - RF ) * T(I,N+2) + RF * ( T(I,N+1) - S ) / T(I,I)
      ENDDO

      R = 0D0
      DO I = 1, N
         S = 0D0
         DO J = 1, N
            S = S + T(I,J) * T(J,N+2)
         ENDDO
         R = R + ABS( T(I,N+1) - S )
      ENDDO
      WRITE(*,*) "Residual =", R

      IF (R.GT.1D+21) STOP "The system of equations is not diagonally &
                            dominant or it is not positive definite!"
      IF (R.GT.1D-6) GOTO 1

      END SUBROUTINE
! ======================================================================

! === G A U S S   E L I M I N A T I O N   B A N D E D   M A T R I X ====
!     KL = Lower band: No. of sub-diagonals
!     KU = Upper band: No. of super-diagonals
!     Example: N = 5, KL = 1, KU = 2
!     X  X  X  0  0 | X
!     X  X  X  X  0 | X
!     0  X  X  X  X | X
!     0  0  X  X  X | X
!     0  0  0  X  X | X
      SUBROUTINE GAUSSB(N,KL,KU,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,N+1)

      DO K = 1, N-1
        UK = T(K,N+1) / T(K,K)
        IJ = 1
        NI = K + KL
        IF ( NI .GT. N ) NI = N
        DO I = K+1, NI
          UI = T(I,K) / T(K,K)
          NJ = K + KU
          IF ( NJ .GT. N ) NJ = N
          DO J = K+1, NJ
            T(I,J) =  T(I,J)  - T(K,J) * UI
          ENDDO
          T(I,N+1) = T(I,N+1) - T(I,K) * UK
        ENDDO
      ENDDO

      DO I = N, 1, -1
        NJ = I + KU
        IF ( NJ .GT. N ) NJ = N
         S = 0D0
         DO J = I+1, NJ
            S = S + T(I,J) * T(J,N+1)
         ENDDO
         T(I,N+1) = ( T(I,N+1) - S ) / T(I,I)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === T H O M A S   A L G O R I T H M   B A N D E D   M A T R I X ======
!     KL = Lower band: No. of sub-diagonals
!     KU = Upper band: No. of super-diagonals
!     If KL = KU = 1 then the solver works
!     similar to TDMA. The system of equations
!     has to be given to the solver in the
!     following compact form:
!     beginning from the left-most column
!     we fill T(:,j) with vectors containing
!     sub-diagonal, diagonal, super-diagonal
!     and finally the RHS (vector b) elements.
!     Example: N = 5, KL = 1, KU = 2
!     2  3  4  0  0 | 5
!     1  2  3  4  0 | 5
!     0  1  2  3  4 | 5
!     0  0  1  2  3 | 5
!     0  0  0  1  2 | 5
!     This system has to be rearranged to:
!     0  2  3  4 | 5
!     1  2  3  4 | 5
!     1  2  3  4 | 5
!     1  2  3  0 | 5
!     1  2  0  0 | 5

      SUBROUTINE THOMAS(N,KL,KU,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,KL+KU+2)

      DO K = 1, N-1
        NI = K + KL
        IF ( NI .GT. N ) NI = N
        DO I = K+1, NI
           U = T(I, K+KL-I+1) / T(K, KL+1)
           IF ( ABS(T(K, KL+1)) .LT. 1D-15 ) &
           WRITE(*,*) 'Check: diagonal element = 0'
           NJ = K + KL + KU - I + 1
           DO J = K+KL-I+2, NJ
              T(I,J) = T(I,J) - T(K, I+J-K) * U
           ENDDO
           T(I, KL+KU+2) = T(I, KL+KU+2) - T(K, KL+KU+2) * U
        ENDDO
      ENDDO

      DO I = N, 1, -1
         K = I + 1
         S = 0D0
         DO J = KL+2, KL+KU+1
            IF (K.GT.N) EXIT
            S = S + T(I,J) * T(K, KL+KU+2)
            K = K + 1
         ENDDO
         T(I, KL+KU+2) = ( T(I, KL+KU+2) - S ) / T(I, KL+1)
      ENDDO

      END SUBROUTINE
! ======================================================================
