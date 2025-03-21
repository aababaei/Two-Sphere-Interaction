!     Subroutines for computing electrostatic interaction between two spheres
!     Developed by Ahmad Ababaei (ahmad.ababaei@imgw.pl) and Antoine Michel (antoine.michel@cea.fr)

      PROGRAM ELECTROFORCES
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PARAMETER ( acu = 1d-9 )
      INTEGER sample

      pi = 4d0*DATAN(1d0)
       e = 1.602176634d-19 ! C

! ============ I N P U T S ==============

        alam = 1.00d-0   ! Radius ratio
          a1 = 1d+0      ! Larger drop radius [μm]
          a2 = alam * a1 ! Smaller drop radius [μm]
          E0 = 0d-1      ! Electric field intensity [V/cm]
         psi = 0d0*pi    ! Electric field angle
          qr = 0.1d0       ! Electric charge ratio
          q1 = 200*e     ! Electric charge [C]
          q2 = qr * q1   ! Electric charge [C]
        epsr = 80d0      ! Dielectric constant / relative permittivity (water = 80, perfect conductor: ∞)

! ============= U N I T S ===============

          a1 = a1 * 1d-4 ! [μm] to [cm]
          a2 = a2 * 1d-4 ! [μm] to [cm]
          q1 = q1 * 2997924580d0 ! [C] to [statC]
          q2 = q2 * 2997924580d0 ! [C] to [statC]

! ========= L O G   D I S T R. ========== 
! Logarithmic distribution of normalized
! gap size ξ = s — 2 in JO84 notation:
      xi_min = 1d-2
      xi_max = 1d+2
      sample = 10
      dlt_xi = DLOG ( xi_max / xi_min ) / DBLE(sample-1)
           s = 2d0 + xi_min
        DO i = 1, sample
           r = s * ( a1 + a2 ) / 2d0

! ============ M E T H O D ==============

!     CALL COULOMB(q1,q2,r,F12)
!     CALL D64(a1,a2,r,q1,q2,E0,psi,F1,F2,acu)
      CALL KCSB14(a1,a2,r,q1,q2,epsr,F12,acu)
!     CALL BPBS21(a1,a2,r,q1,q2,F12,acu)


! ============ O U T P U T ==============
      WRITE(*,*) s-2d0, F1, F2, F12
!     STOP

! ========= L O G   D I S T R. ==========
      s = DLOG ( s - 2d0 ) + dlt_xi
      s = DEXP ( s ) + 2d0
      ENDDO

      END PROGRAM ELECTROFORCES

! ======================= S U B R O U T I N E S ========================

! === Coulomb (1785) ===================================================

      SUBROUTINE COULOMB(q1,q2,r,F12)
      IMPLICIT DOUBLE PRECISION (A-Z)
      F12 = q1 * q2 / r**2
      END SUBROUTINE COULOMB

! ======================================================================

! === Davis (1964) =====================================================
!     Davis, M. H. (1964). Two charged spherical conductors in a uniform electric field: Forces and field strength. The Quarterly Journal of Mechanics and Applied Mathematics, 17(4), 499-511.

      SUBROUTINE D64(a1,a2,r,q1,q2,E0,psi,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,L-Z)

!     eps0 = 8.854187817d-12 * 2.99792458d9**2 * 1d-6 ! esu
      eps0 = 1d0

      IF ( DABS(E0) .LT. 1d-10 ) E0 = 1d-10
      Ex = E0*DSIN(psi)
      Ez = E0*DCOS(psi)

          h = r - ( a1 + a2 )
        eps =  h / a1
       alam = a2 / a1
     coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
     coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
        mu1 = DACOSH(coshal)
        mu2 = DACOSH(coshbe)
        m12 = mu1 + mu2
          c = a1 * DSINH(mu1)

        Qs1 = 2d0 * eps0 * c**2 * E0 * DCOS(psi) * ( S(1,mu2,m12,acu) + S(1,0d0,m12,acu) )
        Qs2 =-2d0 * eps0 * c**2 * E0 * DCOS(psi) * ( S(1,mu1,m12,acu) + S(1,0d0,m12,acu) )
        C11 = 2d0 * eps0 * c * S(0,mu2,m12,acu)
        C12 =-2d0 * eps0 * c * S(0,0d0,m12,acu)
        C22 = 2d0 * eps0 * c * S(0,mu1,m12,acu)
        del = C11 * C22 - C12**2
        P11 = C22 / del
        P12 =-C12 / del
        P22 = C11 / del
         v1 =-( P11*Qs1 + P12*Qs2 ) / ( E0 * c * DCOS(psi) )
         v2 =-( P12*Qs1 + P22*Qs2 ) / ( E0 * c * DCOS(psi) )
         w1 = ( P11*(Q1-Qs1) + P12*(Q2-Qs2) ) / ( E0 * c * DCOS(psi) )
         w2 = ( P12*(Q1-Qs1) + P22*(Q2-Qs2) ) / ( E0 * c * DCOS(psi) )
        p11 = c * P11
        p12 = c * P12
        p22 = c * P22

       F2zo = 0d0
       F2xo = 0d0
        rel = 1d0
         n  = 0d0
      DO WHILE ( rel .GT. acu )
         Yn =-DSQRT(2d0)*(2d0*n+1d0)*DEXP((n+5d-1)*mu2)
         Yn = Yn * ( (2d0*n+1d0)*(DEXP((2d0*n+1d0)*mu1)+1d0) - w2*DEXP((2d0*n+1d0)*mu1) + w1 )
         Yn = Yn / ( DEXP((2d0*n+1d0)*m12) - 1d0 )
        Ynp =-DSQRT(2d0)*(2d0*n+3d0)*DEXP((n+15d-1)*mu2)
        Ynp = Ynp* ( (2d0*n+3d0)*(DEXP((2d0*n+3d0)*mu1)+1d0) - w2*DEXP((2d0*n+3d0)*mu1) + w1 )
        Ynp = Ynp/ ( DEXP((2d0*n+3d0)*m12) - 1d0 )
         Zn = DSQRT(8d0)*(2d0*n+1d0)*DEXP((n+05d-1)*mu2) * (DEXP((2d0*n+1d0)*mu1)-1d0) / ( DEXP((2d0*n+1d0)*m12) - 1d0 )
        Znp = DSQRT(8d0)*(2d0*n+3d0)*DEXP((n+15d-1)*mu2) * (DEXP((2d0*n+3d0)*mu1)-1d0) / ( DEXP((2d0*n+3d0)*m12) - 1d0 )
        F2z = 2d0 * DCOS(psi)**2 * Yn / (2d0*n+1d0) * ( Yn - 2d0*DCOSH(mu2)*(n+1d0)/(2d0*n+3d0)*Ynp )
        F2z = F2z + DSIN(psi)**2 * n*(n+1d0)/(2d0*n+1d0)*Zn * ( Zn - 2d0*DCOSH(mu2)*(n+2d0)/(2d0*n+3d0)*Znp )
        F2z = F2zo + F2z * eps0/4d0*(c*E0)**2
        F2x = (n+1d0)/(2d0*n+1d0)/(2d0*n+3d0) * ( (n+2d0)*Znp*Yn - n*Zn*Ynp )
        F2x = F2xo + F2x * eps0/4d0*(c*E0)**2 * DSIN(2d0*psi) * DSINH(mu2)
        rel = DMAX1( DABS(F2z-F2zo)/DABS(F2z), DABS(F2x-F2xo)/DABS(F2x) )
       F2zo = F2z
       F2xo = F2x
          n = n + 1d0
      ENDDO
      F1z = Ez * ( q1 + q2 ) - F2z
      F1x = Ex * ( q1 + q2 ) - F2x
      F1  = F1z
      F2  = F2z

      END SUBROUTINE D64

      FUNCTION S(m,xi,m12,acu)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER m
       So = 0d0
        S = 0d0
        n = 0d0
      rel = 1d0
      DO WHILE (rel .GT. acu)
         nn1 = 2d0*n + 1d0
           S = So + nn1**m * DEXP(nn1*xi) / ( DEXP(nn1*m12) - 1d0 )
         rel = DABS(S-So)/DABS(S)
          So = S
           n = n + 1d0
      ENDDO
      END FUNCTION

! ======================================================================

! === Khachatourian, Chan, Stace, Bichoutskaia (2014) ==================
!     Khachatourian, A., Chan, H. K., Stace, A. J., & Bichoutskaia, E. (2014). Electrostatic force between a charged sphere and a planar surface: A general solution for dielectric materials. The Journal of chemical physics, 140(7).

      SUBROUTINE KCSB14(a1,a2,r,q1,q2,epsr,F12,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

          k0 = 1d0  ! air relative permittivity
          k1 = epsr
          k2 = epsr
          km = 1d0

           h = r - ( a1 + a2 )
         eps =  h / a1
        alam = a2 / a1
      coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
      coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
         et1 = DACOSH(coshal)
         et2 = DACOSH(coshbe)
           c = a1 * DSINH(et1)

      F12o = 0d0
!     F21o = 0d0
        iN = 10
1     ALLOCATE ( T(2*iN+2,8) )
         T = 0d0

      DO i = 1, iN+1
         n = DBLE(i-1)
        fm = DEXP(-(n-05d-1)*(et1+et2))
        fn = DEXP(-(n+05d-1)*(et1+et2))
        fp = DEXP(-(n+15d-1)*(et1+et2))

         IF (i.GT.1) THEN
            T(2*i-1,2) = -5d-1*n*(km+k1)
            T(2*i-1,3) =  5d-1*n*(km-k1)*fm
            T(2*i  ,1) =  5d-1*n*(km-k2)*fm
            T(2*i  ,2) = -5d-1*n*(km+k2)
         ENDIF
            T(2*i-1,4) = (5d-1+n)*DCOSH(et1)*(km+k1) + 5d-1*DSINH(et1) * (km-k1)
            T(2*i-1,5) = (-(5d-1+n)*DCOSH(et1)+5d-1*DSINH(et1))*(km-k1)* fn
            T(2*i  ,3) = (-(5d-1+n)*DCOSH(et2)+5d-1*DSINH(et2))*(km-k2)* fn
            T(2*i  ,4) = (5d-1+n)*DCOSH(et2)*(km+k2) + 5d-1*DSINH(et2) * (km-k2)
         IF (i.LT.iN+1) THEN
            T(2*i-1,6) = -5d-1*(n+1d0)*(km+k1)
            T(2*i-1,7) =  5d-1*(n+1d0)*(km-k1)*fp
            T(2*i  ,5) =  5d-1*(n+1d0)*(km-k2)*fp
            T(2*i  ,6) = -5d-1*(n+1d0)*(km+k2)
         ENDIF
!           K factor not needed in ESUnits:
            T(2*i-1,8) = DSQRT(2d0)*c*DEXP(-(n+5d-1)*et1)*q1/a1**2
            T(2*i  ,8) = DSQRT(2d0)*c*DEXP(-(n+5d-1)*et2)*q2/a2**2
      ENDDO

      CALL THOMAS(2*iN+2,3,3,T)
!        n = 0:
        fn = DEXP(-5d-1*(et1+et2))
       F12 = fn * 5d-1 * ( -T(1,8) + T(3,8)*DEXP(-et1) ) * T(2,8)
!      F21 = fn * 5d-1 * ( -T(2,8) + T(4,8)*DEXP(-et2) ) * T(1,8)
!        n = 1, N:
      DO i = 1, iN+1
         n = DBLE(i)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*T(2*i-1,8)*DEXP( et1) - (n+5d-1)*T(2*i+1,8) &
                 +  (n+1d0)/2d0*T(2*i+3,8)*DEXP(-et1) ) * T(2*i+2,8)
!      F21 = F21 + fn * ( n/2d0*T(2*i  ,8)*DEXP( et2) - (n+5d-1)*T(2*i+2,8) &
!                +  (n+1d0)/2d0*T(2*i+4,8)*DEXP(-et2) ) * T(2*i+1,8)
      ENDDO
      DEALLOCATE ( T )

!     Force in SI:
!     F12 = -F12/k
!     F21 = -F21/k
!     Force in dynes (K factor not needed in ESUnits):
      F12 = -F12
!     F21 = -F21
      rel = DABS(F12-F12o)/DABS(F12)
      IF ( rel .GT. acu ) THEN
           iN = INT(1.5 * FLOAT(iN)) ! 50% increase
         F12o = F12
!        F21o = F21
         WRITE(*,*) 'n_max, F12, rel = ', iN,F12,rel
         GOTO 1
      ENDIF

      END SUBROUTINE KCSB14
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
         S = 0D0
         DO J = KL+2, KL+KU+1
            K = I + J - ( KL + 1 )
            IF ( K .GT. N ) EXIT
            S = S + T(I,J) * T(K, KL+KU+2)
         ENDDO
         T(I, KL+KU+2) = ( T(I, KL+KU+2) - S ) / T(I, KL+1)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === Banerjee, Peters, Brown, Song (2021) =============================
!     Banerjee, S., Peters, T., Brown, N., & Song, Y. (2021). Exact closed-form and asymptotic expressions for the electrostatic force between two conducting spheres. Proceedings of the Royal Society A, 477(2246), 20200866.
      SUBROUTINE BPBS21(a1,a2,r,q1,q2,F12,acu)
      IMPLICIT NONE
      DOUBLE PRECISION :: a1,a2,r,q1,q2,F12,acu
      DOUBLE PRECISION :: eulc,sep,rad,alp,muc,lam,xc,yc
      DOUBLE PRECISION :: dalpds,dmucds,dlamds
      DOUBLE PRECISION :: dxcdalp,dxcdmuc,dycdalp,dycdmuc
      DOUBLE PRECISION :: c11a,c12a,c22a,digamxc,digamyc
      DOUBLE PRECISION :: sbc11,sbc12,sbc22,sbc11a,sbc22a
      DOUBLE PRECISION :: dc11ds,dc12ds,dc22ds
      DOUBLE PRECISION :: c11l,c12l,c22l
      DOUBLE PRECISION :: c11c,c12c,c22c,digax,digaax,digay,digaay,digaa
      DOUBLE PRECISION :: c11,c22,c12,p11,p22,p12
      DOUBLE PRECISION :: V1,V2,VR,fv,qr,fr,fq
      DOUBLE PRECISION :: ddgadx,ddgady,xc0,yc0,q0,psxc0,psyc0
      DOUBLE PRECISION :: ddga2xda2,ddgaxda,ddga2xdx,ddgaxdx
      DOUBLE PRECISION :: ddga2yda2,ddgayda,ddga2ydy,ddgaydy
      DOUBLE PRECISION :: ddga12da2
      DOUBLE PRECISION :: dc11dsl,dc22dsl,dc12dsl,fvl
      DOUBLE PRECISION :: dc11dsc,dc22dsc,dc12dsc,fvc
      DOUBLE PRECISION :: dc11dsa,dc22dsa,dc12dsa,fva
      INTEGER :: i,ii,j,k
      ! this is implementation of the formulas provided by Banerje et
      ! al. Proceedings Royal Society A, 2021
      ! 'Exact closed-form and asymptotic expressions for the
      ! electrostatic force between two conducting spheres'

      ! define parameters
      eulc=0.577215664901532860606512090082 ! euler constant
      sep =r/(a1+a2)        ! normalized separation
      rad =(a1-a2)/(a1+a2)  ! asymetry
      alp =( (sqrt(sep**2-rad**2)-sqrt(sep**2-1.))/sqrt(1.-rad**2) )**2 !maps the separation distance into [0,1] interval
      muc = -log( sqrt(alp) ) ! another distance variable 
      lam =(1.-rad**2)/(2.*sep)*(1.-alp**2)/alp ! simplifies capacitance expression later on
      xc=1./2.-atanh(rad*tanh(log(sqrt(alp))))/log(alp)
      yc=1./2.+atanh(rad*tanh(log(sqrt(alp))))/log(alp)

      qr=q2/q1 ! charge ratio
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! A. compute the capacitance coefficients 
      ! different formulations are used depending on the separation
      ! distance
      ! 1. Lambert series formulation, at large separations
      IF ( muc.gt. 0.95 ) THEN
       c11l=0.
       c22l=0.
       c12l=0.
       DO i=0,10 ! Banerjee: error of 10^-10
        c11l=c11l+lam*(alp**(DBLE(i)+xc))/(1.-alp**(2.*DBLE(i)+2.*xc))
        c22l=c22l+lam*(alp**(DBLE(i)+yc))/(1.-alp**(2.*DBLE(i)+2.*yc))
       IF(i.gt.0)c12l=c12l-lam*(alp**(DBLE(i)) )/( 1.-alp**(2.*DBLE(i)))
       ENDDO
       c11=c11l
       c22=c22l
       c12=c12l
      ENDIF

      ! 2. Closed form formulation in middle range
      IF ( (muc.gt.0.35) .and. (muc.le.0.95) )THEN
       CALL QDIGAMMA(xc,alp   ,digax) ! implem ok
       CALL QDIGAMMA(xc,alp**2,digaax)
       CALL QDIGAMMA(yc,alp   ,digay)
       CALL QDIGAMMA(yc,alp**2,digaay)
       CALL QDIGAMMA(5.d-1,alp**2,digaa)
       c11c=lam/(4.*muc) * (digaax-2.*digax+dlog((1.+alp)/(1.-alp)) )
       c22c=lam/(4.*muc) * (digaay-2.*digay+dlog((1.+alp)/(1.-alp)) )
       c12c=lam/(4.*muc) * (digaa+dlog(1.-alp**2) )
       c11=c11c
       c22=c22c
       c12=c12c
      ENDIF

      ! 3. Asymptotic form for small separation
      IF ( muc.le. 0.35 ) THEN
       CALL DIGAMMA(xc,digamxc) ! implem ok
       CALL DIGAMMA(yc,digamyc)
       CALL SBERNOUC(muc,xc,yc,sbc11,sbc22,sbc12)
       c11a= lam/(4.*muc)*(dlog(1./muc)-digamxc-sbc11)
       c22a= lam/(4.*muc)*(dlog(1./muc)-digamyc-sbc22)
       c12a=-lam/(4.*muc)*(dlog(1./muc)+eulc   -sbc12)
       c11=c11a
       c22=c22a
       c12=c12a
      ENDIF

      ! Test notes: Implementation yields expected results
      ! lambert series  = closed form at long  separation distances
      ! asymptotic form = closed form at short separation distances

      ! The potential coefficients matrix is defined as the inverse of
      ! the capacitance coefficients matrix  
      p11= 1./(c11*c22-c12*c12)*c22
      p22= 1./(c11*c22-c12*c12)*c11
      p12=-1./(c11*c22-c12*c12)*c12

      ! Use Eq. 48 in Banerjee (2020) to determine the voltage ratio at
      ! constant charge
      xc0=(1.-rad)/2.
      yc0= 1.-xc0
      CALL DIGAMMA (xc0,psxc0)
      CALL DIGAMMA (yc0,psyc0)
      q0 =(eulc+psyc0)/(eulc+psxc0)
      !qr = q0*.6 ! Fix qr value to check implementation => OK
      vr = (c11*qr-c12)/(c22-c12*qr)

      ! end of matrix coefficients computation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! B. determine the force at constant voltage
      ! f_v = F_v/(pi eps0 V_1^2) = dc11/ds + 2v dc12/ds + v^2 dc22/ds

      dalpds=-4.*alp/lam
      dmucds=2./lam
      dlamds=4.*(1.+alp**2)/(1.-alp**2)-lam/sep
      dxcdalp=-(1.-2.*xc)/(4.*alp*muc) + rad*COSH( ((1.-2.*xc)*muc) )**2 /( (1.+alp)**2 *muc )
      dxcdmuc=-2.*alp*dxcdalp
      dycdalp=-dxcdalp
      dycdmuc=-dycdmuc

      ! 1. lambert series capacitance derivatives
      IF ( muc.gt.0.95 ) THEN
       dc11dsl=c11/lam*dlamds
       dc22dsl=c22/lam*dlamds
       dc12dsl=c12/lam*dlamds
       DO i=0,10
        dc11dsl=dc11dsl-4.*(alp**(3.*DBLE(i)+3.*xc)+alp**(DBLE(i)+xc) ) &
               *(DBLE(i)+xc+alp*dxcdalp*log(alp))/(alp**(2.*DBLE(i)+2.*xc)-1.)**2
        dc22dsl=dc22dsl-4.*(alp**(3.*DBLE(i)+3.*yc)+alp**(DBLE(i)+yc) ) &
               *(DBLE(i)+yc+alp*dycdalp*log(alp))/(alp**(2.*DBLE(i)+2.*yc)-1.)**2
        IF(i.gt.0)THEN
         dc12dsl=dc12dsl+4.*( alp**(3.*DBLE(i))+alp**(DBLE(i)) )*DBLE(i)/(alp**(2.*DBLE(i))-1.)**2
        ENDIF
       ENDDO
       fvl=dc11dsl+2.*vr*dc12dsl+vr**2*dc22dsl
       fv =fvl
      ENDIF

      ! 2. closed form derivatives
      IF ( muc.le.0.95 ) THEN !(muc.gt.0.35) .and. (muc.le.0.95) )THEN
        CALL DQDIGAMMADQ (xc,alp**2,ddga2xda2)
        CALL DQDIGAMMADQ (xc,alp   ,ddgaxda  )
        CALL DQDIGAMMADX (xc,alp**2,ddga2xdx )
        CALL DQDIGAMMADX (xc,alp   ,ddgaxdx  )
        CALL DQDIGAMMADX (xc,alp   ,ddgaxdx  )
        CALL DQDIGAMMADQ (yc,alp**2,ddga2yda2)
        CALL DQDIGAMMADQ (yc,alp   ,ddgayda  )
        CALL DQDIGAMMADX (yc,alp**2,ddga2ydy )
        CALL DQDIGAMMADX (yc,alp   ,ddgaydy  )
        CALL DQDIGAMMADQ (5.d-1,alp**2,ddga12da2)
        dc11dsc=-2.*alp/muc*( alp*ddga2xda2-ddgaxda+1./(1.-alp**2) )       &
               -c11*( 2./lam/muc+1./sep-4.*COSH(2.*muc)/SINH(2.*muc)/lam ) &
               -alp*dxcdalp/muc*(ddga2xdx-2.*ddgaxdx)
        dc22dsc=-2.*alp/muc*( alp*ddga2yda2-ddgayda+1./(1.-alp**2) )       &
               -c22*( 2./lam/muc+1./sep-4.*COSH(2.*muc)/SINH(2.*muc)/lam ) &
               -alp*dycdalp/muc*(ddga2ydy-2.*ddgaydy)
        dc12dsc=-2.*alp**2/muc*(ddga12da2-1./(1.-alp**2))                  &
               -c12*( 2./lam/muc+1./sep-4.*COSH(2.*muc)/SINH(2.*muc)/lam )
        fvc=dc11dsc+2.*vr*dc12dsc+vr**2*dc22dsc
        fv=fvc
      ENDIF

      ! 3. asymptotic form derivatives
      ! After testing I noticed that the asymptotic formulations does not yield
      ! accurate results in the event of vr or qr near 1 so I chose to
      ! use the closed form formula instead
      IF ( 1.EQ.0 ) THEN !muc.LE.0.35 )THEN
       CALL SBERNOUCD(muc,xc,yc,sbc11a,sbc11,sbc22a,sbc22,sbc12)
       CALL DDIGAMMADX (xc,ddgadx)
       CALL DDIGAMMADX (yc,ddgady)
       dc11dsa=-1./(2.*muc**2)+c11/lam*(dlamds-2./muc)-dxcdmuc/(2.*muc) *(ddgadx + sbc11a) - sbc11
       dc22dsa=-1./(2.*muc**2)+c22/lam*(dlamds-2./muc)-dycdmuc/(2.*muc) *(ddgady + sbc22a) - sbc22
       dc12dsa= 1./(2.*muc**2)+c12/lam*(dlamds-2./muc)+sbc12
       fva=dc11dsa+2.*vr*dc12dsa+vr**2*dc22dsa
      ENDIF
      fq = (p11+p12*qr)**2/q0*fv
      ! dimensional force
      ! ommit 4*pi*eps0 factor for cgs units
      fq = fq*q0*q1**2 / (a1+a2)**2
      F12=fq
      !WRITE(*,*)rad,qr/q0,'tut',sep-1,muc,fq,(p11+p12*qr)**2/q0*fva

      END SUBROUTINE BPBS21

! ======================================================================

      SUBROUTINE SBERNOUCD (muc,xc,yc,sbc11a,sbc11b,sbc22a,sbc22b,sbc12)
      IMPLICIT NONE
      DOUBLE PRECISION :: muc,xc,yc
      DOUBLE PRECISION :: sbc11a,sbc11b,sbc22a,sbc22b,sbc12,facto
      DOUBLE PRECISION :: n0a,n0b,d0,n111,n112,n221,n222,n121,n122,d1
      DOUBLE PRECISION :: b2kxc ,b2kmxc ,b2kyc ,b2kmyc ,b2k12 ,binnk,b2k
      DOUBLE PRECISION :: b2kxct,b2kmxct,b2kyct,b2kmyct,b2k12t,binnkm
      DOUBLE PRECISION, DIMENSION(:) :: bern(11)
      INTEGER :: k,n,exp2k,exp2km
      ! Sum appearing in equs 32 and 33 in banerjee 2020
      ! calls for bernoulli polynomials and bernoulli's number
      ! needed to compute the capacitance coefficients
      bern(:) =0.
      bern( 1)= 1.    ! careful indexes start at 0
      bern( 2)=-1./2.
      bern( 3)= 1./6.
      bern( 5)=-1./30.
      bern( 7)= 1./42.
      bern( 9)=-1./30.
      bern(11)= 5./66.

      b2kxct =0.
      b2kmxct=0.
      b2kyct =0.
      b2kmyct=0.
      b2k12t =0.

      DO k=1,5 !5 ! sufficient for accurate results according to authors
        b2kxc =0.
        b2kmxc=0.
        b2kyc =0.
        b2kmyc=0.
        b2k12 =0.
        DO n=0,2*k ! bernoulli polynomials
          binnk =facto(2*k  )/facto(n)/facto(2*k  -n) ! binomial coeff 2k
          binnkm=facto(2*k-1)/facto(n)/facto(2*k-1-n) ! binomial coeff 2k-1
                        b2kxc =b2kxc +binnk *bern(2*k+1-n)*(xc**n) ! B_2k(x)
          IF(n.lt.(2*k))b2kmxc=b2kmxc+binnkm*bern(2*k  -n)*(xc**n) ! B_2k-1(x)
                        b2kyc =b2kyc +binnk *bern(2*k+1-n)*(yc**n) ! B_2k(y)
          IF(n.lt.(2*k))b2kmyc=b2kmyc+binnkm*bern(2*k  -n)*(yc**n) ! B_2k-1(y)
                        b2k12 =b2k12 +binnk *bern(2*k+1-n)*(.5**n) ! B_2k(1/2)
        ! implem ok. Bernoulli number when x=0 and right results for x=1
        ENDDO
        b2k= bern(2*k+1)
        n0a=2.**(4*k)
        n0b=2.**(4*k-1)
        d0 =facto(2*k)
        ! expression for the sums in eqs 32 and 33
        b2kxct =b2kxct +n0b*b2kxc *b2k12/d0*muc**(2*k-2)
        b2kmxct=b2kmxct+n0a*b2kmxc*b2k12/d0*muc**(2*k)
        b2kyct =b2kyct +n0b*b2kyc *b2k12/d0*muc**(2*k-2)
        b2kmyct=b2kmyct+n0a*b2kmyc*b2k12/d0*muc**(2*k)
        b2k12t =b2k12t +n0b*b2k   *b2k12/d0*muc**(2*k-2)
      ENDDO
      sbc11a=b2kmxct
      sbc11b=b2kxct
      sbc22a=b2kmyct
      sbc22b=b2kyct
      sbc12 =b2k12t

      END SUBROUTINE SBERNOUCD

! ======================================================================

      SUBROUTINE SBERNOUC (muc,xc,yc,sbc11,sbc22,sbc12)
      IMPLICIT NONE
      DOUBLE PRECISION :: muc,xc,yc,sbc11,sbc22,sbc12,facto
      DOUBLE PRECISION :: n0,d0,n111,n112,n221,n222,n121,n122,d1
      DOUBLE PRECISION :: b2kxc ,b2kyc ,b2k12 ,binnk
      DOUBLE PRECISION :: b2kxct,b2kyct,b2k12t
      DOUBLE PRECISION, DIMENSION(:) :: bern(11)
      INTEGER :: k,n
      ! Sum appearing in equs 20, 21 and 22 in banerjee 2020
      ! calls for bernoulli polynomials and bernoulli's number
      ! needed to compute the capacitance coefficients
      bern(:) =0.
      bern( 1)= 1.    ! careful indexes start at 0
      bern( 2)=-1./2.
      bern( 3)= 1./6.
      bern( 5)=-1./30.
      bern( 7)= 1./42.
      bern( 9)=-1./30.
      bern(11)= 5./66.

      b2kxct=0.
      b2kyct=0.
      b2k12t=0.
      DO k=1,5 ! sufficient for accurate results according to authors
        b2kxc=0.
        b2kyc=0.
        b2k12=0.
        DO n=0,2*k ! bernoulli polynomials
          binnk=facto(2*k)/facto(n)/facto(2*k-n)
          b2kxc=b2kxc+binnk*bern(2*k+1-n)*(xc**n) ! B_2k(x)
          b2kyc=b2kyc+binnk*bern(2*k+1-n)*(yc**n) ! B_2k(y)
          b2k12=b2k12+binnk*bern(2*k+1-n)*(.5**n) ! B_2k(1/2)
        ENDDO
        n0=2.**(4*k-1)
        d0=facto(2*k)*k
        b2kxct=b2kxct+n0*b2kxc*b2k12/d0*muc**(2*k)
        b2kyct=b2kyct+n0*b2kyc*b2k12/d0*muc**(2*k)
        b2k12t=b2k12t+n0*bern(2*k+1)*b2k12/d0*muc**(2*k)
      ENDDO
      sbc11=b2kxct
      sbc22=b2kyct
      sbc12=b2k12t
      END SUBROUTINE

! ======================================================================

      SUBROUTINE DIGAMMA (x,digammax)
      IMPLICIT NONE
      DOUBLE PRECISION :: x,eulc,digammax,inc
      INTEGER :: n
      ! derivative of the logarithm of the gamma function
      eulc=0.577215664901532860606512090082 !Euler's constant
      digammax=-eulc
      DO n=1,10000 ! increase to improve convergence
        inc=(1./DBLE(n)-1./(DBLE(n)+x-1.))
        digammax=digammax+inc
      ENDDO
      END SUBROUTINE DIGAMMA

! ======================================================================

      SUBROUTINE DDIGAMMADX (x,ddigamdx)
      IMPLICIT NONE
      DOUBLE PRECISION :: x,eulc,facto,ddigamdx,inc
      INTEGER :: n
      ! derivative of digamma function with respect to x
      eulc=0.577215664901532860606512090082 !Euler's constant
      ddigamdx=0.
      DO n=0,10000 ! increase to improve convergence
        inc =1./(x+DBLE(n))**2
        ddigamdx=ddigamdx+inc ! formula ok
      ENDDO
      END SUBROUTINE DDIGAMMADX

! ======================================================================

      SUBROUTINE QDIGAMMA (x,q,qdigammax)
      IMPLICIT NONE
      DOUBLE PRECISION :: x,q,eulc,qdigammax
      DOUBLE PRECISION :: qdi0,qdin,inc
      INTEGER :: n
      ! q-digamma function (series representation)
      ! Verified implem :: OK
      qdi0=0.
      qdin=0.
      IF( q.GT.0. .AND. q.LT.1) THEN
        qdi0=-log(1.-q)
        DO n=0,10000 !
          inc= ( q**( DBLE(n)+x) )/( 1.-q**( DBLE(n)+x) )
          qdin=qdin+inc
        ENDDO
        qdigammax=qdi0+log(q)*qdin
      ENDIF
      IF( q.GT.1.) THEN
        qdi0=-log(q-1.)
        DO n=0,10000 !
          inc=-( q**(-DBLE(n)-x) )/( 1.-q**(-DBLE(n)-x) )
          qdin=qdin+inc
        ENDDO
        qdigammax=qdi0+log(q)*(x-.5)+log(q)*qdin
      ENDIF
      END SUBROUTINE QDIGAMMA

! ======================================================================

      SUBROUTINE DQDIGAMMADQ (x,q,qdigammax)
      IMPLICIT NONE
      DOUBLE PRECISION :: x,q,eulc,qdigammax
      DOUBLE PRECISION :: qdi0,qdin,inc1,inc2,inc3
      DOUBLE PRECISION :: ex,den
      INTEGER :: n
      ! derivative of the q-digamma function with respect to q
      qdi0=0.
      inc1=0.
      IF (q.LT.1) THEN
        qdi0=1./(1.-q)
        DO n=0,10000 ! increase to improve convergence
          ex  = DBLE(n)+x
          den =(1.-q**ex)**2
          inc1=inc1+q**(ex-1) * (ex*log(q)+1.-q**ex) / den
        ENDDO
      qdigammax=qdi0+inc1
      ENDIF
      IF (q.GT.1) THEN
        qdi0=-1./(q-1.)
        DO n=0,10000 ! increase to improve convergence
          ex  =-DBLE(n)-x
          den =(1.-q**ex)**2
          inc1=inc1+q**(ex-1) * (ex*log(q)+1.-q**ex) / den
        ENDDO
      qdigammax=qdi0 + x/q - 1./2./q + inc1
      ENDIF
      END SUBROUTINE DQDIGAMMADQ

! ======================================================================

      SUBROUTINE DQDIGAMMADX (x,q,qdigammax)
      IMPLICIT NONE
      DOUBLE PRECISION :: x,q,eulc,qdigammax,qdigammax2
      DOUBLE PRECISION :: qdi0,qdin,inc1,inc2
      DOUBLE PRECISION :: ex,den
      INTEGER :: n
      ! derivative of the q-digamma function with respect to x
      ! function def. differs if |q|<1 or |q|>1
      qdi0=0.
      inc1=0.
      IF (abs(q) .LT. 1) THEN ! implem ok
        DO n=0,10000
          inc1=inc1 + q**( DBLE(n)+x) / ( 1.-q**( DBLE(n)+x) )**2
        ENDDO
        inc1=inc1*log(q)**2
      ENDIF

      IF (abs(q) .GT. 1) THEN ! implem ok
        DO n=0,10000
          inc1=inc1 + q**(-DBLE(n)-x) / ( 1.-q**(-DBLE(n)-x) )**2
        ENDDO
        inc1=inc1*log(q)**2 + log(q)
      ENDIF

      qdigammax=qdi0+inc1
      END SUBROUTINE DQDIGAMMADX

! ======================================================================

      DOUBLE PRECISION FUNCTION FACTO (n)
      IMPLICIT NONE
      INTEGER :: i,n
      DOUBLE PRECISION :: fn
      facto=1.
      IF (n.eq.0) THEN
        facto=1.
      ENDIF
      IF (n.gt.0) THEN
        DO i=1,n
          facto=facto*DBLE(i)
        ENDDO
      ENDIF
      END FUNCTION FACTO

! ======================================================================
