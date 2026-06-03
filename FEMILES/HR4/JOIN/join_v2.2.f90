!version 2.2
!     +-----------------------------------------------------------------
!     | Program:        Join spectra - EMILES
!     | Author:         Elena Ricciardelli  
!     | Last modified:  07/11/2016;  A.Vazdekis 26/12/2023
!     +-----------------------------------------------------------------
! Compile: gfortran join_v2.2.f90 -o join_v2.2.x
! Exec: ./join_v2.2.x < join_v2.2
! Program to combine spectra in different spectral regions. After resampling and smoothing the input spectra
! the spectra are normalized according to the reference spectrum/spectra. 
! At the edges spectra are just re-normalize according to the continuum in the reference spectrum. Spectra in the middle
! are tilted forcing the continuum to follow the straight line joining the red and blue edges of the two reference spectra.
! Nov. 2016: added correction factors  to normalizse the spectra
! Dec. 2023: Adapted to the FUV and E-IRTF. A.Vazdekis
!
!ANY UPDATE OF THE OUTPUT SPECTRA NEED TO BE UPDATED: FIND WHERE BY LOOKING AT 'UPDATED'
!THIS PROGRAM ASSUMES THAT MILES IS SEGMENT NUMBER 2. ANY CHANGE MUST BE UPDATED!!!!
!
module parameters
implicit none
integer, parameter:: dp=kind(0.d0)
integer, parameter:: nskip_spec=0
real(dp), allocatable:: ws(:), fs(:)
real(dp), allocatable:: s1(:), sf(:) 
INTEGER:: nages !# of SSP ages, defined in findSSPfiles
INTEGER:: nZ    !# of SSP metallicities, defined in findSSPfiles
INTEGER:: nimf    !# of IMF slopes, defined in findSSPfiles
INTEGER:: nw
INTEGER, ALLOCATABLE:: iwar(:)
REAL(kind=dp), ALLOCATABLE, DIMENSION(:):: tssp
REAL(kind=dp), ALLOCATABLE, DIMENSION(:):: Zssp
REAL(kind=dp), ALLOCATABLE, DIMENSION(:):: mussp
REAL(dp), parameter:: res_ngsl=5. ! ngsl resolution just before miles starting point
CHARACTER(LEN=100), PARAMETER :: tmp_path='tmp'
character(len=80),  dimension(61) :: cabem
character(len=100):: suffix
integer, parameter:: flag_smooth=0 ! 0=smooth the continuum in the common window (for miles and cat) / 1=don't smooth
!-------------------------- parameters of the broadening routine --------------
        INTEGER NMAXFFT
        PARAMETER (NMAXFFT=2.**18.)     ! UPDATED 2**16 maximum spectrum size for FFT work
!------------------------------------------------------------------------------
        INTEGER NMAXBUFF                   !maximum no. of Buffers in plotsplus
        PARAMETER (NMAXBUFF=6)
!
        CHARACTER*12 CLAVE_RED         !fixed identification for REDUCEME files
        PARAMETER (CLAVE_RED='abcdefghijkl')
!
        INTEGER NBDMAX        !maximum no. of bands allowed for a generic index
        PARAMETER (NBDMAX=198)                       !note that NBDMAX=NWVMAX/2
!
        INTEGER NWVMAX      !2*maximum no. of bands allowed for a generic index
        PARAMETER (NWVMAX=396)                       !see subroutine selindex.f
!
        INTEGER NINDMAX                         !maximum no. of indices defined
        PARAMETER (NINDMAX=100)                      !see subroutine selindex.f
!
        INTEGER NLINMAX                           !maximum no. of lines defined
        PARAMETER (NLINMAX=40)                       !see subroutine sellines.f
!
        INTEGER NPBANDMAX           !maximum no. of points to define a bandpass
        PARAMETER (NPBANDMAX=40)                      !see subroutine selband.f
!
        INTEGER MAX_ID_RED      !maximum number of simultaneous graphic devices
        PARAMETER (MAX_ID_RED=8)
!
        CHARACTER*12 CREDUCEVERSION                            !current version
        PARAMETER (CREDUCEVERSION='REDUCEMEv4.1')
!------------------------------------------------------------------------------
        INTEGER NSCAN                 !image dimension in the spatial direction
        !INTEGER NCHAN              !image dimension in the wavelength direction
        !REAL(DP)::    STWV                     !central wavelength of the first pixel
        !REAL(DP)::     DISP       !dispersion (Angs/pixel) in the wavelength direction
!------------------------------------------------------------------------------
        REAL(DP)::  AIRMASS                                                   !airmass
        REAL(DP)::  TIMEXPOS                                            !exposure time
        CHARACTER*255 OBJECT                       !name of the observed object
        CHARACTER*255 FITSFILE             !file name of the original FITS file
        CHARACTER*255 COMMENT       !comment to be included in the image header
        CHARACTER*20 THISPROGRAM                !current program being executed
! NOTE: if OBJECT, FITSFILE or COMMENT are redimensioned, the changes must
! also be performed in basicred.f and imath.f
!------------------------------------------------------------------------------
        COMMON/BLKRED01/AIRMASS
        COMMON/BLKRED02/TIMEXPOS
        COMMON/BLKRED03/OBJECT
        COMMON/BLKRED04/FITSFILE
        COMMON/BLKRED05/COMMENT
        COMMON/BLKRED06/THISPROGRAM
!------------------------------------------------------------------------------
        CHARACTER*75 INFILEX                    !function
        CHARACTER*75 OUTFILEX                   !function
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        CHARACTER*255   READC               !read character
!------------------------------------------------------------------------------
        REAL(DP), PARAMETER:: PI=3.141592654D0
        REAL(DP), PARAMETER:: C=299792.458D0
!------------------------------------------------------------------------------

end module parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program join
use parameters
implicit none 
integer:: ispec, nspec, nmodo, k, ns, flag_res, flag_cat, im, ia, imu, ii, i, nnew,  il
integer, allocatable, dimension(:):: n, flag, nskip
real(dp):: res_out, dl_out, l0, l00, nu, wmiles0, scalef1
real(dp), allocatable, dimension(:) :: res_in, disp, cfact
real(dp), allocatable, dimension(:,:) :: W, F, lcont
!real(dp), allocatable:: s1(:), sf(:)
character(len=200),  allocatable, dimension(:) :: specname
character(len=200):: output,  fname, string, outdir, fwar, dir
character(len=1):: s, l
character(len=10):: lmet, lage, limf
!character(len=20):: suffix
logical:: lex, lex1


write(*,*) '-------------------------------------------------------'
write(*,*) 'INSERT INPUT VALUES:'
write(*,*) 'Resolution of the ouptut spectra, flag_res (0=FWHM [AA]; 1=sigma [km/s] )'
read(*,*) res_out
write(*,*) 'Is the resolution to be intended as FWHM [AA] or sigma [km/s]? [0/1]'
read(*,*) flag_res
write(*,*) 'Is segment 4 the CAT spectrum? [yes=1/no=0]'
read(*,*) flag_cat
write(*,*) 'Dispersion of the output spectra [AA]'
read(*,*) dl_out
write(*,*) 'Name of the output spectra'
read(*,*) outdir
write(*,*) 'Number of spectra to be joined'
read(*,*) nspec
write(*,*) 'Suffix to be added to the spectrum name'
!read(*,*) suffix

write(*,*) 'Filename of the input spectrum'
write(*,*) 'Number of lines to be skipped in the input file'
write(*,*) 'Initial resolution of the input spectrum [AA/pixel]'
write(*,*) 'Wavelenghts delimitating the regions where to compute the continuum on the left/right sides:(l1/l2/l3/l4[AA])'
write(*,*) 'Has the spectrum to be considered as reference (no scaling will be applied) or not ? [0/1]'

!NSPEC=2

allocate(specname(nspec))
allocate(nskip(nspec))
allocate(res_in(nspec))
allocate(lcont(4, nspec))
allocate(N(nspec))
allocate(W(nspec, nmaxfft)) !increase !?
allocate(F(nspec, nmaxfft))
allocate(flag(nspec))
allocate(disp(nspec))
allocate(cfact(nspec))

nmodo=5

!file listing  the spectra having NaN
open(unit=12, file='WarningNaN.txt')

!-------------  READ, SMOOTH and RESAMPLE the spectra ----------------------------------------------

do ispec=1, nspec
   read(*,*) specname(ispec), nskip(ispec), res_in(ispec),  disp(ispec), (lcont(k,ispec), k=1,4), flag(ispec), cfact(ispec)
enddo

allocate(iwar(nspec))

!check if the first string is file or dir. In the latter case read all the files in the dir
!and store suffix for met and age
string=specname(1)
inquire(file=string, exist=lex) 

if(lex .eqv. .false.) then !string is a path for a dir
   call findSSPfiles(string, nskip(1))
else !string is a file
   nZ=1
   nages=1
   NIMF=1
   allocate(Zssp(nZ))
   allocate(tssp(nages))
   Zssp=0.
   tssp=0.
endif

write(*,*) 'Nages:',nages
write(*,*) 'Nmet:', nZ
write(*,*) 'Nimf:', nimf

do imu=1, nimf
do im=1, nZ
do ia=1, nages
!for Salpeter imf, solar metallicity, t=10Gyr:
!do imu=4,4
!do im=6,6
!do ia =45,45

write(limf, '(F4.2)') mussp(imu)

write(lmet,'(f4.2)') abs(Zssp(im))
if (Zssp(im) < 0 ) then
   s = 'm'
else
   s = 'p'
endif

ii=index(trim(adjustl(lmet)),'.')
lmet=trim(adjustl(s))//trim(adjustl(lmet))//repeat('0',max(0,2-len_trim(adjustl(lmet(ii+1:)))))

write(lage,'(f7.4)') tssp(ia)
lage = repeat( '0', 7-len_trim(adjustl(lage))) // adjustl(lage)


!First read MILES to compute initial wavel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ispec=2 !UPDATED: THIS MIGHT BE CHANGED IF MILES IS NOT SEGMENT 2 OF SPECTRA TO BE JOINED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(lex .eqv. .false.) then
   !fname=trim(adjustl(specname(ispec)))//'Z'//trim(adjustl(lmet))//'T'//trim(adjustl(lage))//trim(adjustl(suffix))
   fname=trim(adjustl(specname(ispec)))//trim(adjustl(limf))//'Z'//trim(adjustl(lmet))//'T'&
        //trim(adjustl(lage))//trim(adjustl(suffix))
else
   fname=specname(ispec)
endif

call spec(fname, nskip(ispec), ns) ! --> ns, ws, fs
!ref wavel array: ws(1) --> ws(ns)
wmiles0=ws(1)

call header(fname, 61) !UPDATED: THIS COULD BE MODIFIED IF THE HEADER ROWS VARY

iwar(:)=0 !iwar=0 ok; iwar=1 --> NaN in one of the spectral bits

do ispec=1, nspec 

if(lex .eqv. .false.) then
   !fname=trim(adjustl(specname(ispec)))//'Z'//trim(adjustl(lmet))//'T'//trim(adjustl(lage))//'_iPp0.00_baseFe'
   !output=trim(adjustl(outdir))//'Z'//trim(adjustl(lmet))//'T'//trim(adjustl(lage))//'_iPp0.00_baseFe'
   fname=trim(adjustl(specname(ispec)))//trim(adjustl(limf))//'Z'//trim(adjustl(lmet))//'T'//&
        trim(adjustl(lage))//trim(adjustl(suffix))
   output=trim(adjustl(outdir))//trim(adjustl(limf))//'Z'//trim(adjustl(lmet))//'T'//&
        trim(adjustl(lage))//trim(adjustl(suffix))
else
   fname=specname(ispec)
   output=outdir
endif

!get dir
il=index(outdir,'/',back=.true.)
dir=outdir(1:il-1)
inquire(file=dir, exist=lex1) 
if(lex1 .eqv. .false.) call system('mkdir '//trim(adjustl(dir)))

if(ispec>1) then
   if(flag(ispec)==flag(ispec-1)) then
      write(*,*) '-----------------------------------------------------------------------------------------------------------'
      write(*,*) 'Check input file!!'
      write(*,*) 'Two contiguous spectra must have different flags [0/1], one of the two should be used for normalization'
      write(*,*) '-----------------------------------------------------------------------------------------------------------'
      stop
   endif
endif

call spec(fname, nskip(ispec), ns)

!write(*,'(I3,1x, A100, 2(1x,F10.4))') ispec, fname, ws(1), ws(ns)

fs(:)=fs(:)*cfact(ispec) !multiply by correction factor

if(res_in(ispec) < res_out) call broad(ns, res_in(ispec), res_out, flag_res, disp(ispec) ) !input-output (ws, fs) declared in the module

if(ispec==1) l00=wmiles0-int((wmiles0-ws(1))/dl_out)*dl_out !UPDATED: IF MILES IS SEGMENT 1 THEN THIS HAS TO BE UPDATED

l0=ws(1) 
nu=(l0-l00)/dl_out
if( nint(nu)/(nu) /= 1)  l0=l00+nint(nu)*dl_out

nnew=ns
call myrebin(ws, fs, ns, nmodo, dl_out, l0, nnew)

if(allocated(ws) .eqv. .true.) deallocate(ws, fs)
allocate(ws(nnew))
allocate(fs(nnew))

!n=nnew
do i=1,nnew
   ws(i)=s1(i)
   fs(i)=sf(i)
enddo
!endif

n(ispec)=nnew
!W(ispec,1:ns)=ws(1:ns)
!F(ispec,1:ns)=fs(1:ns)
W(ispec,1:nnew)=ws(1:nnew) !new
F(ispec,1:nnew)=fs(1:nnew) !new 

if(iwar(ispec).eq. 0 .and. isnan(f(ispec,1)) .eqv. .true. ) then
   iwar(ispec)=1
   if(iwar(ispec-1) == 0) write(12,'(a80,3x,i1,x,es15.7)') fname, iwar(ispec), f(ispec,1)
endif
enddo

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

write(*,*) '------------------------------------------------------------------------------'
write(*,*) 'Writing output on file: ',trim(adjustl(output))

open(unit=10, file=output)

call join_spec(nspec, n, W, F, lcont, flag_cat, flag, 10,scalef1)
close(10)

enddo !age loop
enddo !met loop
enddo !imf loop

close(12)

end program join

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE header(name, nskip)
USE parameters
!------------------------------------------------------------
!read header from miles file
!input: name, nskip
!output: cabem
!ANY UPDATE OF THE OUTPUT SPECTRA NEED TO BE UPDATED HERE: FIND WHERE BY LOOKING AT 'UPDATED'
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER:: nskip, i, ns, errstat
CHARACTER(LEN=*) :: name
CHARACTER(LEN=10) :: date, time

open(unit=3,file=name,IOSTAT=errstat,action='read')
if (errstat /= 0) then    ! check for errors
   write(*,*) 'Internal error'
   write(*,*) 'Could not open the file containing the spectrum', trim(adjustl(name))
   stop
end if

do i=1,nskip
   read(3,'(A80)') cabem(i)
enddo

!cabem(3)='NAXIS1  =              53688  /  Length of axis'
!cabem(9)='OBJECT  = A_1680.20-49999.40     /'
cabem(3)='NAXIS1  =              54266  /  Length of axis'
cabem(9)='OBJECT  = A_1160.00-49999.40     /'
cabem(10)='FILENAME= E.FITS            /  IRAF filename'
call date_and_time(date, time)

cabem(14)='DATE    = 20/12/23          / Date FITS file'
cabem(15)='IRAF-TLM= '//trim(adjustl(time(1:2)))//':'//trim(adjustl(time(3:4)))//':'&
      //trim(adjustl(time(5:6)))//' ('//trim(adjustl(date(7:8)))//'/'&
      //trim(adjustl(date(5:6)))//'/'//trim(adjustl(date(3:4)))//') / Last modification'
!cabem(18)='CRVAL1  =             1680.20 / central wavel pixel 1'
cabem(18)='CRVAL1  =             1160.00 / central wavel pixel 1'

close(3)

END SUBROUTINE HEADER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE spec(name, nskip, ns)
USE parameters
!------------------------------------------------------------
!read file with spectrum
!input: name, nskip
!output: ns, ws, fs
!ANY UPDATE OF THE OUTPUT SPECTRA NEED TO BE UPDATED HERE: FIND WHERE BY LOOKING AT 'UPDATED'
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER:: nskip, i, ns, errstat
CHARACTER(LEN=*) :: name
CHARACTER(LEN=100),ALLOCATABLE,DIMENSION(:) :: sname

open(unit=3,file=name,IOSTAT=errstat,action='read')
if (errstat /= 0) then    ! check for errors
   write(*,*) 'Internal error'
   write(*,*) 'Could not open the file containing the spectrum', trim(adjustl(name))
   stop
end if

call skip(3,nskip)
call crow(3,ns)

if(allocated(ws) .eqv. .true.) deallocate(ws,fs)
allocate(ws(ns))
allocate(fs(ns))

call skip(3,nskip)
do i=1, ns
   read(3,*,iostat=errstat) ws(i), fs(i)
   if (errstat /= 0) then    ! check for errors
      write(*,*) 'Check spectrum, error while reading'
      stop
   end if
enddo

300 format(f16.4,10x,es22.16)  ! UPDATED FROM f16.8 TO f16.4,TO ROUND WAVELENGTHS

close(3)

END SUBROUTINE spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE crow(unit,row)
!  count number of rows in the input file

IMPLICIT NONE
INTEGER, INTENT(IN) :: unit
INTEGER, INTENT(OUT) :: row
INTEGER :: ierr,i
CHARACTER(len=1) :: x

  row=0 
  count:DO
     READ(unit,*,IOSTAT=ierr) x
     IF(ierr/=0) EXIT
     row= row+1
  ENDDO count
  REWIND (unit)
  
END SUBROUTINE crow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE skip(unit,nskip)
!skip nskip lines from file unit
!unit must be already open
IMPLICIT NONE
INTEGER:: unit, nskip, i, errstat
CHARACTER(LEN=1) :: x !line to skip

do i=1, nskip
   read(unit,*, IOSTAT=errstat) x
   if (errstat /= 0) then    ! check for errors
   write(*,*) 'Internal error'
   write(*,*) 'Check header of file unit',unit
   stop
end if
enddo

END SUBROUTINE skip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE cskip(unit, nskip)
IMPLICIT NONE
INTEGER:: unit, nskip
CHARACTER(LEN=1) :: x

nskip=0
do
   read(unit,*) x
   if(x/='#') exit
   nskip=nskip+1
enddo
backspace(unit) 

END SUBROUTINE cskip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine broad(ns, res_in, res_out, flag_res, disp) !generalize for input sigma in km/s 
use parameters
implicit none
integer:: ns, j, nchan, flag_res
integer, parameter:: ncmax=20000
integer, parameter:: nsmax=500
real(dp) :: DISP, STWV, res_in, res_out
real(dp):: A(NCMAX,NSMAX)
real(DP):: S(NCMAX),SS(NCMAX),X(NCMAX)
real(dp):: SIGMA,SIGMAVEC(NCMAX)
real(dp):: WAV(NCMAX)


sigma=dsqrt(res_out**2.d0-res_in**2.d0) 
if(sigma==0.) stop

nchan=ns
stwv=ws(1)

do j=1, ns
wav(j)=ws(j)
s(j)=fs(j)
if(flag_res==0) then
   sigmavec(j)=sigma*c/(2.35d0*wav(j))
else if(flag_res==1) then
   sigmavec(j)=sigma
else
  write(*,*) 'This option is not handled yet: flag_res should be 0 or 1'
endif

enddo

call broaden_elena(s,ss,nchan,stwv,disp,sigmavec,.false.)

do j=1, ns
ws(j)=wav(j)
fs(j)=ss(j)
enddo

end subroutine broad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE BROADEN(SS1,SS2,NCHAN,STWV,DISP,SIGMA,LERR)
!
        use parameters
        IMPLICIT NONE
        INTEGER NCHAN
        REAL(DP):: SS1(NCHAN),SS2(NCHAN)
        REAL(DP):: STWV,DISP
        REAL(DP):: SIGMA(NCHAN)
        LOGICAL LERR
!
        INTEGER I,K,K1,K2,N,INC
        REAL(DP):: W0,W,W1,W2
        REAL(DP):: SUM,FACTOR1,FACTOR2
        REAL(DP):: PIXSIGMA
        REAL(DP):: FINTGAUSS
        REAL(DP):: FINTGAUSSE
        REAL(DP):: FACTOR01,FACTOR02,PIXSIGMA0
!------------------------------------------------------------------------------
! calculamos todas las constantes fuera de los bucles para aumentar velocidad
!------------------------------------------------------------------------------

        IF(LERR)THEN
          DO I=1,NCHAN
            IF(SIGMA(I).GT.0.0)THEN
              FACTOR01=-C*C/(2.d0*SIGMA(I)*SIGMA(I))
              FACTOR02=C/(SQRT(2.d0*PI)*SIGMA(I))
              PIXSIGMA0=SIGMA(I)/(C*DISP)
              W0=REAL(I-1)*DISP+STWV                          !longitud de onda
              SUM=0.0d0
              FACTOR1=FACTOR01/(W0*W0)
              FACTOR2=FACTOR02/W0
              PIXSIGMA=PIXSIGMA0*W0     !numero de pixels que equivalen a sigma
              INC=NINT(6.*PIXSIGMA)               !numero de pixels a cada lado
              K1=I-INC                               !limite inferior: -6 sigma
              IF(K1.LT.1) K1=1                       !ojo con el borde inferior
              K2=I+INC                               !limite superior: +6 sigma
              IF(K2.GT.NCHAN) K2=NCHAN               !ojo con el borde superior
              N=NINT(10./PIXSIGMA) !exigimos 10 intervalos para muestrear sigma
              IF(MOD(N,2).NE.0) N=N+1                 !N debe ser un numero par
              IF(N.LT.10) N=10          !como minimo usamos 10 intervalos/pixel
              DO K=K1,K2  !sumamos entre +-6 sigma (salvo en los bordes, claro)
                W=REAL(K-1)*DISP+STWV      !l.d.o. del centro del pixel K-esimo
                W1=W-DISP/2.d0                 !l.d.o. inferior del pixel K-esimo
                W2=W+DISP/2.d0                 !l.d.o. superior del pixel K-esimo
                SUM=SUM+SS1(K)*SS1(K)*FINTGAUSSE(W1,W2,N,W0,FACTOR1)
              END DO
              SUM=SQRT(SUM)
              SS2(I)=FACTOR2*SUM
            ELSE
              SS2(I)=SS1(I)
            END IF
          END DO
        ELSE
          DO I=1,NCHAN
            IF(SIGMA(I).GT.0.0)THEN
              FACTOR01=-C*C/(2.d0*SIGMA(I)*SIGMA(I))
              FACTOR02=C/(SQRT(2.*PI)*SIGMA(I))
              PIXSIGMA0=SIGMA(I)/(C*DISP)
              W0=REAL(I-1)*DISP+STWV                          !longitud de onda
              SUM=0.0
              FACTOR1=FACTOR01/(W0*W0)
              FACTOR2=FACTOR02/W0
              PIXSIGMA=PIXSIGMA0*W0     !numero de pixels que equivalen a sigma
              INC=NINT(6.d0*PIXSIGMA)               !numero de pixels a cada lado
              K1=I-INC                               !limite inferior: -6 sigma
              IF(K1.LT.1) K1=1                       !ojo con el borde inferior
              K2=I+INC                               !limite superior: +6 sigma
              IF(K2.GT.NCHAN) K2=NCHAN               !ojo con el borde superior
              N=NINT(10./PIXSIGMA) !exigimos 10 intervalos para muestrear sigma
              IF(MOD(N,2).NE.0) N=N+1                 !N debe ser un numero par
              IF(N.LT.10) N=10          !como minimo usamos 10 intervalos/pixel
              DO K=K1,K2  !sumamos entre +-6 sigma (salvo en los bordes, claro)
                W=REAL(K-1)*DISP+STWV      !l.d.o. del centro del pixel K-esimo
                W1=W-DISP/2.d0                 !l.d.o. inferior del pixel K-esimo
                W2=W+DISP/2.d0                 !l.d.o. superior del pixel K-esimo
                SUM=SUM+SS1(K)*FINTGAUSS(W1,W2,N,W0,FACTOR1)
              END DO
              SS2(I)=FACTOR2*SUM
            ELSE
              SS2(I)=SS1(I)
            END IF
          END DO
        END IF
!
      END SUBROUTINE BROADEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Broadening in the first pixels adjusted by using specular pixels for pixel < w0 (or >wf)
        SUBROUTINE BROADEN_ELENA(SS1,SS2,NCHAN,STWV,DISP,SIGMA,LERR)
!
        use parameters
        IMPLICIT NONE
        INTEGER NCHAN
        REAL(DP):: SS1(NCHAN),SS2(NCHAN)
        REAL(DP):: STWV,DISP
        REAL(DP):: SIGMA(NCHAN)
        LOGICAL LERR
!
        INTEGER I,K,K1,K2,N,INC, KK
        REAL(DP):: W0,W,W1,W2
        REAL(DP):: SUM,FACTOR1,FACTOR2
        REAL(DP):: PIXSIGMA
        REAL(DP):: FINTGAUSS
        REAL(DP):: FINTGAUSSE
        REAL(DP):: FACTOR01,FACTOR02,PIXSIGMA0
!------------------------------------------------------------------------------
! calculamos todas las constantes fuera de los bucles para aumentar velocidad
!------------------------------------------------------------------------------

        IF(LERR)THEN
          DO I=1,NCHAN
            IF(SIGMA(I).GT.0.0)THEN
              FACTOR01=-C*C/(2.d0*SIGMA(I)*SIGMA(I))
              FACTOR02=C/(SQRT(2.d0*PI)*SIGMA(I))
              PIXSIGMA0=SIGMA(I)/(C*DISP)
              W0=REAL(I-1)*DISP+STWV                          !longitud de onda
              SUM=0.0d0
              FACTOR1=FACTOR01/(W0*W0)
              FACTOR2=FACTOR02/W0
              PIXSIGMA=PIXSIGMA0*W0     !numero de pixels que equivalen a sigma
              INC=NINT(6.*PIXSIGMA)               !numero de pixels a cada lado
              K1=I-INC                               !limite inferior: -6 sigma
              IF(K1.LT.1) K1=1                       !ojo con el borde inferior
              K2=I+INC                               !limite superior: +6 sigma
              IF(K2.GT.NCHAN) K2=NCHAN               !ojo con el borde superior
              N=NINT(10./PIXSIGMA) !exigimos 10 intervalos para muestrear sigma
              IF(MOD(N,2).NE.0) N=N+1                 !N debe ser un numero par
              IF(N.LT.10) N=10          !como minimo usamos 10 intervalos/pixel
              DO K=K1,K2  !sumamos entre +-6 sigma (salvo en los bordes, claro)
                W=REAL(K-1)*DISP+STWV      !l.d.o. del centro del pixel K-esimo
                W1=W-DISP/2.d0                 !l.d.o. inferior del pixel K-esimo
                W2=W+DISP/2.d0                 !l.d.o. superior del pixel K-esimo
                SUM=SUM+SS1(K)*SS1(K)*FINTGAUSSE(W1,W2,N,W0,FACTOR1)
              END DO
              SUM=SQRT(SUM)
              SS2(I)=FACTOR2*SUM
            ELSE
              SS2(I)=SS1(I)
            END IF
          END DO
        ELSE
          DO I=1,NCHAN
            IF(SIGMA(I).GT.0.0)THEN
              FACTOR01=-C*C/(2.d0*SIGMA(I)*SIGMA(I))
              FACTOR02=C/(SQRT(2.*PI)*SIGMA(I))
              PIXSIGMA0=SIGMA(I)/(C*DISP)
              W0=REAL(I-1)*DISP+STWV                          !longitud de onda
              SUM=0.0
              FACTOR1=FACTOR01/(W0*W0)
              FACTOR2=FACTOR02/W0
              PIXSIGMA=PIXSIGMA0*W0     !numero de pixels que equivalen a sigma
              INC=NINT(6.d0*PIXSIGMA)               !numero de pixels a cada lado
              K1=I-INC                               !limite inferior: -6 sigma
              K2=I+INC                               !limite superior: +6 sigma
              !! modified by elena
              !IF(K1.LT.1) K1=1                       !ojo con el borde inferior
              !IF(K2.GT.NCHAN) K2=NCHAN               !ojo con el borde superior
              !! end OF modification
              N=NINT(10./PIXSIGMA) !exigimos 10 intervalos para muestrear sigma
              IF(MOD(N,2).NE.0) N=N+1                 !N debe ser un numero par
              IF(N.LT.10) N=10          !como minimo usamos 10 intervalos/pixel
              DO K=K1,K2  !sumamos entre +-6 sigma (salvo en los bordes, claro)
                !! modified by elena
                KK=K
                IF(K.LT.1) KK=-(K-I)+I !for pixels below the first one (K<1) I use the flux in the specular pixel
                IF(K.GT.NCHAN) KK=NCHAN-(K-NCHAN)
                W=REAL(K-1)*DISP+STWV      !l.d.o. del centro del pixel K-esimo
                W1=W-DISP/2.d0                 !l.d.o. inferior del pixel K-esimo
                W2=W+DISP/2.d0                 !l.d.o. superior del pixel K-esimo
                SUM=SUM+SS1(KK)*FINTGAUSS(W1,W2,N,W0,FACTOR1)
                !! end modification
              END DO
              SS2(I)=FACTOR2*SUM
            ELSE
              SS2(I)=SS1(I)
            END IF
          END DO
        END IF
!
      END SUBROUTINE BROADEN_ELENA
! Integral de una gaussiana de anchura SIGMA por el metodo de Simpson, entre 
! X1 y X2, con N intervalos (N debe ser par). Utilizamos la formula (1prim) del
! libro "Calculo Numerico Fundamental", Demidovich y Maron, pag. 656 (ojo, hay
! una errata en la formula de sigma2).
        FUNCTION FINTGAUSS(X1,X2,N,X0,FACTOR1)
        USE parameters
        IMPLICIT NONE
        REAL(DP)::FINTGAUSS
        REAL(DP):: X1,X2
        INTEGER N
        REAL(DP):: X0,FACTOR1
!
        REAL(DP):: DX1,DX2
        REAL(DP):: DX0,DFACTOR1
!
        INTEGER I
        REAL(DP)::  X,Y1,Y2
        REAL(DP):: H
        REAL(DP):: SUM1,SUM2
!------------------------------------------------------------------------------
! chequeamos que N es par
        IF(MOD(N,2).NE.0)THEN
          STOP 'FATAL ERROR: N is odd in subroutine FINTGAUSS!'
        END IF
! pasamos a doble precision las variables REAL de entrada
        DX1=DBLE(X1)
        DX2=DBLE(X2)
        DX0=DBLE(X0)
        DFACTOR1=DBLE(FACTOR1)
!
        H=(DX2-DX1)/DBLE(N) !tamaño de cada intervalo
!
        SUM1=0.D0
        DO I=1,N-1,2
          X=DX1+H*DBLE(I)
          SUM1=SUM1+DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        SUM2=0.D0          !sumamos terminos pares salvo el primero y el ultimo
        DO I=2,N-2,2
          X=DX1+H*DBLE(I)
          SUM2=SUM2+DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        Y1=DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))       !funcion en el primer punto
        Y2=DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))       !funcion en el ultimo punto
!
        FINTGAUSS=REAL((Y1+Y2+4.D0*SUM1+2.D0*SUM2)*H/3.D0)            !solucion
!
      END FUNCTION FINTGAUSS
!
!******************************************************************************
! Integral de una gaussiana**2 de anchura SIGMA por el metodo de Simpson, entre 
! X1 y X2, con N intervalos (N debe ser par). Utilizamos la formula (1prim) del
! libro "Calculo Numerico Fundamental", Demidovich y Maron, pag. 656 (ojo, hay
! una errata en la formula de sigma2).
        FUNCTION FINTGAUSSE(X1,X2,N,X0,FACTOR1)
        USE parameters  
        IMPLICIT NONE
        REAL(DP):: FINTGAUSSE
        REAL(DP):: X1,X2
        INTEGER N
        REAL(DP):: X0,FACTOR1
!
        REAL(DP):: DX1,DX2
        REAL(DP):: DX0,DFACTOR1
!
        INTEGER I
        REAL(DP):: X,Y1,Y2
        REAL(DP):: SUM1,SUM2
!------------------------------------------------------------------------------
! chequeamos que N es par
        IF(MOD(N,2).NE.0)THEN
          STOP 'FATAL ERROR: N is odd in subroutine FINTGAUSSE!'
        END IF
! pasamos a doble precision las variables REAL de entrada
        DX1=DBLE(X1)
        DX2=DBLE(X2)
        DX0=DBLE(X0)
        DFACTOR1=DBLE(FACTOR1)
!
        SUM1=0.D0                                     !sumamos terminos impares
        DO I=1,N-1,2
          X=DX1+(DX2-DX1)*DBLE(I)/DBLE(N)
          SUM1=SUM1+ &
          DEXP(DFACTOR1*(X-DX0)*(X-DX0))*DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        SUM2=0.D0          !sumamos terminos pares salvo el primero y el ultimo
        DO I=2,N-2,2
          X=DX1+(DX2-DX1)*DBLE(I)/DBLE(N)
          SUM2=SUM2+ &
          DEXP(DFACTOR1*(X-DX0)*(X-DX0))*DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        Y1=DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))* &
        DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))
        Y2=DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))* &
        DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))
!
        FINTGAUSSE= &
        DBLE((Y1+Y2+4.D0*SUM1+2.D0*SUM2)*(DX2-DX1)/DBLE(3*N))        !solucion
!
   END FUNCTION FINTGAUSSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!SUBROUTINE myrebin(lam, flux, n, nmodo, dl, l0, nnew, s1, sf)
SUBROUTINE myrebin(lam, flux, n, nmodo, dl, l0, nnew)
!adapted from Vazdekis' program lector.f
!change lambda scale from log/ln to linear and vice versa by interpolating (polynomial)
!flux is not conserved in the rebinning
USE parameters
IMPLICIT NONE
INTEGER:: n, nmodo, i, k1, err, nnew, itonto
REAL(KIND=dp), DIMENSION(n) :: lam, flux!, s11, s22, s1, sf
REAL(KIND=dp), DIMENSION(n) :: s11, s22
!real(dp), allocatable, dimension(:) :: s1, sf
REAL(KIND=dp)::disps, dlamb,  dl, l0, dy

!write(*,*) 'Starting resampling using polynomial interpolation...'
!write(*,*) ''

nnew=n
if(nmodo==1) then
   lam(:)=dlog10(lam(:))
else if(nmodo==2) then
   lam(:)=dlog(lam(:))
else if(nmodo==3) then
   lam(:)=10.d0**(lam(:))
else if(nmodo==4) then
   lam(:)=dexp(lam(:))
else if(nmodo==5) then !linear binning with different disp
   nnew=1+int((lam(n)-lam(1))/dl)
else
   write(*,*) 'This option is not handled yet'
   write(*,*) 'nmodo should be 1/2/3/4 '
endif

do i=1, n
   s11(i)=DBLE(lam(i))
   s22(i)=DBLE(flux(i))
enddo

if(nmodo==5) then
   disps=dl
   lam(1)=l0
else
   disps=dble(dabs(lam(n)-lam(1))/(n-1))
endif

dlamb=dble(lam(1))-disps

if(allocated(s1) .eqv. .true.) deallocate(s1, sf)
allocate(s1(nnew))
allocate(sf(nnew))

do i=1,nnew 

dlamb=dble(dlamb+disps)
s1(i)=dlamb
call hunt(s11, n, s1(i), k1)
!write(*,*) 'nnew=1',s11(k1), s22(k1), s1(i)

if(k1==0) k1=1

if(i==1) then
   call polint(s11(k1), s22(k1), 2, s1(i), sf(i), dy)
   itonto=2
   call interp_sgl(s11, s22, 1, s1(i), sf(i), dy)
else if(i==2 .or. k1 == 1) then
   call polint(s11(k1), s22(k1), 3, s1(i), sf(i), dy)
   itonto=3
else if(i == (nnew-1) ) then
   call polint(s11(k1-1),s22(k1-1),3,s1(i),sf(i),dy)
   itonto=3
else if(i == nnew) then
   call polint(s11(k1-1),s22(k1-1),2,s1(i),sf(i),dy)
   itonto=2
else
  call polint(s11(k1-1),s22(k1-1),4,s1(i),sf(i),dy)
   itonto=4
endif
!write(*,*) 'NEWREBIN1',s11(k1), s22(k1), s1(i),sf(i),i,n,nnew,itonto

enddo

END SUBROUTINE myrebin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE hunt(xx,n,x,jlo)
      USE parameters
      IMPLICIT NONE
! subroutine from the Numerical Recipes
      INTEGER jlo,n
      REAL(KIND=dp)::x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
    END subroutine hunt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


      SUBROUTINE polint(xa,ya,n,x,y,dy)
      USE parameters
      IMPLICIT NONE
!    Numerical Recipes subroutine
      INTEGER n,NMAX
      REAL(KIND=dp)::  dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL(KIND=dp)::   den,dif,dift,ho,hp,w,cc(NMAX),d(NMAX)

      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        cc(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=cc(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          cc(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=cc(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
    END SUBROUTINE polint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE interp_sgl ( x, y, npts, x0, y0, error )
!
!  Purpose:
!    To linearly interpolate the value y0 at position x0, given a
!    set of (x,y) measurements organized in increasing order of x.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/17/96    S. J. Chapman        Original code
!
USE parameters
IMPLICIT NONE

! Declare local parameters
!INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=4) ! Precision

! Declare calling arguments:
INTEGER, INTENT(IN) :: npts           ! Dimension of arrays x and y
REAL(KIND=dp),DIMENSION(npts), INTENT(IN) :: x
                                      ! Independent variable x.
REAL(KIND=dp),DIMENSION(npts), INTENT(IN) :: y
                                      ! Dependent variable y.
REAL(KIND=dp),INTENT(IN) :: x0      ! Point to interpolate.
REAL(KIND=dp),INTENT(OUT) :: y0     ! Interpolated value.
REAL(KIND=DP), INTENT(OUT) :: error         ! Error flag: 0 -- No error
                                      !            -1 -- x0 < x(1)
                                      !             1 -- x0 > x(npts)
! Declare local variables:
INTEGER :: i                     ! Index variable
INTEGER :: ibase                 ! Index for interpolation.
REAL(KIND=dp) :: slope                    ! Slope for interpolation.


! Assume that the input data set is in ascending order of x.
! See if the measurement position x0 is smaller or larger
! than any value in the data set.

IF( npts == 1 ) THEN
   y0=y(1)
   error=0
   return
ENDIF

IF ( x0 < x(1) ) THEN
   error = -1
   y0=y(1)+((x0-x(1))/(x(2)-x(1)))*(y(2)-y(1))
ELSE IF ( x0 > x(npts) ) THEN
   error = 1
   y0=y(npts-1)+((x0-x(npts-1))/(x(npts)-x(npts-1)))*(y(npts)-y(npts-1))

ELSE

   ! Point is between x(1) and x(npts).  Find the two points
   ! that it is between.
    DO i = 1, npts-1
       IF ( (x0 >= x(i)) .AND. (x0 <= x(i+1)) ) THEN
            ibase = i        ! Found the point
            EXIT
       END IF
    END DO

    ! Now linearly interpolate point.
    slope = ( y(ibase+1)-y(ibase) ) / ( x(ibase+1)-X(ibase) )
    y0 = slope * ( x0 - x(ibase) ) + y(ibase)

    ! Set error flag to 0.
    error = 0
END IF

END SUBROUTINE interp_sgl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine join_spec(nspec, n, w, f, lcont, flag_cat, flag, unit, scalef1)
use parameters
!ANY UPDATE OF THE OUTPUT SPECTRA NEED TO BE UPDATED HERE: FIND WHERE BY LOOKING AT 'UPDATED'
integer:: nspec, flag_cat, unit, nout, ispec, j, i1, i2, i3, i4
integer, dimension(nspec):: n, flag, i0, if
real(dp) :: a0, b0, af, bf, ff, scalef, cont_nm, cont_ci, scalef1
real(dp), dimension(nspec):: cont1, cont2,  wi, wf
real(dp), dimension(nspec, nmaxfft):: w, f
real(dp), dimension(4, nspec):: lcont
real(dp), allocatable, dimension(:) :: wout, fout

cont1(:)=0.d0
cont2(:)=0.d0
write(*,*) 'Segments to join [0==>(NGSL+MILES+IRTF], [1==>(NGSL+MILES+MIUSC+CAT+IRTF]:',flag_cat
199 format(I3,4(X,F11.4),2(X,I6),2(X,es22.16))
!----------------      compute the continuum in the overlapping regions  -------------------------------------------
DO ispec=1, nspec

!NaN fluxes --> 0.0

do j=1, n(ispec)
   if(isnan(f(ispec,j)) .eqv. .true.) f(ispec,j)=0.
enddo

!for MILES: first smooth miles SSP to 5AA and then compute the continuum
! only for NGLS-MILES joining !!!!

IF(ispec .eq. 2 ) THEN !miles

   if(allocated(ws) .eqv. .true.) deallocate(ws, fs)
   ns=n(ispec)
   allocate(ws(ns))
   allocate(fs(ns))
   ws(1:ns)=w(ispec,1:ns)
   fs(1:ns)=f(ispec,1:ns)

!   if(flag_smooth .eq. 0) call broad(ns, 2.51d0, 5.d0, 0, 0.9d0 ) !UPDATED input-output (ws, fs) declared      if(flag_smooth .eq. 0) call broad(ns, 2.76d0, 5.d0, 0, 0.9d0 ) !UPDATED input-output (ws, fs) declared in the module

   if(lcont(1,ispec) /= 0. .and. lcont(2,ispec) /= 0.) &
        call continuum(lcont(1,ispec), lcont(2,ispec), ws(:), fs(:), ns, cont1(ispec)) !left
    cont_nm=cont1(ispec)
    
ENDIF

IF(flag_cat .eq. 1 .and. ispec .eq. 4) THEN !ONLY IN CASE SEGMENT 4 IS CAT SPECTRUM   
!              CaT-IRTF: smooth CaT to IRTF resolution for continuum computation

   if(allocated(ws) .eqv. .true.) deallocate(ws, fs)
   ns=n(ispec)
   allocate(ws(ns))
   allocate(fs(ns))
   ws(1:ns)=w(ispec,1:ns)
   fs(1:ns)=f(ispec,1:ns)
 
    if(flag_smooth .eq. 0) call broad(ns, 1.5d0, 4.2d0, 0, 0.9d0 ) !input-output (ws, fs) declared in the module
!    if(flag_smooth .eq. 0) call broad(ns, 1.5d0, 3.47d0, 0, 0.9d0 ) !UPDATED (Joining after MILES)input-output (ws, fs) declared in the module

   if(lcont(3,ispec) /= 0. .and. lcont(4,ispec) /= 0.) &
        call continuum(lcont(3,ispec), lcont(4,ispec), ws(:), fs(:), ns, cont2(ispec)) !left
    cont_ci=cont2(ispec) !cont cat irtf

ENDIF

IF(flag_cat .eq. 0 .and. ispec .eq. 2) THEN !ONLY IN CASE SEGMENT 4 IS CAT SPECTRUM   
  !miles to be smoothed to match irtf

   if(allocated(ws) .eqv. .true.) deallocate(ws, fs)
   ns=n(ispec)
   allocate(ws(ns))
   allocate(fs(ns))
   ws(1:ns)=w(ispec,1:ns)
   fs(1:ns)=f(ispec,1:ns)

   if(flag_smooth .eq. 0) call broad(ns, 2.76d0, 3.47d0, 0, 0.9d0 ) !UPDATED input-output (ws, fs) declared in the module

   if(lcont(3,ispec) /= 0. .and. lcont(4,ispec) /= 0.) &
        call continuum(lcont(3,ispec), lcont(4,ispec), ws(:), fs(:), ns, cont2(ispec)) !left
    cont_mi=cont2(ispec) !cont miles irtf

ENDIF 

if(lcont(1,ispec) /= 0. .and. lcont(2,ispec) /= 0.) &
      call continuum(lcont(1,ispec), lcont(2,ispec), w(ispec,:), f(ispec,:), n(ispec), cont1(ispec)) !left

if(lcont(3,ispec) /= 0. .and. lcont(4,ispec) /= 0.) &
      call continuum(lcont(3,ispec), lcont(4,ispec), w(ispec,:), f(ispec,:), n(ispec), cont2(ispec)) !right

!use broadened continuum
if(ispec == 2) cont1(ispec)=cont_nm !use broadened continuum
if(flag_cat .eq. 1 .and. ispec == 4)  cont2(ispec)=cont_ci !ONLY IN CASE SEGMENT 4 IS CAT SPECTRUM
if(flag_cat .eq. 0 .and. ispec == 2)  cont2(ispec)=cont_mi !ONLY IN CASE SEGMENT 4 IS CAT SPECTRUM   
! if(ispec == 4) cont2(ispec)=cont_ci
!else
! if(ispec == 2) cont2(ispec)=cont_mi
!endif

call igrid(w(ispec,:), n(ispec), lcont(1,ispec), i1)
call igrid(w(ispec,:), n(ispec), lcont(2,ispec), i2)
call igrid(w(ispec,:), n(ispec), lcont(3,ispec), i3)
call igrid(w(ispec,:), n(ispec), lcont(4,ispec), i4)

wi(ispec)=0.5d0*(lcont(1,ispec)+lcont(2,ispec))

if(ispec == 2) then !copia miles dal primo pixel
   wi(ispec)=lcont(1,ispec)!+0.9d0
endif

wf(ispec)=0.5d0*(lcont(3,ispec)+lcont(4,ispec))
if(ispec == 1) then
   wf(ispec)=lcont(3,ispec)
endif

if(flag_cat .eq. 1 .and. ispec==4) wf(ispec)=lcont(4,ispec)
if(flag_cat .eq. 1 .and. ispec==5) wi(ispec)=lcont(2,ispec)

if(flag_cat .eq. 0 .and. ispec==2) wf(ispec)=lcont(4,ispec)
if(flag_cat .eq. 0 .and. ispec==3) wi(ispec)=lcont(2,ispec)

call igrid(w(ispec,:), n(ispec), wi(ispec), i0(ispec))
call igrid(w(ispec,:), n(ispec), wf(ispec), if(ispec))

if(ispec == 1 ) i0(ispec)=1 
if(ispec == nspec) if(ispec)=n(ispec)


write(*,199)ispec,lcont(1,ispec),lcont(2,ispec),lcont(3,ispec),lcont(4,ispec),ns,n(ispec),cont1(ispec),cont2(ispec)
ENDDO
write(*,*)''

!-------------------------------------------------------------------
! Write output file 
!---------------------------------------------------------------------
scalef1=0

!write header
do j=1,61
   write(unit,'(A80)') cabem(j)
enddo
200 format(f16.4,3x, es22.16) ! UPDATED from f16.8

!-----------    join the spectra   -------------------
do ispec=1, nspec


if(flag(ispec) == 0) then !ref spectrum

   do j=i0(ispec), if(ispec)
      write(unit,200) w(ispec,j), f(ispec,j)
   enddo

else if(flag(ispec) == 1 .and. ispec==1) then ! edges --> just scaling

   scalef=cont1(ispec+1)/cont2(ispec)
   scalef1=scalef
   if(f(ispec,j) == 0.) scalef=0.
   do j=i0(ispec), if(ispec)-1
      write(unit,200) w(ispec,j), scalef*f(ispec,j)
   enddo

else if(flag(ispec) ==1 .and. ispec==nspec) then ! edges --> just scaling

   scalef=cont2(ispec-1)/cont1(ispec)
   if(f(ispec,j) == 0.) scalef=0.
   do j=i0(ispec)+1, if(ispec)
      write(unit,200) w(ispec,j), scalef*f(ispec,j)
   enddo

else ! middle --> re-normalize 

   a0=(cont2(ispec)-cont1(ispec))/(wf(ispec)-wi(ispec)) !slope of the initial spectrum
   b0=cont1(ispec)-a0*wi(ispec)
   
   af=(cont1(ispec+1)-cont2(ispec-1))/(wf(ispec)-wi(ispec)) !forcing the slope to be equal to the line joining the ref spectra
   bf=cont2(ispec-1)-af*wi(ispec)

   do j=i0(ispec)+1, if(ispec)-1

      ff=af*w(ispec,j)+bf
      f0=a0*w(ispec,j)+b0
      scalef=ff/f0
      if(f(ispec,j) == 0.) scalef=0.
      write(unit,200) w(ispec,j), f(ispec,j)*scalef 

   enddo
endif


ENDDO

end subroutine join_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE continuum(l1, l2, w, f, n, cont)
USE parameters
IMPLICIT NONE
integer:: n, i1, i2, i
real(dp):: l1, l2, cont, summ
real(dp), dimension(n):: w, f


call igrid(w, n, l1, i1)
call igrid(w, n, l2, i2)
cont=sum(f(i1:i2))/(i2-i1+1) !average counts
!cont=sum(f(i1:i2))


END SUBROUTINE continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE igrid(arr, n, x0, iout)
USE parameters
IMPLICIT NONE
INTEGER:: n, iout, i
REAL(KIND=dp) :: x0
REAL(KIND=dp),DIMENSION(n):: arr, da

do i=1, n
da(i)=dabs(arr(i)-x0)
enddo

iout=minloc(da, 1) !check if arr(iout) needs to be smaller than x0
!if(arr(iout) > x0 ) iout=iout+1


END SUBROUTINE igrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE findSSPfiles(dir, nskip)
USE parameters
IMPLICIT NONE
INTEGER:: i, j, k, nrow, ia, im, l, errstat,ii, row, seed, nskip
REAL(kind=dp):: m, a, harvest, mmu
REAL(kind=dp),DIMENSION(5000):: age,met,mu
CHARACTER(LEN=*) ::dir
CHARACTER(LEN=200):: filename, amet, aage, rl, tmpfile, aimf
CHARACTER(LEN=1):: s, skip
LOGICAL:: lex_tmp

!create temporary file with random name
!in order not to overwrite it in case of simultaneous run
seed= 589666
call random_seed
call random_number(harvest)

write(rl,'(i7)') INT(harvest*1.d7)
tmpfile=trim(adjustl(tmp_path))//trim(adjustl(rl))

inquire(file=tmpfile, exist=lex_tmp) !check if tmfile with same name exists
if(lex_tmp .eqv..true.) then
   write(*,*) 'tmp file already in use ',trim(adjustl(tmpfile)) 
   write(*,*) 'Try again...'
   stop
endif


dir=dir(1:len_trim(dir)-4)

!listing files in dir

call system('find '//trim(adjustl(dir))//' -name \*.* >'//trim(adjustl(tmpfile)))
!find ../MODELS_ALPHA/OUT_MILES_UN_iPp0.00_baseFe/ -name \* > tmp

open(unit=1,file=tmpfile, iostat=errstat, action='read')
if (errstat /= 0) then    ! check for errors
   write(*,*) 'Internal error'
   write(*,*) 'Could not open the tmp file', tmpfile
   stop
end if

call crow(1,nrow)

l=len(trim(dir))+8 !old format
im=0
ia=0
age=0
met=0
mu=0

nages=0
nZ=0
nimf=0

list: do i=1,nrow

read(1,'(a200)') filename

aimf=filename(l-3:l) !get substring with imf slope
s=filename(l+2:l+2) !get substring with sign of metallicity
amet=filename(l+3:l+6) !get substring with metallicity value
aage=filename(l+8:l+14)!get substr with age value

if(i == 1) suffix=filename(l+15:)

read(aimf,'(F5.2)') mmu
read(amet,*) m !convert into real
read(aage,'(F7.4)') a

!count # of imf slopes
do j=1, nimf
   if(mmu == mu(j)) goto 200
enddo

nimf=nimf+1
mu(nimf)=mmu
200 continue

if(s=='m') m=-m

!count # of metallicity values
do j=1,nZ
   if(m == met(j)) goto 100
enddo

nZ=nZ+1
met(nZ)=m

100 continue

!count # of age values
do j=1,nages
   if(a == age(j)) cycle list
enddo

nages=nages+1
age(nages)=a

enddo list

allocate(tssp(nages))
allocate(Zssp(nZ))
allocate(mussp(nimf))

mussp(1:nimf)=mu(1:nimf)
tssp(1:nages)=age(1:nages)
call sort(nZ,met)
Zssp(1:nZ)=met(1:nZ)

close(1)

call system('rm '//trim(adjustl(tmpfile)) )

!read last SSP file and get # of rows = # of wavel
open(unit=10,file=filename, iostat=errstat, action='read')
if (errstat /= 0) then    ! check for errors
   write(*,*) 'Internal error'
   write(*,*) 'Could not open the SSP file', filename
   stop
end if

do k=1,nskip
 read(10,*,IOSTAT=errstat) skip

 if(errstat/=0) then
    write(*,*) 'Internal error'
    write(*,*) 'Error in reading header of SSP file: ', filename
    stop
 end if
enddo

call crow(10,nw) 
write(*,*) 'Number of lines in the spectra:',nw
close(10)


END SUBROUTINE findSSPfiles


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE sort(nn,arr)
USE parameters
IMPLICIT NONE
INTEGER,INTENT(IN)::nn
REAL(kind=dp),INTENT(INOUT),DIMENSION(nn)::arr
!INTEGER,INTENT(OUT),DIMENSION(nn)::ii
INTEGER::x1,x2,ii0
REAL(kind=dp)::aa


  sorting:DO x1=2,nn
     aa=arr(x1)
     !ii0=x1
     DO x2=x1-1,1,-1
        IF(arr(x2)<=aa) EXIT
        arr(x2+1)=arr(x2)
        !ii(x2+1)=x2
     ENDDO
     arr(x2+1)=aa
     !ii(x2+1)=ii0
  ENDDO sorting


END SUBROUTINE sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
