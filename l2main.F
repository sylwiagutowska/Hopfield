SUBROUTINE L2MAIN(cform,nnlo,coord,zz,so,NSPIN1,fnamehelp,helpfiles,vresp_write)
!                                                                       
  USE defs; USE parallel; USE param
  USE atspdt
  USE char; USE com
  USE charp; USE chard; USE charf
  USE lo; USE lohelp
  USE struk
  USE xa; USE xa3; USE xdos
  USE ORB
  USE density_matrix
  IMPLICIT REAL*8 (A-H,O-Z)
#ifdef Parallel
  INCLUDE 'mpif.h'
#endif

  integer          nnlo,nnnlo
  CHARACTER *4     cform                                           
  CHARACTER *5     coord
  CHARACTER *10    BNAME
  CHARACTER *180    VECFN,VECFND,FNAME,fnamehelp,fnamehelp2
  CHARACTER*250    parseline
  LOGICAL          SO, more_kpoints, doos, vresp_write, helpfiles
  REAL*8 ZZ(*),s_kvec,t_kvec,z_kvec
  complex*16 alyl,blyl,dum_hclm
  COMPLEX*16    H_ALM(nume,(LMAX2+1)*(LMAX2+1)), &
                H_BLM(nume,(LMAX2+1)*(LMAX2+1)), &
                H_cLM(nume,(LMAX2+1)*(LMAX2+1),nloat)
  complex*16 dmat_p(2:4,2:4,nume),dmat_d(5:9,5:9,nume),dmat_f(10:16,10:16,nume)
  complex*16 densmat(2:16,2:16)
  complex*16 sdensmat(2:16,2:16)
  COMPLEX*16       YL((LMAX2+1)*(LMAX2+1))
  COMPLEX*16       PHSHEL,CFAC,IMAG1                     
  COMPLEX*16,ALLOCATABLE  :: alm(:,:),blm(:,:),clm(:,:,:)
  COMPLEX*16,ALLOCATABLE  :: alm_buf(:,:),blm_buf(:,:)
!
!bk   93/08/26: parameters, variables, and commons 
!     for force calculation
  complex*16  ayp,byp &
       ,aalm(3,nume,(LMAX2+1)*(LMAX2+1)) &
       ,bblm(3,nume,(LMAX2+1)*(LMAX2+1)) &
       ,cclm(3,nume,(LMAX2+1)*(LMAX2+1),nloat) &
       ,tuu(ngau),tdd(ngau),tud(ngau),tdu(ngau) &
       ,tuu12(ngau),tuu21(ngau),tuu22(ngau) &
       ,tud21(ngau),tdu12(ngau)
  COMPLEX*16               :: aalm_buf(3,nume,(LMAX2+1)*(LMAX2+1)) &
       ,bblm_buf(3,nume,(LMAX2+1)*(LMAX2+1))
  complex*16, allocatable :: aalm_buf_tmp(:,:,:), bblm_buf_tmp(:,:,:)
  integer     lv(ngau),lpv(ngau),mv(ngau),mpv(ngau)
  logical, allocatable   :: forcea(:,:)
  logical, allocatable   :: forcea_buf(:,:)
  logical     force,forout,forb_log1
  REAL*8                   :: vRHOLM(NRAD,NCOM),tc_buf(nrad),vtc_buf(nrad)
  REAL*8                   :: s12(nloat),se12(nloat),s21(nloat),se21(nloat)
  REAL*8                   :: s22(nloat,nloat)
  REAL*8                   :: vs12(nloat),vse12(nloat),vs21(nloat),vse21(nloat)
  REAL*8                   :: vs22(nloat,nloat)

  real*8      k(3),krot(3),fequ(0:3)             ! &
!        ,ftot(0:3,ndif),fvdrho(0:3,ndif) &
!        ,fsph(0:3,ndif),fnsp(0:3,ndif) &
!        ,fpsi(0:3,ndif),fsph2(0:3,ndif) 
        
  REAL*8, ALLOCATABLE      :: ftot(:,:),fvdrho(:,:),fsph(:,:),fnsp(:,:)
  REAL*8, ALLOCATABLE      :: fpsi(:,:),fsph2(:,:),forb(:,:)
  real*8, allocatable      :: fsph_buf(:,:),fnsp_buf(:,:),fsph2_buf(:,:)
  real*8, allocatable      :: forb_buf(:,:),fvdrho_buf(:,:)
!                                                                       
  COMMON /RADFU/  RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2), &
       RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)                 
  COMMON /UHELP/  UA(NRAD),UB(NRAD),UBA(NRAD),UAB(NRAD)
  common /xmean/  xwteh(0:3),xwtel(0:3),xwt1h(0:3),xwt1l(0:3) 
  REAL(8)     :: rhoefg(NRAD,5,23)
  INTEGER     :: lefg(3,5,23),mefg(3,5,23),nefg(5)
  INTEGER     :: lq,lqmax
  common /mkef/delef,ts2
!para begin
  COMMON /PROC/    VECFN,VECFND
  COMMON /IPROC/   IPROC
  logical done_init_xa3
  integer, allocatable :: intebuf_lm(:,:),intebuf(:)
  real*8, allocatable  :: vRHObuf(:,:),RHObuf(:,:)
  character  line*200,text_myid_atm*100
  character, allocatable ::lines(:)*200
!para end

  real*8     h_ekin(iblock),h_k(3,iblock),h_al(iblock),h_bl(iblock)
  complex*16 h_yl((LMAX2+1)*(LMAX2+1),iblock)
  complex*16 h_alyl((LMAX2+1)*(LMAX2+1),iblock)
  complex*16 h_blyl((LMAX2+1)*(LMAX2+1),iblock)
  complex*16, allocatable:: h_ablyl_hk(:,:)

  logical single_kpoint
  real*8 testmax(2),testmax1(2)


  allocate (jatu_eece(nat),latu_eece(nat))
  jatu_eece=0
  latu_eece=0
  natu=0

  allocate (ftot(0:3,ndif),fvdrho(0:3,ndif))
  allocate (fsph(0:3,ndif),fnsp(0:3,ndif))
  allocate (fpsi(0:3,ndif),fsph2(0:3,ndif),forb(0:3,ndif)) 
  allocate (forcea(0:3,nat))
  allocate (fvdrho_buf(0:3,ndif),fsph_buf(0:3,ndif),fnsp_buf(0:3,ndif),fsph2_buf(0:3,ndif))
  allocate (forb_buf(0:3,ndif))
  allocate (forcea_buf(0:3,nat))

  eferm=ef+delef

  IF(MODUS.EQ.'QTL  ')  ef=9.0

!bk   93/09/07: force setup
!     1. existence of tape19,70,71 
!     2. set label force 
  force=.false.
  read(2,'(//)',end=771)
  read(19,'(//)',end=771)
  if(modus.eq.'FOR ')   force=.true.
  IF(myid.EQ.0) write(6,*) '    FORCE-CALCULATION:',force
771 continue

  forout=.false.
  TWOPI=TWO*PI
  TEST1=0.0D0                                                        
  ETOT=0.0D0                                                        
  XWT=0.0                                                           
  SQRT2=SQRT(2.0D0)                                                 
  SQRT3=SQRT(3.D0)                                                  
  SQFP=SQRT(4.D0*PI)                                                
!                                                                       
  CIN=1.d0/137.0359895d0**2                       
  IF (.NOT.REL) CIN=4.0*1.0D-22                                     

  NK=0
      DO JK=1,NKPT                                                   
        NB(JK)=0                                                       
     ENDDO
  forb_log=.false.
  forb_log1=.false.       ! for correctly reading case.vorbup/dn
  read(7,*,end=555,err=555) nmod, nsp, natorb, bext
  do iatorb=1,natorb
    read(7,*,end=555,err=555) iatom(iatorb), nlorb(iatorb)
    do ilorb=1,nlorb(iatorb)
      read(7,*,end=555,err=555) lorb(ilorb,iatorb)
      do m1=-lorb(ilorb,iatorb),lorb(ilorb,iatorb)
        do m2=-lorb(ilorb,iatorb),lorb(ilorb,iatorb)
          read(7,*,end=555,err=555) rval,cval
          vorb(iatorb,lorb(ilorb,iatorb),m1,m2)=cmplx(rval,cval)
        enddo
      enddo
    enddo
  enddo
  forb_log1=.true.
 555  continue


  READ(18,2032) ISCF                                                     
  IF(myid.EQ.0) THEN
     IF(MODUS.EQ.'TOT  '.or.MODUS.EQ.'FOR  '.or.MODUS.EQ.'CLM  ')then
        WRITE(8,787) ISCF                        
        WRITE(8,*) '   NORM OF CLM(R) =          '                        
        WRITE(8,*)                                                 
     ENDIF
     !_vresp
     if (vresp_write)  then       
        WRITE(28,787) ISCF                        
        WRITE(28,*) '   NORM OF CLM(R) =          '                        
        WRITE(28,*)  
     endif                                                      
     !_vresp
     WRITE(6,800)
  ENDIF
  LFIRST=1 
  nnlo=0


  ALLOCATE(alm((lmax2+1)*(lmax2+1),nume),blm((lmax2+1)*(lmax2+1),nume), &
       clm((lmax2+1)*(lmax2+1),nume,nloat))
  ALLOCATE(alm_buf((lmax2+1)*(lmax2+1),nume),blm_buf((lmax2+1)*(lmax2+1),nume))
  CALL init_charp(nume)
  CALL init_chard(nume)
  CALL init_charf(nume)

  IF (MODUS.EQ.'EFG ') THEN                                         
     LMX=3
  ELSE                                                              
     LMX=LMAX2                                                    
  ENDIF


  forcea(0:3,1:nat)=.false.
  forcea_buf(0:3,1:nat)=.false.
  jatom_old=0
  kkk_old=-1
  single_kpoint=.false.



!     ---------------------------------
!     START LOOP FOR ALL ATOMS
!     ---------------------------------
  non_equiv_loop: do jatom_pe=1,nat,npe_atm

     jatom=jatom_pe+myid_atm
     densmat=0.d0
!     densmatd=0.d0


     if (myid_atm.ne.0.and.myid_vec.eq.0) then
        write(text_myid_atm,*) myid_atm
        fname='tmp21'
!        open(unit=21,file=trim(fname)//"_proc_"//trim(adjustl(text_myid_atm)),status='scratch')
        open(unit=21,status='scratch')
        open(unit=22,status='scratch')
     endif

     if(helpfiles) then
       close(31)
       call open_helpfile(fnamehelp,jatom & 
                          ,fnamehelp2)
       open(31,file=fnamehelp2,status='unknown',form='formatted')
     endif

!                                                                      
     lfirst=1
     do i=1,jatom-1
        lfirst=lfirst + mult(i)
     enddo

     if((iproc.gt.0.and.modus.eq.'QTL  ').or.(iproc.gt.0.and.modus.eq.'EFG  ')) then
        close(10)
        call mknam(fname,vecfn,1)
        open(10,file=fname,status='old',form='unformatted',err=950)
     endif

    call cputim(t1c)
    call walltim(t1w)

!                                                                       
! calculate radial functions U(r), UE(R), ...                       
     call atpar (rel,nat,jatom,lfirst,cform,zz(jatom))               


!.......for momentum program #write to *.radwf file
     IF(MODUS.EQ.'ALM  ') then
        write(23,4645)jatom,jri(jatom),r0(jatom),dx(jatom),rmt(jatom)
        do l=0,lmax2
           write(23,*) l
           if (l.le.lomax) then
              write(23,4646) (RRAD1(jrj,l),RRAD2(jrj,l),RADE1(jrj,l),RADE2(jrj,l), &
                   a1lo(jrj,1,l),b1lo(jrj,1,l),a1lo(jrj,2,l),b1lo(jrj,2,l), &
                   a1lo(jrj,3,l),b1lo(jrj,3,l),jrj=1,jri(jatom))
           else
              write(23,4647) (RRAD1(jrj,l),RRAD2(jrj,l),RADE1(jrj,l),RADE2(jrj,l), &
                    jrj=1,jri(jatom))
           endif
        enddo
     endif

!................................

     do i=jatom_old+1,jatom
        read(5,1003) ( (lm(j,jlm),j=1,2), jlm=1,ncom )                    
     enddo
     jatom_old=jatom
! neg L MEANS NEGATIVE SPHERICAL HARMONIC COMB (SEE KURKI-SUONIO)   
     ncomu=0
     DO JLM=2,NCOM
        if(abs(lm(1,jlm)).eq.1) forcea(0,jatom)=.true.
        if((lm(1,jlm).eq.1).and.(lm(2,jlm).eq.1)) forcea(1,jatom)=.true.
        if((lm(1,jlm).eq.-1).and.(lm(2,jlm).eq.1)) forcea(2,jatom)=.true.
        if((lm(1,jlm).eq.1).and.(lm(2,jlm).eq.0)) forcea(3,jatom)=.true.
        IF(LM(1,JLM).EQ.0) GOTO 198                                       
     ENDDO
198  LMMAX=JLM-1                                                       

     IF(myid_vec.EQ.0)  write(6,*) ' atom',jatom,' ncomu',ncomu,' lmmax',lmmax
!                                                                       
     IF(myid_vec.EQ.0) THEN
        WRITE(6,1005) LMMAX,((LM(J,JLM),J=1,2), JLM=1,lmmax)             
        WRITE(21,13) JATOM,IATNR(JATOM),(POS(I,lfirst),I=1,3),MULT(JATOM),ZZ(jatom),aname(jatom)
        if(LMMAX.ne.49) then
            WRITE(21,1005)LMMAX,((LM(J,JLM),J=1,2), JLM=1,lmmax) 
        else
            WRITE(21,1006)LMMAX
        endif            
     ENDIF

     IMAX=JRI(JATOM)
     DO I=1,IMAX
        R(I)=R0(JATOM)*EXP((I-1)*DX(JATOM))                            
     ENDDO

     vrholm(1:IMAX,1:LMMAX)=zero
     rholm(1:IMAX,1:LMMAX)=zero

     IF (MODUS.EQ.'EFG  ') rhoefg(1:jri(jatom),:,:)=zero

     xwt1(0:21)=0.0                                                       
     xwt1l(0:3)=0.0                                                       
     xwt1h(0:3)=0.0                                                       
     xwteh(0:3)=0.0                                                       
     xwtel(0:3)=0.0                                                       

     kkk=0                                                             

!bk   93/09/07: read from tape70,71
!     2. tape71 --> non-spherical matrixelements for [yu91]-force
     if (force) then
        do
           read(2,*) ja
           if (ja.gt.jatom) goto 900
           DO ih=1,ngau
              read(2,'(6i3,4e19.12,/,6(3x),4e19.12)',end=669)  &
                   lpv(ih),ll,lv(ih),mpv(ih),mm,mv(ih), &
                   tuu(ih),tdd(ih),tud(ih),tdu(ih)
              if (ll.eq.0) goto 669
              IF ((LPV(IH) .LE. LOMAX).OR.(LV(IH) .LE. LOMAX)) then 
                 read(2,'(6(3x),6e19.12,/,6(3x),4e19.12)',end=669) &
                      tuu21(ih),tuu12(ih),tuu22(ih),tud21(ih),tdu12(ih)
              ELSE
                 tuu21(ih)=(0.0d0,0.0d0)
                 tuu12(ih)=(0.0d0,0.0d0)
                 tuu22(ih)=(0.0d0,0.0d0)
                 tud21(ih)=(0.0d0,0.0d0)
                 tdu12(ih)=(0.0d0,0.0d0)
              ENDIF
           ENDDO
           goto 910
669        ihmx=ih-1
           write(6,*) 'NGAU for Atom',ja,ihmx
           if (ja.eq.jatom) exit
        enddo
     endif

!.....READ IN K POINT AND BASIS VECTORS                                 
!para begin
! ensure to mimick the proper vector file
! if running on multiple processors
     iloop=0
788  continue
     if((IPROC.GT.0.AND.MODUS.EQ.'QTL  ').or.(IPROC.GT.0.AND.MODUS.EQ.'EFG  ')) then
        iloop=iloop+1
        write(6,*)'ILOOP ',iloop
        close(10)
        call mknam(FNAME,VECFN,ILOOP)
        write(6,*)'opening ',FNAME
        open(10,FILE=FNAME,STATUS='old',FORM='unformatted',ERR=950)
     endif

     nnlo=nlo+nlon+nlov
     n_pw=nmat-(nlo+nlon+nlov)
     ibpp_max=(n_pw-1)/(iblock*npe_vec)+1
     CALL init_xa3(ibpp_max*iblock,nlo+nlon+nlov,nume)

     if (myid_atm.eq.0.and.myid_vec.eq.0) then
        rewind(10)
        do i=1,nat
           read(10) emist
           read(10) emist
        enddo
     endif

! kpoint loop begin     
4    CONTINUE

    if(kkk.eq.0 .and. kkk_old.eq.1) single_kpoint=.true.

    call cputim(t1c)
    call walltim(t1w)

    if(.not.single_kpoint) then 
       CALL read_vec(more_kpoints,nemin,nemax,kkk,n,jatom,ne,etot,so,nspin1,s_kvec,t_kvec,z_kvec,bname,trc,trw)
    else
       if (kkk.eq.0) then
          kkk=1
          more_kpoints=.true.
       else
          more_kpoints=.false.
       endif
    endif
     kkk_old=kkk


     IF(.NOT.more_kpoints) GOTO 998

     IF(myid_vec.EQ.0) THEN
       if(helpfiles) WRITE(31,205) s_kvec,t_kvec,z_kvec,n,ne,bname
     ENDIF

!  write weights on file 17 if case.inso exists
     if(jatom.eq.1.and.myid_vec.eq.0)then
        read(4,*,err=504,end=504)
        rewind 4
        write(17,501)kkk,ne
        write(17,503)(weigh(kkk,iii),iii=1,ne)
504     continue
     endif

!.....CALCULATE  BESSELFUNCTIONS                                        
!     YLM ARE CALCULATED FOR EACH ENERGY SEPARATLY TO CONSERVE CM-SPACE 
!     CALL HARMON(N-(nlo+nlon+nlov),bkx,bky,bkz,lmax2,fj,dfj,rmt(jatom))

     isize=MIN(ibpp*iblock,n-(nlo+nlon+nlov)-myid_vec*ibpp*iblock)

     CALL HARMON(isize,bkx,bky,bkz,lmax2,fj,dfj,rmt(jatom))

!......for momentum program
     IF(MODUS.EQ.'ALM  ') then
           WRITE(24,2055) s_kvec,t_kvec,z_kvec,n,ne,bname
           write(24,*) jatom,nemin,nemax,' jatom,nemin,nemax'
     endif


! zero koeff. for charge analysis
     CALL zero_charp(nemin,nemax)
     CALL zero_chard(nemin,nemax)
     CALL zero_charf(nemin,nemax)
     CALL zerotc_lohelp(nemin,nemax)
     TC100(0:lmax2,nemin:nemax)=zero                                           
     TCA100(0:lmax2,nemin:nemax)=zero
     TCB100(0:lmax2,nemin:nemax)=zero
     CALL zero_xdos(nemin,nemax)
     dmat_p=0.d0
     dmat_d=0.d0
     dmat_f=0.d0

! CALCULATE ALM, BLM                                                
     fac=4.0d0*pi*rmt(jatom)**2/SQRT(vol)                              
     latom=lfirst-1                                                    

     equiv_atom_loop: DO mu=1,mult(jatom)


        latom=latom+1
        alm_buf(1:(lmax2+1)*(lmax2+1),nemin:nemax)=zeroc
        blm_buf(1:(lmax2+1)*(lmax2+1),nemin:nemax)=zeroc
        clm(1:(lmax2+1)*(lmax2+1),nemin:nemax,1:nloat)=zeroc
        IF (force) THEN
!           aalm_buf(1:3,nemin:nemax,1:(lmax2+1)*(lmax2+1))=zeroc
!           bblm_buf(1:3,nemin:nemax,1:(lmax2+1)*(lmax2+1))=zeroc
           aalm_buf_tmp(1:(lmax2+1)*(lmax2+1),nemin:nemax,1:3)=zeroc
           bblm_buf_tmp(1:(lmax2+1)*(lmax2+1),nemin:nemax,1:3)=zeroc
           cclm(1:3,nemin:nemax,1:(lmax2+1)*(lmax2+1),1:nloat)=zeroc
        ENDIF


        blocked_loop: DO ii=1,iblock*ibpp,iblock
           ibb=MIN(iblock,n-(nlo+nlon+nlov)-ii+1-(myid_vec*ibpp*iblock))
           if(ibb.le.0) EXIT
           i3=0
           DO i=ii,ii+iblock-1
              IF(myid_vec*ibpp*iblock+i-1.GT.n-(nlo+nlon+nlov)) EXIT
              i3=i3+1
              bk(1)=bkx(i)
              bk(2)=bky(i)
              bk(3)=bkz(i)

              CALL ROTATE (bk,rotij(1,1,latom),bkrot)

              bk(1)=bkrot(1)*br1(1,1)+bkrot(2)*br1(1,2) &
                   +bkrot(3)*br1(1,3)   
              bk(2)=bkrot(1)*br1(2,1)+bkrot(2)*br1(2,2) &
                   +bkrot(3)*br1(2,3)   
              bk(3)=bkrot(1)*br1(3,1)+bkrot(2)*br1(3,2) &
                   +bkrot(3)*br1(3,3)   
              IF(force.AND.forcea(0,jatom)) THEN 
                 k(1)=kx(i)
                 k(2)=ky(i)
                 k(3)=kz(i)
                 CALL ROTATE(k,rotij(1,1,latom),krot)
                 k(1)=krot(1)*br1(1,1) &
                      +krot(2)*br1(1,2)+krot(3)*br1(1,3)   
                 k(2)=krot(1)*br1(2,1) &
                      +krot(2)*br1(2,2)+krot(3)*br1(2,3)   
                 k(3)=krot(1)*br1(3,1) &
                      +krot(2)*br1(3,2)+krot(3)*br1(3,3)   
                 CALL ROTATE (k,rotloc(1,1,jatom),krot)  
                 h_k(1,i3)=krot(1)                   
                 h_k(2,i3)=krot(2)                   
                 h_k(3,i3)=krot(3)                   
              ENDIF

              CALL ROTATE (bk,rotloc(1,1,jatom),bkrloc)

              CALL YLM (bkrloc,lmax2,yl)
              arg1=bkrot(1)*pos(1,lfirst)*twopi
              arg2=bkrot(2)*pos(2,lfirst)*twopi
              arg3=bkrot(3)*pos(3,lfirst)*twopi
              argt=(bkx(i)*tauij(1,latom)+bky(i)*tauij(2,latom)+bkz(i)*tauij(3,latom))*twopi

              phshel=EXP( imag*(arg1+arg2+arg3+argt) )
              index=0
              DO  l=0,lmx
                 maxx=2*l+1
                 DO  m=1,maxx
                    index=index+1
                    h_yl(index,i3)=CONJG(yl(index))*phshel
                 ENDDO
              ENDDO
           ENDDO

           index=0
           rmt2=1.D0/(rmt(jatom)**2)
           DO l=0,lmx
              i3=0
              DO i=ii,ii+iblock-1
                 IF(myid_vec*ibpp*iblock+i.GT.n-(nlo+nlon+nlov)) EXIT
                 i3=i3+1
                 if(lapw(l)) then
                    h_al(i3)=dfj(l,i)*pe(l)-fj(l,i)*dpe(l) 
                    h_bl(i3)=fj(l,i)*dp(l)-dfj(l,i)*p(l)
                 ELSE
!                    h_al(i3) = fj(l,i)/p(l)/rmt(jatom)**2
                    h_al(i3) = rmt2*fj(l,i)/p(l)
                    h_bl(i3) = 0.d0
                 ENDIF
              ENDDO
              maxx=2*l+1
              DO m=1,maxx
                 index=index+1
                 i3=0
                 DO i=ii,ii+iblock-1
                    IF(myid_vec*ibpp*iblock+i.GT.n-(nlo+nlon+nlov)) EXIT
                    i3=i3+1
                    h_alyl(index,i3)=h_al(i3)*h_yl(index,i3)
                    h_blyl(index,i3)=h_bl(i3)*h_yl(index,i3)
                 ENDDO
              ENDDO
           ENDDO
!_REAL           lda=2*(LMAX2+1)*(LMAX2+1)
!_COMPLEX           lda=(LMAX2+1)*(LMAX2+1)
           ldc=lda
!           ldb=nmat
           ldb=ibpp_max*iblock

!_REAL           CALL dgemm('N','N',2*index,nemax-nemin+1,ibb,1.d0, &
!_REAL            h_alyl,lda,a(ii,nemin),ldb,1.d0,alm_buf(1,nemin),ldc)

!_REAL           CALL dgemm('N','N',2*index,nemax-nemin+1,ibb,1.d0, &
!_REAL            h_blyl,lda,a(ii,nemin),ldb,1.d0,blm_buf(1,nemin),ldc)

!_COMPLEX           CALL zgemm('N','N',index,nemax-nemin+1,ibb,(1.d0,0.d0), &
!_COMPLEX            h_alyl,lda,a(ii,nemin),ldb,(1.d0,0.d0),alm_buf(1,nemin),ldc)
!_COMPLEX           CALL zgemm('N','N',index,nemax-nemin+1,ibb,(1.d0,0.d0), &
!_COMPLEX            h_blyl,lda,a(ii,nemin),ldb,(1.d0,0.d0),blm_buf(1,nemin),ldc)

              maxindex=(lmx+1)**2
           IF (force.AND.forcea(0,jatom)) THEN 

!_REAL        lda=2*maxindex
!_COMPLEX     lda=maxindex
              ldb=ibpp_max*iblock
              ldc=lda
              
              do i_h_k=1,3
                 do index=1,maxindex              
                    do i3=1,ibb
                       h_ablyl_hk(index,i3)=h_alyl(index,i3)*h_k(i_h_k,i3)
                    enddo
                 enddo
!                 aalm_buf_tmp=0.0d0
!_REAL           CALL dgemm('N','N',2*maxindex,nemax-nemin+1,ibb,1.d0, &
!_REAL                   h_ablyl_hk,lda,a(ii,nemin),ldb,1.d0,aalm_buf_tmp(1,nemin,i_h_k),ldc)
!_COMPLEX        CALL zgemm('N','N',maxindex,nemax-nemin+1,ibb,(1.d0,0.d0), &
!_COMPLEX                h_ablyl_hk,lda,a(ii,nemin),ldb,(1.d0,0.d0),aalm_buf_tmp(1,nemin,i_h_k),ldc)
!                 do num=nemin,nemax
!                    do index=1,maxindex
!                       aalm_buf(i_h_k,num,index)=aalm_buf(i_h_k,num,index)+aalm_buf_tmp(index,num)
!                    enddo
!                 enddo
              enddo

              do i_h_k=1,3
                 do index=1,maxindex              
                    do i3=1,ibb
                       h_ablyl_hk(index,i3)=h_blyl(index,i3)*h_k(i_h_k,i3)
                    enddo
                 enddo
                 
!                 aalm_buf_tmp=0.0d0 ! Why?, use 0.0D0 inside dgemm ?
!_REAL           CALL dgemm('N','N',2*maxindex,nemax-nemin+1,ibb,1.d0, &
!_REAL                   h_ablyl_hk,lda,a(ii,nemin),ldb,1.d0,bblm_buf_tmp(1,nemin,i_h_k),ldc)
!_COMPLEX        CALL zgemm('N','N',maxindex,nemax-nemin+1,ibb,(1.d0,0.d0), &
!_COMPLEX                h_ablyl_hk,lda,a(ii,nemin),ldb,(1.d0,0.d0),bblm_buf_tmp(1,nemin,i_h_k),ldc)
!                 do num=nemin,nemax
!                    do index=1,maxindex
!                       bblm_buf(i_h_k,num,index)=bblm_buf(i_h_k,num,index)+aalm_buf_tmp(index,num)
!                    enddo
!                 enddo
              enddo

           ENDIF
        ENDDO blocked_loop

        if(force)then
                do i_h_k=1,3
                 do num=nemin,nemax
                    do index=1,maxindex
                       aalm_buf(i_h_k,num,index)=aalm_buf_tmp(index,num,i_h_k) 
                       bblm_buf(i_h_k,num,index)=bblm_buf_tmp(index,num,i_h_k)
                    enddo
                 enddo
                enddo
        endif

      ldc=(LMAX2+1)*(LMAX2+1)

        IF(myid_vec.NE.0) CYCLE
#else
        alm(1:ldc,nemin:nemax)=alm_buf(1:ldc,nemin:nemax)
        blm(1:ldc,nemin:nemax)=blm_buf(1:ldc,nemin:nemax)
        IF(force.AND.forcea(0,jatom)) THEN
           aalm(1:3,nemin:nemax,1:ldc)=aalm_buf(1:3,nemin:nemax,1:ldc)
           bblm(1:3,nemin:nemax,1:ldc)=bblm_buf(1:3,nemin:nemax,1:ldc)
        ENDIF
#endif

        if (nlo.ne.0) call lomain (nemin,nemax,lfirst,latom,n,jatom, &
             alm,blm,clm,aalm,bblm,cclm,force,forcea(0,jatom))

        index=0
        DO l=0,lmax2                                                   
           maxx=2*l+1
           cfac=imag**l
           cfac=cfac*fac
           DO m=1,maxx
              index=index+1
              DO num=nemin,nemax
                 alm(index,num)=alm(index,num)*cfac
                 blm(index,num)=blm(index,num)*cfac
!if(jatom.eq.2 .and. l.eq.2 .and. num.eq.50 ) then
!write(998,200) mu,num,index,alm(index,num)
! 200   format('mmu=',i1,' num=',i3,' m=',i2,2f10.4,3x,2f10.4,3x,f10.4)
!endif
              ENDDO
              DO jlo=1,ilo(l)
                 DO num=nemin,nemax
                    clm(index,num,jlo)=clm(index,num,jlo)*cfac
                 ENDDO
              ENDDO
              IF (force.AND.forcea(0,jatom)) THEN
                 DO num=nemin,nemax
                    aalm(1,num,index)=aalm(1,num,index)*cfac
                    aalm(2,num,index)=aalm(2,num,index)*cfac
                    aalm(3,num,index)=aalm(3,num,index)*cfac
                    
                    bblm(1,num,index)=bblm(1,num,index)*cfac
                    bblm(2,num,index)=bblm(2,num,index)*cfac
                    bblm(3,num,index)=bblm(3,num,index)*cfac
                 ENDDO
                 DO jlo=1,ilo(l)
                    DO num=nemin,nemax
                       cclm(1,num,index,jlo)=cclm(1,num,index,jlo)*cfac
                       cclm(2,num,index,jlo)=cclm(2,num,index,jlo)*cfac
                       cclm(3,num,index,jlo)=cclm(3,num,index,jlo)*cfac
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
!... Umspeichern von H_ALM auf ALM und H_BLM auf BLM
        h_clm=0.0d0
        CALL c_transpose(alm,(LMAX2+1)*(LMAX2+1),(LMAX2+1)*(LMAX2+1), &
             nemin,nemax,h_alm,nume)
        CALL c_transpose(blm,(LMAX2+1)*(LMAX2+1),(LMAX2+1)*(LMAX2+1), &
             nemin,nemax,h_blm,nume)
        DO jlo=1,nloat
           call c_transpose(clm(1,1,jlo),(LMAX2+1)*(LMAX2+1),(LMAX2+1)*(LMAX2+1), &
                nemin,nemax,h_clm(1,1,jlo),nume)
        ENDDO

!......for momentum program
        IF(MODUS.EQ.'ALM  ') then
           WRITE(24,*) mu,'  ATOM'
           DO NUM=NEMIN,NEMAX
              write(24,*) num,WEIGHT(NUM),' NUM, weight'
              INDEX=0
              DO L=0,LMAX2
                 DO M=-l,l
                    INDEX=INDEX+1
                    write(24,4893)l,m,index,alm(INDEX,NUM),blm(INDEX,NUM),(clm(INDEX,NUM,jlo),jlo=1,3)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

!..........................................
!                                                                       
!.... MAIN FORCE CALCULATIONS
!bk   93/08/26: 
        IF (force.AND.forcea(0,jatom).and.myid_vec.eq.0) CALL FOMAI1(latom,jatom, &
             nemin,nemax,lmx,ihmx,ALM,blm,clm, &
             aalm,bblm,cclm,tuu,tdd,tud, &
             tdu,tuu12,tuu21,tuu22,tud21,tdu12,lv,lpv, &
             mv,mpv,forout,fsph,fnsp,fsph2,forb)

!.....C(L,M) TO BE CALCULATED                                           
        LQMAX=0          
        IF(MODUS.EQ.'QTL ') LMMAX=1
        LQ=0

        
        ilm_loop: DO ILM=1,LMMAX
           LL=IABS(LM(1,ILM))                                                
           MM=LM(2,ILM)                                                      
           IF (modus.EQ.'EFG  ') THEN
              IF(ll.eq.2) THEN                            
                 lq=lq+1
                 lqmax=lq 
                 nefg(lq)=0
              ELSE if (ll.ne.0) then
                 CYCLE
              ENDIF
           ENDIF
!     SEE KURKI-SUONI  FACTOR FOR  LM+ (1,0), LM-(0,-1), M.NE.2N *(-1)  
           IMAG1=(1.0,0.0)                                                   
!     conjg fix
!      IF(LM(1,ILM).LT.0) IMAG1=-IMAG                                    
   
              IF(LM(1,ILM).LT.0) IMAG1= IMAG
              IF(MOD(MM,2).EQ.1) IMAG1=-IMAG1                     
              l_sum: DO L=0,lmx
                    L1=L+1
                    MMAX=2*L+1
                    lp_sum: DO LP=0,lmx
                          LP1=LP+1
                          IF(NOTRI(LL,L,LP).LT.0) CYCLE
                          MPMAX=2*LP+1
                          IMAX=JRI(JATOM)
                          DO I=1,IMAX
                             UA(I)=RRAD1(I,L)*RRAD1(I,LP)+ &
                                  CIN*RRAD2(I,L)*RRAD2(I,LP)    
                             UB(I)=RADE1(I,L)*RADE1(I,LP)+ &
                                  CIN*RADE2(I,L)*RADE2(I,LP)    
                             UAB(I)=RRAD1(I,L)*RADE1(I,LP)+ &
                                  CIN*RRAD2(I,L)*RADE2(I,LP)   
                             UBA(I)=RRAD1(I,LP)*RADE1(I,L)+ &
                                  CIN*RRAD2(I,LP)*RADE2(I,L)   
                          ENDDO
                          
                          DO jlo=1,ilo(l)
                             IF(loor(jlo,l)) THEN
                                DO I=1,IMAX
                                   U21(I,jlo)=a1lo(I,jlo,L)*RRAD1(I,LP) &
                                        +CIN*b1lo(I,jlo,L)*RRAD2(I,LP)    
                                   Ue21(I,jlo)=a1lo(I,jlo,L)*RADE1(I,LP) &
                                        +CIN*b1lo(I,jlo,L)*RADE2(I,LP)    
                                ENDDO
                             ENDIF
                          ENDDO
                          DO jlop=1,ilo(lp)
                             IF(loor(jlop,lp)) THEN
                                DO i=1,imax
                                   u12(i,jlop)=rrad1(i,l)*a1lo(i,jlop,lp) &
                                        +cin*rrad2(i,l)*b1lo(i,jlop,lp)    
                                   ue12(i,jlop)=rade1(i,l)*a1lo(i,jlop,lp) &
                                        +cin*rade2(i,l)*b1lo(i,jlop,lp) 
                                ENDDO
                             ENDIF
                          ENDDO
                          DO jlo=1,ilo(l)
                             DO jlop=1,ilo(lp) 
                                IF(loor(jlo,l).AND.loor(jlop,lp)) THEN
                                   DO i=1,imax
                                      u22(i,jlop,jlo)=a1lo(i,jlo,l)*a1lo(i,jlop,lp) &
                                           +cin*b1lo(i,jlo,l)*b1lo(i,jlop,lp) 
                                   ENDDO
                                ENDIF
                             ENDDO
                          ENDDO
                          !     
                          !.....M SUM                           
                          suma(nemin:nemax) = zero
                          sumb(nemin:nemax) = zero
                          sumab(nemin:nemax)= zero
                          sumba(nemin:nemax)= zero
                          sum21(nemin:nemax,1:ilo(l)) = zero
                          sume21(nemin:nemax,1:ilo(l))= zero
                          sum12(nemin:nemax,1:ilo(lp)) = zero
                          sume12(nemin:nemax,1:ilo(lp))= zero
                          sum22(nemin:nemax,1:ilo(lp),1:ilo(l))=zero

                          
                          M=-L1 
                          MP=-LP1
                          m_loop: DO MS=1,MMAX
                             M=M+1
                             MP=-LP1
                             mp_loop: DO MPS=1,MPMAX
                                MP=MP+1
                                MTEST=-M+MM+MP
                                IF(MTEST.NE.0) CYCLE
                                LY=L*L1+M+1
                                LPY=LP*LP1+MP+1
                                GNT=GAUNT(L,LL,LP,M,MM,MP)
                                IF (MODUS.EQ.'EFG  '.AND.LL.EQ.2) THEN
                                   NEFG(LQ)=NEFG(LQ)+1
                                   LEFG(1,LQ,NEFG(LQ))=LM(1,ILM)
                                   LEFG(2,LQ,NEFG(LQ))=L
                                   LEFG(3,LQ,NEFG(LQ))=LP
                                   MEFG(1,LQ,NEFG(LQ))=MM
                                   MEFG(2,LQ,NEFG(LQ))=M
                                   MEFG(3,LQ,NEFG(LQ))=MP
                                END IF
                                DO NUM=NEMIN,NEMAX
                                   SA=(h_ALM(num,LY)*dconjg(h_ALM(num,LPY)) &
                                        *IMAG1)*GNT
                                   SB=(h_BLM(num,LY)*dconjg(h_BLM(num,LPY)) &
                                        *IMAG1)*GNT
                                   SAB=(h_ALM(num,LY)*dconjg(h_BLM(num,LPY)) &
                                        *IMAG1)*GNT
                                   SBA=(h_BLM(num,LY)*dconjg(h_ALM(num,LPY)) &
                                        *IMAG1)*GNT
                                   DO jlo=1,ilo(l)
                                      S21(jlo)=(h_cLM(num,LY,jlo)*dconjg(h_ALM(num,LPY)) &
                                           *IMAG1)*GNT
                                      Se21(jlo)=(h_cLM(num,LY,jlo)*dconjg(h_BLM(num,LPY)) &
                                           *IMAG1)*GNT
                                   ENDDO
                                   DO jlop=1,ilo(lp)
                                      S12(jlop)=h_ALM(num,LY)*dconjg(h_cLM(num,LPY,jlop))*IMAG1*GNT
                                      Se12(jlop)=h_BLM(num,LY)*dconjg(h_cLM(num,LPY,jlop))*IMAG1*GNT
                                   ENDDO
                                   DO jlo=1,ilo(l)
                                      DO jlop=1,ilo(lp)
                                         S22(jlop,jlo)=(h_cLM(num,LY,jlo)*dconjg(h_cLM(num,LPY,jlop)) &
                                              *IMAG1)*GNT
                                      ENDDO
                                   ENDDO
                                   
                                   SUMA(NUM)=SUMA(NUM)+SA
                                   SUMB(NUM)=SUMB(NUM)+SB
                                   SUMAB(NUM)=SUMAB(NUM)+SAB
                                   SUMBA(NUM)=SUMBA(NUM)+SBA
                                   DO jlo=1,ilo(l)
                                      SUM21(NUM,jlo)=SUM21(NUM,jlo)+S21(jlo)
                                      SUMe21(NUM,jlo)=SUMe21(NUM,jlo)+Se21(jlo)
                                   ENDDO
                                   DO jlop=1,ilo(lp)
                                      SUM12(NUM,jlop)=SUM12(NUM,jlop)+S12(jlop)
                                      SUMe12(NUM,jlop)=SUMe12(NUM,jlop)+Se12(jlop)
                                   ENDDO
                                   DO jlo=1,ilo(l)
                                      DO jlop=1,ilo(lp)
                                         SUM22(NUM,jlop,jlo)=SUM22(NUM,jlop,jlo)+S22(jlop,jlo)
                                      ENDDO
                                   ENDDO
                                   !                                 
                                   IF (MODUS.EQ.'EFG  '.AND.LL.EQ.2) THEN
                                      DO JR=1,JRI(JATOM)
                                         Q=SA*UA(JR)+SB*UB(JR)+SAB*UAB(JR)+ &
                                              SBA*UBA(JR)
                                         DO jlo=1,ilo(l)
                                            Q=Q+s21(jlo)*u21(jr,jlo) + se21(jlo)*ue21(jr,jlo)
                                         ENDDO
                                         DO jlop=1,ilo(lp)
                                            Q=Q+S12(jlop)*U12(jr,jlop) + se12(jlop)*ue12(jr,jlop)
                                         ENDDO
                                         DO jlo=1,ilo(l)
                                            DO jlop=1,ilo(lp)
                                               Q=Q+s22(jlop,jlo)*u22(jr,jlop,jlo)
                                            ENDDO
                                         ENDDO
                                         RHOEFG(JR,LQ,NEFG(LQ))= &
                                              RHOEFG(JR,LQ,NEFG(LQ))+      &
                                              Q/MULT(JATOM)*WEIGHT(NUM)
                                      ENDDO
                                   ENDIF
                                ENDDO
            

                          IF(LL.eq.0) call csplit(nemin,nemax,l,jatom,mu, &
                               alm,blm,clm,coord,dmat_p,dmat_d,dmat_f) 

                          sa=0.0d0
                          sb=0.0d0
                          sab=0.0d0
                          sba=0.0d0
                          s12=0.0d0
                          se12=0.0d0
                          s21=0.0d0
                          se21=0.0d0
                          s22=0.0d0
                          !_vresp
                          vsa=0.0d0
                          vsb=0.0d0
                          vsab=0.0d0
                          vsba=0.0d0
                          vs12=0.0d0
                          vse12=0.0d0
                          vs21=0.0d0
                          vse21=0.0d0
                          vs22=0.0d0
                          NEDo=nemax-nemin+1
                          sa =ddotK(NEDo,suma(nemin),1,weight(nemin),1)
                          sb =ddotK(NEDo,sumb(nemin),1,weight(nemin),1)
                          sab=ddotK(NEDo,sumab(nemin),1,weight(nemin),1)
                          sba=ddotK(NEDo,sumba(nemin),1,weight(nemin),1)
                          DO jlo=1,ilo(l)
                                s21(jlo )=ddotk(NEDo,sum21(nemin,jlo ),1,weight(nemin),1)
                                se21(jlo)=ddotk(NEDo,sume21(nemin,jlo),1,weight(nemin),1)
                          ENDDO
                          DO jlop=1,ilo(lp)
                                s12(jlop )=ddotk(NEDo,sum12(nemin,jlop ),1,weight(nemin),1)
                                se12(jlop)=ddotk(NEDo,sume12(nemin,jlop),1,weight(nemin),1)
                          ENDDO
                          DO jlo=1,ilo(l)
                                DO jlop=1,ilo(lp)
                                   s22(jlop,jlo)=ddotk(NEDo,sum22(nemin,jlop,jlo),1,weight(nemin),1)
                                ENDDO
                           ENDDO


                          DO num=nemax,nemin,-1 !nemin,nemax 
                           if(abs(weight(num)).gt.1D-18)then
                             facv = (-e(num))*weight(num)
                             vsa=vsa+suma(num)  *facv
                             vsb=vsb+sumb(num)  *facv
                             vsab=vsab+sumab(num)  *facv
                             vsba=vsba+sumba(num)  *facv
                             DO jlo=1,ilo(l)
                                vs21(jlo)=vs21(jlo)+sum21(num,jlo)  *facv
                                vse21(jlo)=vse21(jlo)+sume21(num,jlo)  *facv
                             ENDDO
                             DO jlop=1,ilo(lp)
                                vs12(jlop)=vs12(jlop)+sum12(num,jlop)  *facv
                                vse12(jlop)=vse12(jlop)+sume12(num,jlop)  *facv
                             ENDDO
                             DO jlo=1,ilo(l)
                                DO jlop=1,ilo(lp)
                                   vs22(jlop,jlo)=vs22(jlop,jlo)+sum22(num,jlop,jlo)*facv
                                ENDDO
                             ENDDO
                           endif
                          ENDDO
                          


                          IMAX=JRI(JATOM)
                          DO I=1,IMAX
                             TC_buf(i)=SA*UA(I)+SB*UB(I)+SAB*UAB(I)+SBA*UBA(I)
                             vTC_buf(i)=vSA*UA(I)+vSB*UB(I)+vSAB*UAB(I)+vSBA*UBA(I)
!                         ENDDO
                          DO jlo=1,ilo(l)
!                            DO i=1,imax
                                tc_buf(i)=tc_buf(i)+S21(jlo)*U21(I,jlo)+Se21(jlo)*Ue21(I,jlo)
                                vtc_buf(i)=vtc_buf(i)+vS21(jlo)*U21(I,jlo)+vSe21(jlo)*Ue21(I,jlo)
!                            ENDDO
                          ENDDO
                          DO jlop=1,ilo(lp)
!                            DO i=1,imax
                                tc_buf(i)=tc_buf(i)+S12(jlop)*U12(I,jlop)+Se12(jlop)*Ue12(I,jlop)
                                vtc_buf(i)=vtc_buf(i)+vS12(jlop)*U12(I,jlop)+vSe12(jlop)*Ue12(I,jlop)
!                            ENDDO
                          ENDDO
                          DO jlo=1,ilo(l)
                             DO jlop=1,ilo(lp)
!                               DO i=1,imax
                                   tc_buf(i)=tc_buf(i)+s22(jlop,jlo)*u22(i,jlop,jlo)
                                   vtc_buf(i)=vtc_buf(i)+vs22(jlop,jlo)*u22(i,jlop,jlo)
!                               ENDDO
                             ENDDO
                          ENDDO
!                         DO i=1,imax
                             RHOLM(I,ILM)=RHOLM(I,ILM)+TC_buf(i)/dble(MULT(JATOM))
                             vRHOLM(I,ILM)=vRHOLM(I,ILM)+vTC_buf(i)/dble(MULT(JATOM))
                          ENDDO
                          CALL cputim(time4)
                          time_rad=time_rad+time4-time3
                       endif                      !EECE
                    ENDDO lp_sum
                 endif                      !EECE
              ENDDO l_sum
           endif                      !EECE
        ENDDO ilm_loop



        CALL cputim(time2)
        CALL walltim(time2_w)
        time_ilm=time_ilm+time2-time1
        time_ilm_w=time_ilm_w+time2_w-time1_w
     ENDDO equiv_atom_loop


     IF(myid_vec.EQ.0) call psplit(jatom,nemin,nemax,test1,kkk,EQBAD,jatombad,lbad,helpfiles,dmat_p,dmat_d,dmat_f,densmat) !,densmatd)
     
     GOTO 4
!  kpoint loop end

 998     CONTINUE
     CALL cputim(time1)
     CALL walltim(time1_w)

     IF(myid_vec.NE.0) CYCLE

         IF((MODUS.EQ.'QTL  '.AND.IPROC.GT.0).or. &
              (MODUS.EQ.'EFG  '.AND.IPROC.GT.0)) THEN
            IF(ILOOP.LT.IPROC) goto 788
         ENDIF

         DO ILM=1,LMMAX  
            DIV=XWT                                                           
            OVRDIV=1.0                                                        
            LL=LM(1,ILM)                                                      
            MM=LM(2,ILM)                                                      
            IF(LL.NE.0) GOTO 56                                               
            CALL CHARGEL2(R,OVRDIV,SQFP,RHOLM(1,ILM),DX(JATOM),JRI(JATOM),TC)
            WRITE(6,207)  JATOM,TC                                            
            WRITE(21,720) jatom,JATOM,TC,rmt(jatom)
            IF(ISPLIT(JATOM).EQ.1) THEN                                    
               WRITE(6,250) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,250) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.2) THEN                               
               WRITE(6,251) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,251) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.3) THEN                               
               WRITE(6,252) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,252)jatom,JATOM, JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.-2) THEN                              
               WRITE(6,253) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,253) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.4) THEN                               
               WRITE(6,254) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,254) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.5) THEN                               
               WRITE(6,255) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,255) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.6) THEN                               
               WRITE(6,256) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,256) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.8) THEN                               
               WRITE(6,257) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
               WRITE(21,257) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,14)
            ELSE IF(ISPLIT(JATOM).EQ.15) THEN                               
               WRITE(6,258) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                    ,(XWT1(I),I=7,21)
               WRITE(21,258) jatom,JATOM,JATOM,(XWT1(I),I=0,3) &
                 ,(XWT1(I),I=7,21)
            ENDIF
               do i=0,3
               if(xwt1l(i).lt.0.00005d0) then
                   XWT1l(I)=0.00001d0
                   XWTel(I)=XWT1l(I)*10.00001d0
               endif
               if(xwt1h(i).lt.0.00005d0) then
                   XWT1h(I)=0.00001d0
                   XWTeh(I)=XWT1h(I)*10.00001d0
               endif
               enddo
               write(6,260)                                                   
               WRITE(6,261) jatom,(XWT1l(I),XWTel(I)/XWT1l(I),I=0,3)
               write(6,262)                                                  
               WRITE(6,263) jatom,(XWT1h(I),XWTeh(I)/XWT1h(I),I=0,3)
               write(21,260)                                                   
               WRITE(21,261) jatom,(XWT1l(I),XWTel(I)/XWT1l(I),I=0,3)
               write(21,262)                                                  
               WRITE(21,263) jatom,(XWT1h(I),XWTeh(I)/XWT1h(I),I=0,3)
!
               write(6,*)
               write(6,557) '    ',1
               DO I=2,4
                 WRITE(6,567)(REAL(DensMAT(I,J)),J=2,4)
               enddo
               write(6,556)
               DO I=2,4
                 WRITE(6,567)(dimag(DensMAT(I,J)),J=2,4)
               enddo
               write(6,557) '    ',2
               DO I=5,9
                 WRITE(6,567)(REAL(DensMAT(I,J)),J=5,9)
               enddo
               write(6,556)
               DO I=5,9
                 WRITE(6,567)(dimag(DensMAT(I,J)),J=5,9)
               enddo
               write(6,557) '    ',3
               DO I=10,16
                 WRITE(6,567)(REAL(DensMAT(I,J)),J=10,16)
               enddo
               write(6,556)
               DO I=10,16
                 WRITE(6,567)(dimag(DensMAT(I,J)),J=10,16)
               enddo
556            format(' Density matrix, imag part')
557            format(' Density matrix ',a4,' block, real part.  L=',i2)
567            format(7x,7f9.5)
               If(nspin1.eq.2) then    ! only in spinpolarized mode
               if(lorb_atom(jatom,1)) then
                       CALL SYMMETRIZATION(jatom,1,sdensmat,DensMAT)
                 write(21,*) 
                 write(21,557) '    ',1
                 DO I=2,4
                   WRITE(21,567)(REAL(sDensMAT(I,J)),J=2,4)
                 enddo
                 write(21,556)
                 DO I=2,4
                   WRITE(21,567)(dimag(sDensMAT(I,J)),J=2,4)
                 enddo
                 write(6,*) 
                 write(6,557) '    ',1
                 DO I=2,4
                   WRITE(6,567)(REAL(sDensMAT(I,J)),J=2,4)
                 enddo
                 write(6,556)
                 DO I=2,4
                   WRITE(6,567)(dimag(sDensMAT(I,J)),J=2,4)
                 enddo
                 write(22,7100) jatom
                 write(22,7102) 1,0.d0,0.d0,0.d0
                 do i=2,4
                   write(22,7101)(sdensmat(i,j),j=2,4)
                 enddo
               endif                 
               if(lorb_atom(jatom,2)) then
                       CALL SYMMETRIZATION(jatom,2,sdensmat,DensMAT)
                 write(21,*) 
                 write(21,557) '    ',2
                 DO I=5,9
                   WRITE(21,567)(REAL(sDensMAT(I,J)),J=5,9)
                 enddo
                 write(21,556)
                 DO I=5,9
                   WRITE(21,567)(dimag(sDensMAT(I,J)),J=5,9)
                 enddo
                 write(6,*) 
                 write(6,557) '    ',2
                 DO I=5,9
                   WRITE(6,567)(REAL(sDensMAT(I,J)),J=5,9)
                 enddo
                 write(6,556)
                 DO I=5,9
                   WRITE(6,567)(dimag(sDensMAT(I,J)),J=5,9)
                 enddo
                 write(22,7100) jatom
                 write(22,7102) 2,0.d0,0.d0,0.d0
                 do i=5,9
                   write(22,7101)(sdensmat(i,j),j=5,9)
                 enddo
               endif                 
               if(lorb_atom(jatom,3)) then
                       CALL SYMMETRIZATION(jatom,3,sdensmat,DensMAT)
                 write(21,*) 
                 write(21,557) '    ',3
                 DO I=10,16
                   WRITE(21,567)(REAL(sDensMAT(I,J)),J=10,16)
                 enddo
                 write(21,556)
                 DO I=10,16
                   WRITE(21,567)(dimag(sDensMAT(I,J)),J=10,16)
                 enddo
                 write(6,*) 
                 write(6,557) '    ',3
                 DO I=10,16
                   WRITE(6,567)(REAL(sDensMAT(I,J)),J=10,16)
                 enddo
                 write(6,556)
                 DO I=10,16
                   WRITE(6,567)(dimag(sDensMAT(I,J)),J=10,16)
                 enddo
                 write(22,7100) jatom
                 write(22,7102) 3,0.d0,0.d0,0.d0
                 do i=10,16
                   write(22,7101)(sdensmat(i,j),j=10,16)
                 enddo
               endif                 
               endif                 
7100              format(i5,' atom density matrix')
7101              format(2(2es20.12,2x))
7102              format(i5,3f10.6, ' L, Lx,Ly,Lz in global orthogonal system')


!
 56         IF(MM.EQ.0) GOTO 502                                              
            IMAX=JRI(JATOM)                                                   
            DO I=1,IMAX                                                    
!_vresp
               vRHOLM(I,ILM)=vRHOLM(I,ILM)*SQRT2
!_vresp
               RHOLM(I,ILM)=RHOLM(I,ILM)*SQRT2
            enddo
 502        continue
         enddo

         call cputim(time3)
         call walltim(time3_w)

         if (myid_atm.eq.0) then
            if(MODUS.EQ.'TOT  '.or.MODUS.EQ.'FOR  '.or.MODUS.EQ.'CLM  ')then
               WRITE(8,1990) JATOM                                           
               WRITE(8,2001) LMMAX-ncomu                                      
               DO ILM=1,LMMAX  
                  LL=LM(1,ILM)                                                 
                  MM=LM(2,ILM)       
                  WRITE(8,2011) LL,MM
                  WRITE(8,2022) ( RHOLM(I,ILM), I=1,IMAX )
                  WRITE(8,2031)
               end DO
               WRITE(8,2030)                                                  
            endif
            if(vresp_write) then
               WRITE(28,1990) JATOM                                            
               WRITE(28,2001) LMMAX-ncomu                                     
               DO ILM=1,LMMAX  
                  LL=LM(1,ILM)                                                 
                  MM=LM(2,ILM)       

                  WRITE(28,2011) LL,MM                                         
                  WRITE(28,2022) ( vRHOLM(I,ILM), I=1,IMAX )
                  WRITE(28,2031)                                               
               end DO
               WRITE(28,2030)      
            endif                                               
         endif

#ifdef Parallel

         intebuf(1)=jatom
         intebuf(2)=LMMAX-ncomu
         intebuf(3)=LMMAX
         intebuf(4)=IMAX

         if (myid_atm.eq.0) then
            do ipid=1,npe_atm-1
               if (intebuf(1).eq.nat) exit
               call mpi_recv( intebuf, 4, mpi_integer, ipid, mpi_any_tag, mpi_vec_comm,statusmpi, ierr)                                
                              
               allocate(intebuf_lm(2,intebuf(3)),RHObuf(1:intebuf(4),&
                 1:intebuf(3)),vRHObuf(1:intebuf(4),1:intebuf(3)))

               call mpi_recv( intebuf_lm, 2*intebuf(3), mpi_integer, &
                 ipid, mpi_any_tag, mpi_vec_comm,statusmpi, ierr)           
               if(vresp_write) then
               call mpi_recv( vRHObuf, intebuf(4)*intebuf(3), mpi_double_precision, &
                 ipid, mpi_any_tag, mpi_vec_comm,statusmpi, ierr)
               endif
               if(MODUS.EQ.'TOT  '.or.MODUS.EQ.'FOR  '.or.MODUS.EQ.'CLM  ')then
                  call mpi_recv( RHObuf, intebuf(4)*intebuf(3), mpi_double_precision, &
                    ipid, mpi_any_tag, mpi_vec_comm,statusmpi, ierr)
                  WRITE(8,1990) intebuf(1)                             
                  WRITE(8,2001) intebuf(2)                             
                  DO ILM=1,intebuf(3)  
                     IF(MODUS1.eq.'EECE ') then
                        if(iabs(intebuf_LM(1,ilm)).gt.2*latu_eece(intebuf(1))) cycle
                     endif
                     LL=intebuf_LM(1,ILM)                                
                     MM=intebuf_LM(2,ILM)       
                     WRITE(8,2011) LL,MM
                     WRITE(8,2022) ( RHObuf(I,ILM), I=1,intebuf(4) )
                     WRITE(8,2031)
                  end DO
                  WRITE(8,2030)                                         
               endif

            if(vresp_write) then
               WRITE(28,1990) intebuf(1)                          
               WRITE(28,2001) intebuf(2)                       
               DO ILM=1,intebuf(3)
                  IF(MODUS1.eq.'EECE ') then
                     if(iabs(intebuf_LM(1,ilm)).gt.2*latu_eece(intebuf(1))) cycle
                  endif
                  LL=intebuf_LM(1,ILM)                                       
                  MM=intebuf_LM(2,ILM)       
                  WRITE(28,2011) LL,MM                                 
                  WRITE(28,2022) ( vRHObuf(I,ILM), I=1,intebuf(4) )
                  WRITE(28,2031)                                       
               end DO
               WRITE(28,2030)                                     
               endif
               deallocate(intebuf_lm,rhobuf,vrhobuf)

            enddo
         else

            allocate(intebuf_lm(2,LMMAX),RHObuf(1:IMAX,1:LMMAX),vRHObuf(1:IMAX,1:LMMAX))

            intebuf_lm(1:2,1:LMMAX)=lm(1:2,1:LMMAX)
            RHObuf(1:IMAX,1:LMMAX)=RHOLM(1:IMAX,1:LMMAX)
            vRHObuf(1:IMAX,1:LMMAX)=vRHOLM(1:IMAX,1:LMMAX)

            call mpi_send( intebuf, 4, mpi_integer, 0, 0, mpi_vec_comm, ierr)
            call mpi_send( intebuf_lm, 2*intebuf(3), mpi_integer, 0, 0, mpi_vec_comm, ierr)                                
            if(vresp_write) then
            call mpi_send( vRHObuf, intebuf(4)*intebuf(3), mpi_double_precision, 0, 0, mpi_vec_comm, ierr)
            endif
            if(MODUS.EQ.'TOT  '.or.MODUS.EQ.'FOR  '.or.MODUS.EQ.'CLM  ')then
               call mpi_send( RHObuf, intebuf(4)*intebuf(3), mpi_double_precision, 0, 0, mpi_vec_comm, ierr)
            endif

            deallocate(intebuf_lm,rhobuf,vrhobuf)

         endif

         deallocate(intebuf)


#endif
         call cputim(time4)
         call walltim(time4_w)
         time_writeclm=time_writeclm+time4-time3
         time_writeclm_w=time_writeclm_w+time4_w-time3_w
!
!bk   93/09/13: calculate integral of Veff*grad(rholm)
!     + working with real spherical y_lmp, p=+/-
 
         if (force.and.myid_vec.eq.0) call FOMAI2(jatom,nr,lmmax,fvdrho)

         IF (MODUS.EQ.'EFG   '.and.lqmax.ne.0.and.myid_vec.eq.0) CALL efgsplit(jatom,lqmax,lefg,mefg,nefg,rhoefg,r)          

         IF(LMMAX.GT.1) THEN 
            DO 7370 ILM1=1,LMMAX                                            
               IF (LM(1,ILM1).EQ.2.AND.LM(2,ILM1).EQ.0) THEN 
                  DO I=1,IMAX                                                
                     RHOLM(I,ILM1)=RHOLM(I,ILM1)/R(I)**3
!                     write(78,'(2f13.7)') r(i),rholm(i,ilm1)
                  enddo
                  LABEL=1 
                  DO 7375 ILM2=ILM1+1,LMMAX
                     IF (LM(1,ILM2).EQ.2.AND.LM(2,ILM2).EQ.2) THEN
                        DO I=1,IMAX
                           RHOLM(I,ILM2)=RHOLM(I,ILM2)/R(I)**3
                        enddo
                        LABEL=2
                        DO 7377 ILM3=ILM2+1,LMMAX
                          IF (LM(1,ILM3).EQ.-2.AND.LM(2,ILM3).EQ.2) THEN
                              DO I=1,IMAX
                                 RHOLM(I,ILM3)=RHOLM(I,ILM3)/R(I)**3
                              enddo
                              LABEL=3
                              GOTO 7372
                           ENDIF
 7377                   CONTINUE
                        GOTO 7372
                     ENDIF
 7375             CONTINUE
                  GOTO 7372
               ENDIF
 7370       CONTINUE                                                        
            GOTO 7371                                                         
 7372       CONTINUE 
!            DO 7374 I1=0,100,20
            EFGFACT=0.02997925
            FAC=0.8D0*PI*324.14D0
            DO 7374 I1=0,jri(jatom)-5,1000
               CALL CHARGEL2(R,1.D0,1.D0,RHOLM(1,ILM1),DX(JATOM) &
                    ,JRI(JATOM)-I1,TC) 
!               CALL rel_CHARGE(R,1.D0,1.D0,RHOLM(1,ILM1),DX(JATOM) &
!                    ,JRI(JATOM)-I1,rTC,zz(jatom)) 
!               CALL rel_CHARGE2(R,1.D0,1.D0,RHOLM(1,ILM1),DX(JATOM) &
!                    ,JRI(JATOM)-I1,r2TC,zz(jatom)) 
               TC=TC*FAC
               EFG20=TC
!               efg20r=rtc*fac
!               efg20r2=r2tc*fac
               IF (LABEL.GT.1) THEN
                  CALL CHARGEL2(R,1.D0,1.D0,RHOLM(1,ILM2),DX(JATOM) &
                       ,JRI(JATOM)-I1,TC) 
                  TC=TC*FAC
                  EFG22=TC
                  IF (LABEL.EQ.3) THEN
                     CALL CHARGEL2(R,1.D0,1.D0,RHOLM(1,ILM3),DX(JATOM) &
                          ,JRI(JATOM)-I1,TC)
                     TC=TC*FAC
                     EFG2M=TC
                     QXX=(-EFG20/SQRT(3.D0)+EFG22)*SQRT(15.D0/4.D0/PI)
                     QYY=(-EFG20/SQRT(3.D0)-EFG22)*SQRT(15.D0/4.D0/PI)
                     QZZ=(EFG20*2.D0/SQRT(3.D0))*SQRT(15.D0/4.D0/PI)
                     QXY=EFG2M*SQRT(15.D0/4.D0/PI)
!       EFG in 10**21 V/(m*m)
                     QXX=-QXX*EFGFACT
                     QXY=-QXY*EFGFACT
                     QYY=-QYY*EFGFACT
                     QZZ=-QZZ*EFGFACT
!         QP=.5*(QXX+QYY)
!         QM=.5*(QXX-QYY)
!         VXX=QP+SQRT(QM*QM+QXY*QXY)      
!         VYY=QP-SQRT(QM*QM+QXY*QXY)
!	  VZZ=QZZ 
                     IF (I1.EQ.0) WRITE(6,3019)
                     IF (I1.EQ.0) WRITE(21,3019)
                     WRITE(6,3121) jatom,QXX,QXY,QYY,QZZ, &
                          R(JRI(JATOM)-I1)
                     IF (I1.EQ.0) WRITE(21,3121) jatom,QXX,QXY,QYY,QZZ, &
                          R(JRI(JATOM)-I1)      
                  ELSE
                     VXX=(-EFG20/SQRT(3.D0)+EFG22)*SQRT(15.D0/4.D0/PI)
                     VYY=(-EFG20/SQRT(3.D0)-EFG22)*SQRT(15.D0/4.D0/PI)
                     VZZ=(EFG20*2.D0/SQRT(3.D0))*SQRT(15.D0/4.D0/PI)
                     IF (I1.EQ.0) WRITE(6,3017)
!C       EFG in 10**21 V/(m*m)
                     VXX=-VXX*EFGFACT
                     VYY=-VYY*EFGFACT
                     VZZ=-VZZ*EFGFACT
                     IF (I1.EQ.0) WRITE(21,3017)
                     WRITE(6,3021) jatom,VXX,VYY,VZZ,R(JRI(JATOM)-I1)
                     IF (I1.EQ.0) WRITE(21,3021) jatom,VXX,VYY,VZZ &
                          ,R(JRI(JATOM)-I1)      
                  ENDIF
               ELSE
                  EFG20=EFG20*2.D0/SQRT(3.D0)*SQRT(15.D0/4.D0/PI)
!C       EFG in 10**21 V/(m*m)
                  EFG20=-EFG20*EFGFACT
                  WRITE(6,7207)  JATOM,EFG20,R(JRI(JATOM)-I1)
!                  WRITE(6,7208)  JATOM,EFG20r,R(JRI(JATOM)-I1)
!                  WRITE(6,7209)  JATOM,EFG20r2,R(JRI(JATOM)-I1)
                  IF (I1.EQ.0) WRITE(21,7720) jatom,JATOM,EFG20, &
                       R(JRI(JATOM)-I1)
               ENDIF
 7374       CONTINUE                                                          
         ENDIF                                                             
 7207    FORMAT(1H ,' EFG INSIDE SPHERE',I5,1H:,F12.6,5X,' UP TO R = ',     &
              F10.7)                                                          
!!$ 7208    FORMAT(1H ,' EFGr INSIDE SPHERE',I5,1H:,F12.6,5X,' UP TO R = ',   &
!!$              F10.5)                                                       
!!$ 7209    FORMAT(1H ,' EFGr2 INSIDE SPHERE',I5,1H:,F12.6,5X,' UP TO R = ',  &
!!$              F10.5)                                                       
 7720    FORMAT(':VZZ',i3.3,':',1X,'EFG INSIDE SPHERE ',I3,' = ', &
              F12.6,5X,' UP TO R =',F10.5)
 7371    CONTINUE                                                          
         REWIND 10

         CALL cputim(time3) 
         call walltim(time3_w)
         time_writeefg=time_writeefg+time3-time4
         time_writeefg_w=time_writeefg_w+time3_w-time4_w
 
#ifdef Parallel
!fix for scf and dmat files for all atoms
        if (myid_atm.ne.0) then
            rewind 21
            iline=0
            do
               read(21,'(a)',IOSTAT=ios) line
               if (ios.ne.0) exit
               iline=iline+1
            enddo
            allocate(lines(iline))
            rewind 21
            do i=1,iline
               read(21,'(a)') line
               lines(i)=line
            enddo
           
            jatom_21=jatom
            close(21)
!dmat
 
#endif


      ENDDO non_equiv_loop


#ifdef Parallel

     CALL MPI_BCAST(nnlo,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)

     if (force.and.myid_vec.eq.0) then

        fvdrho_buf(0:3,ndif)=0.0d0
        fsph_buf(0:3,ndif)=0.0d0
        fnsp_buf(0:3,ndif)=0.0d0
        fsph2_buf(0:3,ndif)=0.0d0
        forb_buf(0:3,ndif)=0.0d0    
#endif


      WRITE(6,881)  XWT                                                 
      WRITE(21,881)  XWT                                                
      IF(MODUS1.ne.'EECE ') then
        WRITE(21,725) ETOT + ts2                                           
        WRITE(6,882) ETOT  + ts2     
      endif                                          
!                                                                       
!******************************************************************* *PH*
!....WRITE Output for the BR1 matrix and the ROTIJ Matrix in file case.rotlm:
!     they are read by LectureROTIJ in lecture.f in SRC_elnes/ or SRC_mdff/.
!    These lines are in SRC_lapw2lm/l2main.frc: surch *PH*
!    The number of the output file is 22. his extention must be .rotlm
      IF(MODUS.EQ.'QTL ') then
         DO JR=1, 3
            WRITE(922, 101) BR1(1,JR), BR1(2,JR), BR1(3,JR)
         ENDDO
         INDEX = 0
         DO LATOM=1, NAT
            DO LATEQ=1, MULT(LATOM)
               INDEX = INDEX + 1
               WRITE(922, 100) LATOM, LATEQ, INDEX
               DO JR=1, 3
                  WRITE(922, 101) (ROTIJ(JC,JR,INDEX),JC=1,3)
               ENDDO
            ENDDO
         ENDDO
 100   FORMAT('inequivalent atomnumber ',I3,' number ',I2,' total ',I4)
 101   FORMAT(3F10.5)
      endif
!********************************************************************* *PH*

!     output of forces
      if (force) call FOMAI3(forout,forcea,ftot,fvdrho,fsph,fnsp, &
                        fpsi,fequ,fsph2,forb)

      REWIND 10                                                         
      IF(TEST1.GT.15.d0) THEN                                              
       IF(myid.eq.0) WRITE(6,1001)  TEST1,EQBAD, jatombad,lbad          
       IF(myid.eq.0) WRITE(21,1001) TEST1,EQBAD, jatombad,lbad             
      else IF(TEST1.GT.2.d0) THEN                                              
       IF(myid.eq.0) WRITE(6,1002)  TEST1,EQBAD, jatombad,lbad             
       IF(myid.eq.0) WRITE(21,1002) TEST1,EQBAD, jatombad,lbad               
      ENDIF                                                          
777  CONTINUE

      if (force) then
         deallocate(h_ablyl_hk,aalm_buf_tmp)
      endif

      CALL fini_charp
      CALL fini_chard
      CALL fini_charf
      CALL fini_lohelp
      CALL fini_xa
      CALL fini_xa3

      IF(MODUS.NE.'QTL  '.AND.TEST1.LT.2.0D0) then
        CALL fini_xdos
        RETURN
      else IF(MODUS.NE.'QTL  '.AND.TEST1.LT.15.0D0) then     
        IF(myid.eq.0) write(21,9001) test1,EQBAD, jatombad,lbad
        IF(myid.eq.0) write(21,9003) 
        CALL fini_xdos
        RETURN
      else if(MODUS.eq.'QTL  ') then                    
        if (myid.eq.0) CALL OUTP(eferm,so,nspin1,fnamehelp,doos)
        return
      else                    
        IF(myid.eq.0) write(21,9001) test1,EQBAD, jatombad,lbad
        IF(myid.eq.0) write(21,9003) 
        CALL fini_xdos
        if(modus.eq.'ALM  ') return
        CALL OUTERR('l2main','QTL-B.GT.15., Ghostbands, check scf files')

        STOP 'L2main - QTL-B Error'
      endif

 9001   format(':WARN : QTL-B value eq.',f7.2,' in Band of energy ',F9.5,&
 2      x,'ATOM=',i5,2x,'L=',i3)
 9003   format(':WARN : You should change the E-parameter for this atom and L-value in case.in1 (or try the -in1new switch)')
!
!        Error messages
!
 900  CALL OUTERR('l2for','ja.ne.jatom.')
      call abort_parallel
      STOP 'L2for - Error'
 910  CALL OUTERR('l2for','ngau too small.')
      call abort_parallel
      STOP 'L2for - Error'
! 920  CALL OUTERR('l2main','NEFG.GT.23.')
!      STOP 'L2main - Error'
 950  CALL OUTERR('l2main','error reading parallel vectors')
      call abort_parallel
      STOP 'L2main - Error'
!
 205  FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I5,1X,A10)                      
!     
 207  FORMAT(1H0,' TOTAL CHARGE INSIDE SPHERE',I5,1H:,F12.6)            
 250  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,PZ,PXY',/,':QTL',i3.3,':',12F7.4)               
 251  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,      ','D-EG,D-T2G ',/,':QTL',i3.3,':',12F7.4)  
 252  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,      ','D-Z2,D-XY,X2Y2,D-XZ,YZ ',/, &
           ':QTL',i3.3,':',12F7.4)
 253  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,PZ,PXY,,','D-Z2,D-XY,D-X2Y2,D-XZ,YZ ',/, &
           ':QTL',i3.3,':',12F7.4)                          
 254  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,PZ,PXY, ','D-Z2,D-XY,X2Y2,D-XZ,YZ ',/, &
           ':QTL',i3.3,':',12F7.4)                            
 255  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,      ','D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ ',/, &
           ':QTL',i3.3,':',12F7.4)                        
 256  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,PX,PY,PZ,',' ',/, &
           ':QTL',i3.3,':',12F7.4)
 257  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,PX,PY,PZ,','D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ ',/, &
           ':QTL',i3.3,':',12F7.4)                        
 258  FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3, &
           ' S,P,D,F,PX,PY,PZ,D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ,F00,F11,F22', &
           ',F33,','F1M,F2M,F3M ',/,':QTL',i3.3,':',19F7.4)
 260  FORMAT(8x,'Q-s-low E-s-low   Q-p-low E-p-low   Q-d-low E-d-low', &
          '   Q-f-low E-f-low')
 261  FORMAT(':EPL',i3.3,':',4(2F8.4,2x))
 262  FORMAT(8x,'Q-s-hi  E-s-hi    Q-p-hi  E-p-hi    Q-d-hi  E-d-hi ', &
          '   Q-f-hi  E-f-hi ')
 263  FORMAT(':EPH',i3.3,':',4(2F8.4,2x))
501   FORMAT(2i5,' k-point, nband')
503   FORMAT(5e17.9)
 720  FORMAT(/,':CHA',i3.3,':',1X,'TOTAL VALENCE CHARGE INSIDE SPHERE ',I3, &
           ' = ',F8.4,4x,'(RMT=',f8.4,' )')         
 725  FORMAT(/,':SUM  :',1X,'SUM OF EIGENVALUES =  ',F20.9,/)            
 787  FORMAT('    VALENCE CHARGE DENSITY   IN  MT SPHERES',5X,I3,        &
           ' ITERATION')
 800  FORMAT(/,10X,68(1H-),/,11X,'W A V E F U N C T I O N S',        &
           '   AND   C H A R G E S   IN   S P H E R E S',/,            &
           10X,68(1H-))                                         
 881  FORMAT(/,':CHA  :',' TOTAL VALENCE CHARGE INSIDE UNIT CELL =',F15.6) 
 882  FORMAT(//,3X,'SUM OF EIGENVALUES:'/24X,F15.6/)                    
!QDMFT              
 8882  FORMAT(//,3X,'KINETIC ENERGY with DMFT correction:',16X,F15.6/)
 8883  FORMAT(//,3X,'HUBBARD ENERGY(included in SUM OF EIGENVALUES):',5X,F15.6/)
!END QDMFT
 1001 FORMAT(//,'   QTL-B VALUE .EQ. ',F10.5,' in Band of energy ',F9.5,&
                2x,'ATOM=',i5,2x,'L=',i3,/, &
                '    Check for ghostbands or EIGENVALUES BELOW XX messages',/, &
                '    Adjust your Energy-parameters for this ATOM and L (or use -in1new switch), check RMTs  !!!',//)
 1002 FORMAT(//,'   QTL-B VALUE .EQ. ',F10.5,'   in Band of energy ',F10.5,&
                 3x,'ATOM=',i5,3x,'L=',i3,/, &
                '    Most likely no ghostbands, but adjust Energy-parameters for this ATOM and L',//)        
 1003 FORMAT(251(I3,I2))                                                 
 1005 FORMAT(/,7X,'LMMAX',I3,/,7x,'LM= ',17(i3,I2),(/,7x,18(i3,i2)))          
 1006 FORMAT(/,7X,'LMMAX',I3,/)          
 1990 FORMAT(3X,'ATOMNUMBER   ',I3,5X,10A4)                             
 2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
 2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
 2022 FORMAT(3X,4ES19.12)                                                 
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
 2032 FORMAT(49X,I3,//)  
2055  FORMAT(/,1X,' K-POINT:',3F14.10,1X,I5,I4,2X,A10)
 3017 FORMAT(/,22X,'VXX',9X,'VYY',9X,'VZZ',7X,'UP TO R')
 3018 FORMAT(/'  EFG COMPONENTS BEFORE TRANSFORMATION [10**21 V/(M*M)]:' &
           ,/, &
           15X,'QXX (-V20/2. + V22)     =',F10.3,/, &
           15X,'QYY (-V20/2. - V22)     =',F10.3,/, &
           15X,'QZZ ( V20)              =',F10.3,/, &
           15X,'QXY ( V2M/2. * SQRT(3.))=',F10.3)  
 3019 FORMAT(/,22X,'QXX',9X,'QXY',9X,'QYY',9X,'QZZ',  &
           7X,'UP TO R')
 3020 FORMAT(/'  EFG COMPONENTS AFTER TRANSFORMATION [10**21 V/(M*M)]:' &
           ,/, &
           15X,'VXX                     =',F10.3,/, &
           15X,'VYY                     =',F10.3,/, &
           15X,'VZZ                     =',F10.3,/)
 3121 FORMAT(':VZZ',i3.3,':',8x,4F12.5,F12.3)
 3021 FORMAT(':VZZ',i3.3,':',8x,3F12.5,F12.3)
 3022 FORMAT(  15X,'QXX=',F10.3,/,15X,'QYY=',F10.3,/,15X,'QZZ=',F10.3 &
           ,15X,'QXY=',F10.3,' UP TO R= ',F10.5)
4645  FORMAT(2i4,3f15.10)
4646  FORMAT(10e20.10)
4647  FORMAT(4e20.10)
4893  FORMAT(3i4,5(2e16.8,2x))  
   13 FORMAT(////,':POS',i3.3,':',1x,'ATOM ',I4,1X,'X,Y,Z =', &
             3F8.5,2X,'MULT=',I2,'  ZZ=',f7.3,2x,a10)

      END                                                               
