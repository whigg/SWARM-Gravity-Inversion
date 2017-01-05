        Program L1Cmpi
        implicit none

        INCLUDE 'mpif.h'
        INTEGER err, rank, size
        integer istart,iend,rankp

        integer year,month,day
        integer years,months,days,yeare,monthe,daye
        real*8 hr,TJDs,TJDe,TJD
        integer nd,i5
        character yearc*4,monthc*4,dayc*4
        integer errsys, nzs2, nres, n10
        
        CALL MPI_INIT(err)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err)
        
        hr=0d0
        years=2002;months=4;days=1
        yeare=2014;monthe=5;daye=31
!        years=2014;months=4;days=1
!        yeare=2014;monthe=4;daye=30
        call JULDAT (years,months,days,hr,TJDs)
        call JULDAT (yeare,monthe,daye,hr,TJDe)
        
        nd=anint(TJDe-TJDs+1)

!        rankp=INT(nd/dble(size))+1
!        istart=rank*rankp+1
!        iend=istart+rankp-1
!        if(iend.gt.nd)iend=nd

        nzs2=nd
        rankp=INT(nzs2/dble(size))
        nres=nzs2-size*rankp !residual
        if(rank==0)then
          istart=1
        else
          istart=rank*rankp+min(rank-1,nres)+1
        endif
        if(1.le.rank.and.rank.le.nres)then
          n10=1
        else
          n10=0
        endif
        iend=istart+rankp+n10-1
        if(rank == size-1) then
                iend=nzs2
        endif
!        write(*,*)'rank,istart,iend=',rank,istart,iend
!        write(*,*)'nd, rankp=',nd,rankp


        do i5=istart,iend
!       nl=iend-istart+1
!       il=i5-istart+1
          TJD=TJDs+i5-1d0
          call CALDAT (TJD,Year,Month,day,Hr)
          write(yearc,'(I4)')year
          write(monthc,'(I4)')month
          write(dayc,'(I4)')day

!        write(*,*)'rank=',rank,istart,iend,i5,nd,rankp,yearc,monthc,dayc
        write(*,'(I4,I3,I3,8I5)') year,month,day,
     .   rank,istart,iend,i5,nd,rankp, iend-istart+1,iend-i5+1

!       call system('date',errsys) !SIGSEGV, segmentation fault occurred
!        call system('./testd1.scr '//yearc//' '//monthc//' '//dayc)
!       call system('cd ../ && pwd && ls')
!        call system('./L1B2L1C_day_woFIT.scr '
        call system('./bin/l1b2l1c.scr '
     .               //yearc//' '//monthc//' '//dayc//' 1')

        enddo !i5

        CALL MPI_FINALIZE(err)        
        stop

        end


         SUBROUTINE JULDAT (I,M,K,H,TJD)
*
*     THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND
*     TIME.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT TIME VALUE
*     CAN BE IN ANY UT-LIKE TIME SCALE (UTC, UT1, TT, ETC.) - OUTPUT
*     JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
*     VAN FLANDERN.
*
*          I      = YEAR (IN)
*          M      = MONTH NUMBER (IN)
*          K      = DAY OF MONTH (IN)
*          H      = UT HOURS (IN)
*          TJD    = JULIAN DATE (OUT)
*
*
      DOUBLE PRECISION H,TJD

*     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
      JD = K-32075+1461*(I+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12
     .     -3*((I+4900+(M-14)/12)/100)/4
      TJD = JD - 0.5D0 + H/24.D0

      RETURN
      END



      SUBROUTINE CALDAT (TJD,I,M,K,H)
*
*     THIS SUBROUTINE COMPUTES CALENDAR DATE AND TIME, GIVEN JULIAN
*     DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE TIME SCALE
*     (UTC, UT1, TT, ETC.) - OUTPUT TIME VALUE WILL HAVE SAME BASIS.
*     OUTPUT CALENDAR DATE WILL BE GREGORIAN.  ALGORITHM BY FLIEGEL AND
*     VAN FLANDERN.
*
*          TJD    = JULIAN DATE (IN)
*          I      = YEAR (OUT)
*          M      = MONTH NUMBER (OUT)
*          K      = DAY OF MONTH (OUT)
*          H      = UT HOURS (OUT)
*
*
      DOUBLE PRECISION TJD,H,DJD,DMOD

      DJD = TJD + 0.5D0
      JD = DJD
      H = DMOD (DJD,1.D0) * 24.D0
*     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
      L = JD + 68569
      N = 4*L/146097
      L = L - (146097*N+3)/4
*     I=YEAR, M=MONTH, K=DAY
      I = 4000*(L+1)/1461001
      L = L - 1461*I/4 + 31
      M = 80*L/2447
      K = L - 2447*M/80
      L = M / 11
      M = M + 2 - 12*L
      I = 100*(N-49) + I + L

      RETURN
      END
