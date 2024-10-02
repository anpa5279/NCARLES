      subroutine print(lu,it,iz_strt,iz_end)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      write(lu,4000)
 4000 format(30X,' --- SOLUTION ---')
      write(lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1),
     +              zol,hol,ucfl, vcfl, wcfl,
     +              t_ref, x_hurr, y_hurr, u_10, v_10, cd_10,
     +              tau_x, tau_y, turb_la, cpou10,
     +              stokess, stokesw, fcor, fcor_h
 4100 format(' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/,
     +       10x,'ZTOP = ',e15.6,
     +       ',  TSFC = ',e15.6,
     +       ',  UW = ',e15.6,',  VW = ',e15.6,/,10x,
     +       'WT = ',e15.6,',  ZL =',e15.6,
     +       ',  HL = ',e15.6,/,10x,'U_cfl = ',e15.6,
     +       ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,/,10x,
     +       'Theta Ref = ',e15.6,/,10x,
     +       'X hurr = ',e15.6,',  Y hurr = ',e15.6,/,10x,
     +       'U_10 = ',e15.6,',  V_10 = ',e15.6,
     +       ',  Cd_10 = ',e15.6,/,10x,
     +       'Tau_x = ',e15.6,',  Tau_y = ',e15.6,/,10x,
     +       'La_t = ',e15.6,',  Wave age = ',e15.6,/,10x,
     +       'Stokes amp = ',e15.6,',  Stokes wave number = ',e15.6,/,
     +       10x,'f cor = ',e15.6,',  f cor_h = ',e15.6)
      write(lu,4200)
 4200 format(//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ',
     +       //,2x,'IZ',4x,'T_MEAN',7x,
     + 'DIVG',8X,'LE_KE',6X,'SGS_KE',7X,'LE_WT',6X,
     + 'SGS_WT',7X,'SHRZ',8X,'BUOY')
      do 19999 iz=iz_end,iz_strt,-1
         write(lu,4300)iz,txym(iz,1)-t_ref,divz(iz),
     +              englez(iz),eavg(iz),wtle(iz,1),
     +              wtsb(iz,1),shrz(iz),buyz(iz)
 4300    format(1X,I3,e12.4,7e12.4)
19999 continue
      write(lu,4400)tsfcc(1),wtsfc(1)
 4400 format('  SURFACE VALUE: TXYM=',F8.2,'               WTSB=',E9.2)

         write(lu,4500) stokess,udrift,vdrift
 4500    format(/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4,
     +          ' VDRIFT = ',e12.4)

      write(lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz),
     +       uwsb(iz),vwle(iz),vwsb(iz),iz=iz_strt,iz_end)
 4600 format(//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',10x,
     +          ' UWSB',10x,' VWLE',10x,' VWSB'
     +       ,/,(1x,i4,6(3x,e13.6)))
      if(ivis .eq. 1) then
         write(lu,4800) xksurf, nmatch, viscon, vise
 4800    format(//,' XKSURF = ',e15.6,' NMATCH = ',i4,/,
     +             ' VISCON = ',e15.6,' VISE = ',e15.6)
         write(lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
 4700    format(//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
      endif
      write(lu,5100) (iz,txym(iz,2),txym(iz,3),
     +              txym(iz,4),txym(iz,5),txym(iz,6),
     +              txym(iz,7),iz=iz_strt,iz_end)
 5100 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,'2-MEAN',10x,
     +          '3-MEAN',10x, '4-MEAN', 10x, '5-MEAN',10x,
     +          '6-MEAN',10x 
     +       ,/,(1x,i4,6(3x,e13.6)))
      write(lu,5200) (iz,txym(iz,8),txym(iz,9),txym(iz,10),txym(iz,11),txym(iz,12),iz=iz_strt,iz_end)
 5200 format(//,' IZ',5x,' SCALAR-7 MEAN',8x,'8-MEAN',10x,
     +          '9-MEAN',10x, '10-MEAN', 10x, '11-MEAN', 10x
     +       ,/,(1x,i4,5(3x,e13.6)))


      return
      end
