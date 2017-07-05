	subroutine bsbn(tempe,dense)

C       AAK (8/12/09): This is very similar to the program given in
C       appendix E.1 of Gordon and Sorochenko (2002). Radio
C       Recombination lines: their physics and applications.

C     GENERAL BN PROGRAM.  RADIO RECOMBINATION LINES FROM H REGIONS AND   
C 1   COLD INTERSTELLAR CLOUDS: COMPUTATION OF THE BN FACTORS.            
C 2   M. BROCKLEHURST, M. SALEM.                                          
C REF. IN COMP. PHYS. COMMUN. 13 (1977) 39                                
C JOB MS13 2355 BN PROGRAM TEST                                           
C ROUTE PRINTER WESTCAM,POST WCAV,NOTIFY                                  
C MSGLEVEL=1                                                              
C LIMSTORE 175K,COMP 75 SECS,PRINTER 5000                                 
C //ONE EXEC FTG1CLG,REGG=175K,LISTC=SOURCE,MAPL=MAP  COMPILE, LOAD, GO   
C //FORT.SYSIN DD *  SOURCE PROGRAM CARDS                                 
C                      WIDE TEMPERATURE RANGE BN PROGRAM                
C                      *********************************                
C                                                                       
C  PROGRAM FOR CALCULATION OF COEFFICIENTS OF DEPARTURE FROM            
C  THERMODYNAMIC EQUILIBRIUM, BN, FOR HYDROGENIC ATOMS, AT ELECTRON     
C  TEMPERATURES FROM 10K TO 20 000K                                     
C                                                                       
C  INPUT CARDS -                                                        
C  ** ALPHANUMERIC TITLE CARD.                                          
C  ** NORMALLY, THE FOLLOWING 6 CARDS (SEE WRITEUP FOR EXPLANATION) -   
C 75  2  4                                                              
C 30 31 32 33 34 35 37 39 41 43 46 49 52 55 58 61 64 68 72 76 80 84 88  
C 92 97102107112117122127132138144150156162168174180187194201208215222  
C230238246254262270279288297306315325335345355365375386397408419430441  
C452463474485496507                                                     
C 75 72 69 66                                                           
C                                                                       
C  ** ONE CARD WITH RADIATION TEMPERATURE AND EMISSION MEASURE OF       
C  BACKGROUND RADIATION FIELD.
C  FORMAT 2E10.3. IF EMISSION MEASURE READ  
C  IS GE 10**10, IT IS TAKEN TO BE INFINITE. THIS CARD IS READ BY       
C  FUNCTION COR(N,ISW). IF THIS SUBPROGRAM IS REPLACED BY THE USER, ANY 
C  CARDS READ BY  COR  WHEN CALLED WITH ISW=0 SHOULD BE PLACED HERE.    
C  (BLANK CARD = NO FIELD).                                             
C                                                                       
C  ** (ONE CARD FOR EACH CASE TO BE CALCULATED) TEMPERATURE, DENSITY,   
C  CASE (THIN=A - 1, THICK=B - 2), PRINT CYCLE (PRINT EVERY K-TH LEVEL),
C  NPLO AND NPHI, WHERE OUTPUT IS TO BE PUNCHED FROM N=NPLO TO NPHI     
C (NO PUNCHED OUTPUT IF NPHI=0), ALPHANUMERIC LABEL TO BE PUNCHED IN    
C  COLUMNS 77-79 OF CARD OUTPUT.                                        
C  FORMAT 2E10.5, 4I5, 37X, A3 (I.E., LABEL IN COLS. 78-80)             
C                                                                       
C  ** A BLANK CARD (WHICH ENDS EXECUTION)                               
C                                                                      
C                                                                       
	common/ncont/nplo,nphi
      COMMON /EXPDAT/ CXP(707),MAXN                                     
      COMMON /FITDAT/ AFIT(4,4),IVAL(4),NFIT                            
      COMMON /PARMS/ DENS,T,ITM                                         
      COMMON /TDEP/ TE32,TE12,CTE                                       
      DIMENSION CO(75), DVAL(507), IND(75,2), IPIV(75), KBOUT(507),     
     1 KCOUT(507), MVAL(75), SK(75,75) , VAL(507)             
C                                                                       
C  MANY MACHINES DO NOT REQUIRE THE FOLLOWING  DOUBLE PRECISION         
C  STATEMENT IN MAIN OR SUBPROGRAMS                                     
C                                                                       
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
C                                                                       
      DOUBLE PRECISION AFIT,ARG,CO,COR,CTE,CX,CXP,D,DENS,DVAL,H,HH,RATIO
     1,SK,T,T1,TE12,TE32,VAL,X,XXXI                               
C     Added by AAK in an attempt to get code to work
      DOUBLE PRECISION TEMPE, DENSE

      DATA NONE,NBL/1H0,1H /,LPPG/45/                                   

	open(unit=3,file='bsbn.inp',form='formatted',status='old')
	open(unit=10,file='bsbn.out',form='formatted',status='unknown')

c     first block of bsbn.inp
      READ (3,*) IC,IR,NFIT  
c     second block of bsbn.inp                                      
      READ (3,*) (MVAL(I),I=1,IC)                                  
c     third block of bsbn.inp
      READ (3,*) (IVAL(I),I=1,NFIT)                                

      T1=0.D0                                                           
      MAXN=MVAL(IC)                                                     
C       The COR function reads in the next two values in bsbn.inp which
C       are the tbck and ebck. if either is set to zero then it doesn't
C       calculate the correction to the rate of population of level n
C       due to a radiation field."
      H=COR(0,0)                                                        
      H=COLRAT(0,0,0.D0,0.D0)                                           
      t=tempe
      dens=dense

C last line of bsbn.in
 10   READ (3,*,END=65) NMIN,ICYC,NPLO,NPHI,LABEL

	close(3)
c     Write the input parameters to the output file
	write(10,*)'t_e    n_e       nmin icyc nplo nphi'
        write(10,*)t,dens,NMIN,ICYC,NPLO,NPHI
	write(10,150)

      IF (T.LE.0.D0.OR.DENS.LE.0.D0) STOP                               
      IF (NMIN.LE.0) NMIN=2                                             
      NPLO=MAX0(NPLO,MVAL(1))                                           
      NPHI=MIN0(NPHI,MVAL(IC))                                          
      ND=NPHI-NPLO+1                                                    
      ICYC=MAX0(1,ICYC)                                                 
      NPAGE=1                                                           
      NLINE=0                                                           
      IF (T.EQ.T1) GO TO 30                                             
      ITM=1                                                             
      IF (T.GE.1000.D0) ITM=3                                           
      TE12=DSQRT(T)                                                     
      TE32=T*TE12                                                       
      CTE=15.778D4/T                                                    
      DO 20 I=1,707                                                     
      CX=0.D0                                                           
      ARG=CTE/DFLOAT(I**2)                                              
      IF (ARG.LE.165.D0) CX=DEXP(-ARG)                                  
 20   CXP(I)=CX 

c     AAK--  COLION(N,IONZ, T, QI) THIS SUBROUTINE COMPUTES THE COLLISIONAL
C       IONIZATION RATE, QI, FROM LEVEL N FOR IONS OF EFFECTIVE CHARGE
C       IONZ AT ELECTRON TEMPERATURE T. WHEN CALLED WITH N=0, COLION
C       COMPUTES AND STORES QUANTITIES WHICH DEPEND ONLY UPON
C       TEMPERATURE AND EFFECTIVE CHARGE. IT IS ASSUMED THAT T AND IONZ
C       REMAIN CONSTANT UNTIL THE NEXT CALL WITH N=0.
c     AAK

       
      CALL COLION (0,1,T,H)      

c       I think the point of COLION called with N=0 is to set up EXPX
c       and CONS from looking at the code. You need to compile with the
c       -noautomatic option in ifort to get this to work.

c     RADCAL -- RADIATIVE CASCADE COEFFICIENTS AND COLLISIONAL RATE
      CALL RADCOL (T,MVAL,IC,NMIN)                                         

C     REDUCE
c                                                                       
C   GIVEN A SET OF INTEGERS                                             
C         M(IT),IT=1,IC,SUCH THAT-                                      
C         1) M(IT+1)=M(IT)+1  FOR IT.LE.IA                              
C         WHERE IA.GE.1  AND                                            
C         2) (M(IT+1) - M(IT)).GT.1 FOR IT.GE.IA,                       
C   AND GIVEN A FUNCTION SUBPROGRAM  BK                                 
C         WHICH CALCULATES THE                                          
C         ELEMENTS OF A LARGE                                           
C         M(IC)*M(IC) MATRIX,                                           
C   THIS SUBROUTINE USES LAGRANGE                                       
C         INTERPOLATION OF ORDER                                        
C         2*(IR+1) TO CALCULATE A                                       
C         SMALLER IC*IC MATRIX SK                                       
C   REQUIRES A FUNCTION SUBPROGRAM                                      
C         PHI                                                           
C   IR MUST BE .LE. (IA-1)         
 30   CALL REDUCE (MVAL,IC,IR,SK)                                       
 
C AAK -- RHS
C                                                                       
C  COMPUTES THE RIGHT HAND SIDE OF EQUATIONS (2.7) OF BROCKLEHURST,     
C  MNRAS 148, 417 (1970).                                               
C  
      CALL RHS (CO,MVAL,IC)                                             

C AAK -- JMD: No fucking clue what this routine does. CALLS MATINV and HELPME

      CALL JMD (SK,CO,MVAL,IC)                                          
C AAK -- MATINV: 
C     SOLVES SIMULTANEOUS EQUATIONS IF L=1                              
C                                                                       
C     INVERTS MATRIX A IF L=0                                           
C                                                                       
C     SOLUTIONS ARE RETURNED IN B                                       
C                                   
      CALL MATINV (SK,IC,CO,1,D,IRROR,75,IPIV,IND)                      

C AAK -- INTERP
C                                                                       
C  COMPUTES SOLUTIONS AT ALL  N  FROM THOSE OBTAINED AT THE CONDENSED   
C  POINTS                                                               
C 
      CALL INTERP (MVAL,CO,VAL,DVAL,IC,2)                               

      J=MVAL(1)                                                         
      K=MVAL(IC)                                                        
      IPUN=0                                                            

C     go from Hnplo to Hnphi in increments of icyc
 	do 60, i=nplo,nphi,icyc
c     DO 60 I=J,K                                                       
      RATIO=DVAL(I)/VAL(I)                                              
      H=T*DFLOAT(I)**3*RATIO/3.158D5                                    
      HH=1.D0-H                                                         
      XXXI=DEXP(15.778D4/(DFLOAT(I*I)*T))*VAL(I)                        
      X=-HH*XXXI*100.D0/T                                               
c     IF (MOD(I,ICYC).NE.0) GO TO 50                                    
      NSKIP=NBL                                                         
      IF (MOD(NLINE,5).EQ.0) NSKIP=NONE                                 
      IF (MOD(NLINE,LPPG).NE.0) GO TO 40                                
      NPAGE=NPAGE+1                                                     
 40   NLINE=NLINE+1                                                     
 	write(10,161)I,VAL(I),HH,DVAL(I),RATIO,H,X
 50   IF (NPHI.EQ.0) GO TO 60                                           
      IF (I.LT.NPLO) GO TO 60                                           
      IF (I.GT.NPHI) GO TO 60                                           
      IPUN=IPUN+1                                                       
      KBOUT(IPUN)=VAL(I)*1.D4                                           
      KCOUT(IPUN)=DLOG10(RATIO)*1.D4                                    
 60   CONTINUE                                                          
      T1=T                                                              
      IF (NPHI.NE.0) CALL BCPCH (KBOUT,KCOUT,T,DENS,NMIN,LABEL,NPLO,NPHI
     1,ND)                                                              
      GO TO 10                                                          
   65 continue
	close(10)
C                                                                       
 11	format(a)
 70   FORMAT (10A8)                                                     
 80   FORMAT (1X,23I3)                                                  
 90   FORMAT (20X,10A8/20X,20(4H****)///)             
 100  FORMAT (7H MVAL (,I3,10H VALUES) -/(1X,24I5))                     
 110  FORMAT (7H0IVAL (,I3,10H VALUES) -/(1X,24I5))                     
 120  FORMAT (5H0IR =,I3)                                               
 130  FORMAT (2G10.3,4I5,37X,A3)                                        
 140  FORMAT (14H1TEMPERATURE =,F6.0,14H K,  DENSITY =,1PG10.3,15H/CM**3
     1,  NMIN =,I3,1H.,51X,4HPAGE,I3)                                   
 150  FORMAT (3H0 N,9X,2HBN,16X,4HBETA,14X,6HDBN/DN,10X,12HD(LN(BN))/DN,
     17X,6H1-BETA,13X,4HZETA)                                           
 160  FORMAT (A1,I3,1P6G18.6)                                           
 161  FORMAT (I3,1P6G18.6)                                           

C	return
      END                                                               
C                                                                       
      SUBROUTINE BCPCH (KB,KC,T,DENS,NMIN,LABEL,NPLO,NPHI,NDIM)         
C                                                                       
C  PUNCHES BN AND CN VALUES IN STANDARD FORMAT                          
C                                                                       
      DIMENSION KB(NDIM), KC(NDIM), K1(19), K2(12)                      
      DOUBLE PRECISION DENL,DENS,T,TL                                   
      EQUIVALENCE (K1(1),K2(1))                                         
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      INTEGER ALPHA(35)                                                 
C  SOME COMPILERS REQUIRE  (ALPHA(I),I=1,35)  IN DATA STATEMENT         
      DATA ALPHA/1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1HA,1HB,1HC,1HD,1HE
     1,1HF,1HG,1HH,1HI,1HJ,1HK,1HL,1HM,1HN,1HO,1HP,1HQ,1HR,1HS,1HT,1HU,1
     2HV,1HW,1HX,1HY,1HZ/                                               
      N=NPHI-NPLO+1                                                     
      TL=DLOG10(T)                                                      
      DENL=DLOG10(DENS)                                                 
      NB=N/19                                                           
      IF (19*NB.NE.N) NB=NB+1                                           
      NC=N/12                                                           
      IF (12*NC.NE.N) NC=NC+1                                           
      IPOINT=0                                                          
      DO 20 I=1,NB                                                      
      DO 10 J=1,19                                                      
      IPOINT=IPOINT+1                                                   
      K=0                                                               
      IF (IPOINT.LE.N) K=KB(IPOINT)                                     
      K1(J)=K                                                           
 10   CONTINUE                                                          
 20   CONTINUE                                                          
      IPOINT=0                                                          
      DO 40 I=1,NC                                                      
      DO 30 J=1,12                                                      
      IPOINT=IPOINT+1                                                   
      K=0                                                               
      IF (IPOINT.LE.N) K=KC(IPOINT)                                     
      K2(J)=K                                                           
 30   CONTINUE                                                          
      NBI=NB+I+1                                                        
 40   CONTINUE                                                          
      RETURN                                                            
C                                                                       
 50   FORMAT (2F10.6,3I5,41X,A3,A1)                                     
 60   FORMAT (19I4,A3,A1)                                               
 70   FORMAT (12I6,4X,A3,A1)                                            
 80   FORMAT (////27H *** OUTPUT PUNCHED ON UNIT,I3)                    
      END                                                               
C                                                                       
      BLOCK DATA                                                        
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      COMMON /GAUNTS/ A1(50),A2(50),A3(50),A4(50),A5(50),A6(50),A7(50),A
     18(50),A9(50),A10(50),A11(50),A14(50),A17(50),A20(50),A25(50),A30(5
     20),A40(50),A50(50),A100(50),A150(50),A225(50),A500(50),IXV(12)    
      COMMON /RCMB/ SV0A(33),SV0B(33),SV0C(33),SV1A(33),SV1B(33),SV1C(33
     1),SV2A(33),SV2B(33),SV2C(33)                                      
      COMMON /GAUSS/ VALUE(12)                                          
      DOUBLE PRECISION VALUE                                            
C                                                                       
C  NOTE - SOME COMPILERS REQUIRE AN IMPLIED DO LOOP WHEN THE VALUES OF  
C  A WHOLE ARRAY ARE SET IN A DATA STATEMENT. THE FOLLOWING STATEMENTS  
C  WILL NEED CHANGING TO (A1(I),I=1,50), ETC. (THE ANS STANDARD FORTRAN 
C  FORM IS EXCESSIVELY LENGTHY FOR LARGE ARRAYS.)                       
C                                                                       
C  RADIATIVE GAUNT FACTORS                                              
C                                                                       
      DATA A1/.7166,.7652,.7799,.7864,.7898,.7918,.7931,.7940,.7946,.795
     11,.7954,.7957,.7959,.7961,.7963,.7964,.7965,.7966,.7966,.7967,.796
     28,.7968,.7968,.7969,.7969,.7969,.7970,.7970,.7970,.7970,.7970,.797
     31,.7971,.7971,.7971,.7971,.7971,.7971,.7971,.7971,.7972,.7972,.797
     42,.7972,.7972,.7972,.7972,.7972,.7972,.7972/                      
      DATA A2/.7566,.8217,.8441,.8549,.8609,.8647,.8672,.8690,.8702,.871
     12,.8720,.8726,.8730,.8734,.8737,.8740,.8742,.8744,.8746,.8747,.874
     29,.8750,.8751,.8751,.8752,.8753,.8753,.8754,.8755,.8755,.8755,.875
     36,.8756,.8756,.8757,.8757,.8757,.8757,.8758,.8758,.8758,.8758,.875
     48,.8759,.8759,.8759,.8759,.8759,.8759,.8759/                      
      DATA A3/.7674,.8391,.8653,.8784,.8861,.8910,.8944,.8968,.8986,.900
     10,.9011,.9019,.9026,.9032,.9037,.9041,.9044,.9047,.9049,.9052,.905
     24,.9055,.9057,.9058,.9059,.9060,.9061,.9062,.9063,.9064,.9064,.906
     35,.9066,.9066,.9067,.9067,.9067,.9068,.9068,.9068,.9069,.9069,.906
     49,.9070,.9070,.9070,.9070,.9070,.9071,.9071/                      
      DATA A4/.7718,.8467,.8750,.8896,.8984,.9041,.9081,.9110,.9132,.914
     19,.9163,.9173,.9182,.9190,.9196,.9201,.9205,.9209,.9213,.9215,.921
     28,.9220,.9222,.9224,.9226,.9227,.9228,.9230,.9231,.9232,.9233,.923
     33,.9234,.9235,.9235,.9236,.9237,.9237,.9238,.9238,.9238,.9239,.923
     49,.9240,.9240,.9240,.9240,.9241,.9241,.9241/                      
      DATA A5/.7741,.8507,.8804,.8960,.9055,.9118,.9162,.9195,.9220,.924
     10,.9255,.9268,.9278,.9287,.9294,.9300,.9306,.9310,.9314,.9318,.932
     21,.9324,.9326,.9329,.9331,.9332,.9334,.9335,.9337,.9338,.9339,.934
     30,.9341,.9342,.9343,.9344,.9344,.9345,.9345,.9346,.9347,.9347,.934
     48,.9348,.9348,.9349,.9349,.9349,.9350,.9350/                      
      DATA A6/.7753,.8531,.8837,.8999,.9099,.9167,.9215,.9251,.9278,.930
     10,.9317,.9331,.9343,.9352,.9361,.9368,.9374,.9379,.9384,.9388,.939
     22,.9395,.9398,.9400,.9403,.9405,.9407,.9408,.9410,.9412,.9413,.941
     34,.9415,.9416,.9417,.9418,.9419,.9420,.9420,.9421,.9422,.9422,.942
     43,.9423,.9424,.9424,.9425,.9425,.9426,.9426/                      
      DATA A7/.7761,.8547,.8858,.9025,.9130,.9200,.9251,.9289,.9318,.934
     12,.9360,.9376,.9389,.9399,.9408,.9416,.9423,.9429,.9434,.9439,.944
     23,.9447,.9450,.9453,.9455,.9458,.9460,.9462,.9464,.9466,.9467,.946
     38,.9470,.9471,.9472,.9473,.9474,.9475,.9476,.9477,.9477,.9478,.947
     49,.9479,.9480,.9480,.9481,.9481,.9482,.9482/                      
      DATA A8/.7767,.8558,.8873,.9044,.9151,.9224,.9277,.9317,.9348,.937
     13,.9393,.9409,.9423,.9434,.9444,.9453,.9460,.9467,.9472,.9477,.948
     22,.9486,.9489,.9493,.9496,.9498,.9501,.9503,.9505,.9507,.9509,.951
     30,.9512,.9513,.9514,.9515,.9517,.9518,.9519,.9519,.9520,.9521,.952
     42,.9522,.9523,.9524,.9524,.9525,.9525,.9526/                      
      DATA A9/.7771,.8565,.8884,.9058,.9167,.9242,.9297,.9338,.9370,.939
     16,.9417,.9434,.9449,.9461,.9472,.9481,.9489,.9496,.9502,.9507,.951
     22,.9517,.9520,.9524,.9527,.9530,.9533,.9535,.9537,.9539,.9541,.954
     33,.9545,.9546,.9548,.9549,.9550,.9551,.9552,.9553,.9554,.9555,.955
     46,.9557,.9557,.9558,.9559,.9559,.9560,.9561/                      
      DATA A10/.7773,.8571,.8892,.9068,.9179,.9256,.9312,.9354,.9388,.94
     114,.9436,.9454,.9470,.9482,.9494,.9503,.9512,.9519,.9526,.9531,.95
     237,.9541,.9545,.9549,.9553,.9556,.9559,.9561,.9564,.9566,.9568,.95
     370,.9572,.9573,.9575,.9576,.9577,.9579,.9580,.9581,.9582,.9583,.95
     484,.9585,.9585,.9586,.9587,.9588,.9588,.9589/                     
      DATA A11/.7775,.8575,.8898,.9076,.9188,.9267,.9324,.9367,.9402,.94
     129,.9452,.9470,.9486,.9500,.9511,.9521,.9530,.9538,.9545,.9551,.95
     256,.9561,.9566,.9570,.9573,.9577,.9580,.9583,.9585,.9588,.9590,.95
     392,.9594,.9595,.9597,.9599,.9600,.9601,.9603,.9604,.9605,.9606,.96
     407,.9608,.9609,.9610,.9610,.9611,.9612,.9612/                     
      DATA A14/.7779,.8583,.8910,.9091,.9207,.9288,.9347,.9393,.9429,.94
     159,.9483,.9503,.9520,.9535,.9548,.9559,.9569,.9578,.9585,.9592,.95
     298,.9604,.9609,.9614,.9618,.9622,.9625,.9629,.9632,.9634,.9637,.96
     339,.9642,.9644,.9646,.9647,.9649,.9651,.9652,.9654,.9655,.9656,.96
     457,.9659,.9660,.9661,.9662,.9662,.9663,.9664/                     
      DATA A17/.7781,.8587,.8916,.9099,.9217,.9300,.9361,.9408,.9446,.94
     176,.9502,.9523,.9541,.9557,.9571,.9582,.9593,.9602,.9611,.9618,.96
     225,.9631,.9637,.9642,.9646,.9651,.9655,.9658,.9662,.9665,.9668,.96
     370,.9673,.9675,.9677,.9679,.9681,.9683,.9685,.9686,.9688,.9689,.96
     491,.9692,.9693,.9694,.9695,.9696,.9697,.9698/                     
      DATA A20/.7782,.8590,.8920,.9105,.9224,.9307,.9370,.9418,.9456,.94
     188,.9514,.9536,.9555,.9571,.9586,.9598,.9609,.9619,.9628,.9636,.96
     243,.9649,.9655,.9661,.9666,.9670,.9675,.9679,.9682,.9686,.9689,.96
     392,.9694,.9697,.9699,.9702,.9704,.9706,.9708,.9709,.9711,.9713,.97
     414,.9716,.9717,.9718,.9719,.9721,.9722,.9723/                     
      DATA A25/.7784,.8592,.8924,.9110,.9230,.9315,.9378,.9428,.9467,.95
     100,.9527,.9550,.9569,.9586,.9601,.9614,.9626,.9637,.9646,.9655,.96
     262,.9669,.9676,.9682,.9687,.9692,.9697,.9701,.9705,.9709,.9712,.97
     315,.9718,.9721,.9724,.9727,.9729,.9731,.9733,.9735,.9737,.9739,.97
     441,.9742,.9744,.9745,.9747,.9748,.9749,.9751/                     
      DATA A30/.7784,.8594,.8926,.9113,.9234,.9319,.9383,.9433,.9474,.95
     107,.9534,.9558,.9578,.9595,.9611,.9625,.9637,.9648,.9657,.9666,.96
     274,.9682,.9689,.9695,.9701,.9706,.9711,.9715,.9720,.9724,.9727,.97
     331,.9734,.9737,.9740,.9743,.9745,.9748,.9750,.9752,.9754,.9756,.97
     458,.9760,.9762,.9763,.9765,.9766,.9768,.9769/                     
      DATA A40/.7785,.8595,.8928,.9116,.9237,.9324,.9389,.9440,.9480,.95
     114,.9542,.9567,.9587,.9606,.9622,.9636,.9649,.9660,.9670,.9680,.96
     288,.9696,.9703,.9710,.9716,.9722,.9727,.9732,.9737,.9741,.9745,.97
     349,.9752,.9756,.9759,.9762,.9765,.9768,.9770,.9773,.9775,.9777,.97
     479,.9781,.9783,.9785,.9787,.9788,.9790,.9791/                     
      DATA A50/.7785,.8596,.8929,.9117,.9239,.9326,.9391,.9443,.9484,.95
     118,.9547,.9571,.9592,.9611,.9627,.9642,.9655,.9666,.9677,.9687,.96
     296,.9704,.9711,.9718,.9724,.9730,.9736,.9741,.9746,.9751,.9755,.97
     359,.9763,.9766,.9770,.9773,.9776,.9779,.9781,.9784,.9786,.9789,.97
     491,.9793,.9795,.9797,.9799,.9801,.9803,.9804/                     
      DATA A100/.7785,.8597,.8931,.9119,.9242,.9329,.9395,.9447,.9489,.9
     1523,.9553,.9578,.9599,.9619,.9635,.9651,.9664,.9676,.9687,.9698,.9
     2707,.9716,.9724,.9731,.9738,.9744,.9750,.9756,.9761,.9766,.9771,.9
     3775,.9779,.9783,.9787,.9791,.9794,.9797,.9801,.9803,.9806,.9809,.9
     4812,.9814,.9817,.9819,.9821,.9823,.9825,.9827/                    
      DATA A150/.7786,.8597,.8931,.9119,.9242,.9330,.9396,.9448,.9490,.9
     1525,.9554,.9579,.9601,.9620,.9637,.9652,.9666,.9678,.9690,.9700,.9
     2710,.9718,.9726,.9734,.9741,.9747,.9754,.9759,.9765,.9770,.9775,.9
     3779,.9783,.9787,.9791,.9795,.9799,.9802,.9805,.9808,.9811,.9814,.9
     4817,.9819,.9822,.9824,.9827,.9829,.9831,.9833/                    
      DATA A225/.7786,.8597,.8931,.9120,.9243,.9330,.9396,.9448,.9490,.9
     1525,.9554,.9580,.9602,.9621,.9638,.9653,.9667,.9680,.9691,.9701,.9
     2711,.9720,.9728,.9735,.9742,.9749,.9755,.9761,.9766,.9772,.9776,.9
     3781,.9785,.9790,.9793,.9797,.9801,.9804,.9808,.9811,.9814,.9817,.9
     4819,.9822,.9825,.9827,.9829,.9832,.9834,.9836/                    
      DATA A500/.7786,.8597,.8931,.9120,.9243,.9330,.9396,.9448,.9490,.9
     1525,.9555,.9580,.9602,.9621,.9639,.9654,.9668,.9680,.9692,.9702,.9
     2712,.9720,.9729,.9736,.9743,.9750,.9756,.9762,.9768,.9773,.9778,.9
     3782,.9787,.9791,.9795,.9799,.9802,.9806,.9809,.9812,.9815,.9818,.9
     4821,.9824,.9827,.9829,.9831,.9834,.9836,.9838/                    
      DATA IXV/11,14,17,20,25,30,40,50,100,150,225,507/                 
C                                                                       
      DATA VALUE/.4975936099985107D0,.4873642779856548D0,.46913727600136
     164D0,.4432077635022005D0,.4100009929869515D0,.3700620957892772D0,.
     23240468259684878D0,.2727107356944198D0,.2168967538130226D0,.157521
     33398480817D0,.0955594337368082D0,.0320284464313028D0/             
C                                                                       
C  DATA FOR CALCULATION OF HYDROGENIC RECOMBINATION COEFFICIENTS        
      DATA SV0A/.06845,.07335,.07808,.08268,.08714,.09148,.09570,.09982,
     1.10385,.10778,.11163,.1209,.1297,.1382,.1462,.1540,.1615,.1687,.17
     257,.1824,.1889,.1953,.2015,.2133,.2245,.2352,.2454,.2552,.2646,.27
     336,.2823,.2906,.2987/                                             
      DATA SV0B/.3140,.3284,.3419,.3547,.3668,.3783,.3892,.3996,.4096,.4
     1191,.4413,.4615,.4799,.4968,.5124,.5269,.5404,.5530,.5648,.5759,.5
     2864,.5963,.6146,.6311,.6461,.6598,.6724,.6840,.6947,.7047,.7140,.7
     3226,.7384/                                                        
      DATA SV0C/.7524,.7649,.7761,.7862,.7955,.8039,.8117,.8188,.8254,.8
     1399,.8521,.8626,.8716,.8795,.8865,.8927,.8982,.9032,.9077,.9118,.9
     2156,.9223,.9279,.9328,.9370,.9408,.9441,.9471,.9498,.9522,.9544,2*
     3.9544/                                                            
      DATA SV1A/.00417,.00444,.00469,.00493,.00516,.00538,.00558,.00578,
     1.00597,.00614,.00631,.0067,.0070,.0073,.0076,.0078,.0080,.0082,.00
     283,.0085,.0086,.0087,.0087,.0088,.0088,.0088,.0087,.0086,.0084,.00
     382,.0080,.0077,.0074/                                             
      DATA SV1B/.0068,.0060,.0053,.0044,.0035,.0025,.0016,+.0005,-.0005,
     1-.0016,-.0044,-.0072,-.0101,-.0130,-.0160,-.0190,-.0220,-.0250,-.0
     2279,-.0308,-.0337,-.0366,-.0422,-.0478,-.0532,-.0586,-.0638,-.0689
     3,-.0739,-.0788,-.0835,-.0882,-.0972/                              
      DATA SV1C/-.1058,-.1141,-.1220,-.1295,-.1368,-.1438,-.1506,-.1571,
     1-.1634,-.1784,-.1923,-.2053,-.2174,-.2289,-.2397,-.2499,-.2596,-.2
     2689,-.2777,-.2862,-.2944,-.3099,-.3242,-.3376,-.3502,-.3622,-.3736
     3,-.3844,-.3948,-.4047,-.4147,2*-.4147/                            
      DATA SV2A/-.00120,-.00131,-.00142,-.00153,-.00164,-.00175,-.00186,
     1-.00197,-.00207,-.00218,-.00228,-.0025,-.0028,-.0030,-.0033,-.0035
     2,-.0038,-.0040,-.0043,-.0045,-.0047,-.0050,-.0052,-.0056,-.0061,-.
     30065,-.0070,-.0074,-.0078,-.0082,-.0086,-.0091,.0095/             
      DATA SV2B/-.0103,-.0110,-.0118,-.0126,-.0133,-.0141,-.0148,-.0155,
     1-.0162,-.0170,-.0187,-.0204,-.0221,-.0237,-.0253,-.0269,-.0284,-.0
     2299,-.0314,-.0329,-.0344,-.0359,-.0388,-.0416,-.0444,-.0471,-.0497
     3,-.0523,-.0549,-.0575,-.0600,-.0625,-.0674/                       
      DATA SV2C/-.0722,-.0768,-.0814,-.0859,-.0904,-.0947,-.0989,-.1031,
     1-.1072,-.1174,-.1272,-.1367,-.146,-.155,-.1638,-.1723,-.1807,-.189
     2,-.197,-.2049,-.2127,-.228,-.243,-.257,-.272,-.285,-.299,-.312,-.3
     325,-.337,3*-.35/                                                  
      END                                                               
C                                                                       
      FUNCTION BK (N,NDASH,IS)                                          
C                                                                       
C  CALLS APPROPRIATE ROUTINES FOR CALCULATION OF ATOMIC DATA FOR        
C  ARRAY SK (IN MAIN)                                                   
C                                                                       
C  N=INITIAL LEVEL, NDASH=FINAL LEVEL, IS=SUBSCRIPT WHICH IDENTIFIES    
C  VALUE OF N IN CONDENSED MATRIX                                       
C                                                                       
      COMMON /EXPDAT/ CXP(707),MAXN                                     
      COMMON /PARMS/ DENS,T,ITM                                         
      COMMON /RCRATS/ RADTOT(75),COLTOT(75)                             
      COMMON /TDEP/ TE32,TE12,CTE                                       
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,AN,ANDASH,BK,C,COLTOT,COR,CTE,CX,CXP,DENS,G,RAD
     1TOT,RT,T,TE12,TE32,TEMP                                           
      IF (N-NDASH) 20,10,50                                             
C                                                                       
C  NDASH=N                                                              
 10   CALL COLION (N,1,T,RT)                                            
      BK=-RADTOT(IS)-(COLTOT(IS)+RT)*DENS                               
      IF (N.LE.20) RETURN                                               
      BK=BK+COR(N,3)                                                    
      RETURN                                                            
C                                                                       
C  NDASH GT N                                                           
 20   CALL RAD (A,N,NDASH,G)                                            
      C=COLRAT(N,NDASH,T,TE12)                                          
      AN=N                                                              
      ANDASH=NDASH                                                      
      IF (NDASH.GT.707) GO TO 30                                        
      CX=CXP(NDASH)                                                     
      IF (CX.LT.1.D-30) GO TO 30                                        
      CXN=CXP(N)                                                        
      IF (CXN.LT.1.D-30) GO TO 30                                       
      TEMP=CXN/CX                                                       
      GO TO 40                                                          
 30   TEMP=DEXP(-CTE*(1.D0/AN**2-1.D0/ANDASH**2))                       
 40   TEMP=(ANDASH/AN)**2*TEMP                                          
      BK=(A+DENS*C)*TEMP                                                
      IF (N.LE.20) RETURN                                               
      IF (NDASH.NE.N+1) RETURN                                          
      BK=BK+COR(N,1)                                                    
      RETURN                                                            
C                                                                        
C  NDASH LT N                                                           
 50   C=COLRAT(NDASH,N,T,TE12)                                          
      BK=C*DENS                                                         
      IF (NDASH.NE.N-1) RETURN                                          
      IF (N.LE.20) RETURN                                               
      BK=BK+COR(N,2)                                                    
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION CAPPA (T,TL,T32,F)                                       
C                                                                       
C  COMPUTES FREE-FREE ABSORPTION COEFFICIENT  CAPPA, GIVEN ELECTRON     
C  TEMPERATURE (T), LOG(T) (TL), T**1.5 (T32), AND FREQUENCY IN GHZ (F) 
C                                                                       

C     ADDED BY AAK IN AN ATTEMPT TO COMPILE CODE
C      DOUBLE PRECISION T

      V=0.6529+.6666667*ALOG10(F)-TL                                    
      IF (V.GT.-2.6) GO TO 10                                           
      ALIV=-1.1249*V+0.3788                                             
      GO TO 30                                                          
 10   IF (V.GT.-0.25) GO TO 20                                          
      ALIV=-1.2326*V+0.0987                                             
      GO TO 30                                                          
 20   ALIV=-1.0842*V+0.1359                                             
 30   CAPPA=4.646*EXPM1(.047993*F/T)*EXP(2.302585*ALIV)/(F**2.33333*T32)
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION COLGL (N,NDASH,T,TE12)                                   
C                                                                       
C  CALCULATES COLLISION RATES FROM LEVEL  N  TO HIGHER LEVEL  NDASH AT  
C  ELECTRON TEMPERATURE T. TE12 IS SQRT(T). USES GAUSS-LAGUERRE         
C  INTEGRATION OF CROSS-SECTIONS (FUNCTION CROSS) OVER MAXWELL          
C  DISTRIBUTION. FUNCTION COLGL IS USED FOR VALUES OUTSIDE THE REGION   
C  OF VALIDITY OF FUNCTION COLRAT.                                      
C                                                                       
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DIMENSION XGL(10), WGL(10)                                        
      DOUBLE PRECISION T,TE12                                           
C     DATA NGL/10/                                                      
C     DATA XGL/.1377935,.7294545,1.808343,3.401434,5.552496,8.330153,   
C    A 11.84379,16.27926,21.99659,29.92070/                             
C     DATA WGL/.3084411,.4011199,.2180683,6.208746E-2,9.501517E-3,      
C    A 7.530084E-4,2.825923E-5,4.249314E-7,1.839565E-9,9.911827E-13/    
      DATA XGL/.4157746,2.294280,6.289945,7*0./                         
      DATA WGL/.7110930,.2785177,1.038926E-2,7*0./                      
      DATA NGL/3/                                                       
      BETA=1.58D5/T                                                     
      EN2=N*N                                                           
      END2=NDASH*NDASH                                                  
      DE=1./EN2-1./END2                                                 
      COLGL=0.                                                          
      DO 10 I=1,NGL                                                     
      E=XGL(I)/BETA+DE                                                  
      COLGL=COLGL+WGL(I)*CROSS(N,NDASH,E)*E*BETA                        
 10   CONTINUE                                                          
      COLGL=COLGL*6.21241E5*SNGL(TE12)*(EN2/END2)                       
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLION (N,IONZ,T,QI)                                   
C                                                                       
C  COMPUTES COLLISIONAL IONIZATION RATE, QI, FROM LEVEL  N  FOR IONS    
C  OF EFFECTIVE CHARGE  IONZ  AT ELECTRON TEMPERATURE  T.               
C                                                                       
C  WHEN CALLED WITH  N=0, COLION COMPUTES AND STORES QUANTITIES WHICH   
C  DEPEND ONLY UPON TEMPERATURE AND EFFECTIVE CHARGE. IT IS ASSUMED     
C  THAT  T  AND  IONZ  REMAIN CONSTANT UNTIL THE NEXT CALL WITH  N=0.   
C                                                                       
      COMMON /TDEP/ TE32,TE12,CTE                                       
      DIMENSION EXPX(507)    
C     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION CTE,QI,T,TE12,TE32                               
      QI=0.D0                                                           
      IF (N.NE.0) GO TO 20                                              
C  INITIALIZE                                                           
      CONS=DFLOAT(IONZ*IONZ)*CTE                                        
CC AAK modified to do levels down to 10. lowercase below are additions.
      nlow=9
      nlow1=nlow+1     
C      DO 10 I=21,507  AAK
      do 10 i=nlow1, 507                                                    
 10   EXPX(I)=EXP(-CONS/FLOAT(I*I))
      RETURN                                                            
c 20   IF (N.LE.20) RETURN  AAK                                            
 20   if (n.le.507) return
      X=CONS/DFLOAT(N*N)                                                
      DXP=EXPX(N)                                                       
      IF (N.GT.507) DXP=EXP(-X)                                         
      IF (X.LE.1.) GO TO 30                                             
      E1=((.250621+X*(2.334733+X))/(1.681534+X*(3.330657+X)))*DXP/X     
      GO TO 40                                                          
 30   E1=-.57721566+X*(.9999193+X*(-.24991055+X*(.5519968E-1+X*(-.976004
     1E-2+X*.107857E-2))))-ALOG(X)                                      
 40   EI=DXP*(1.666667-.6666667*X)/X+E1*(.6666667*X-1.)-0.5*E1*E1/DXP   
      QI=5.444089*EI/SNGL(TE32)                                         
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION COLRAT (N,NP,T,TE12)                                     
C                                                                       
C  CALCULATES RATE OF COLLISIONS FROM LEVEL  N  TO HIGHER LEVEL  NP     
C  AT ELECTRON TEMPERATURE  T.  TE12 IS  SQRT(T). SETS RATE=0 FOR  NP-N 
C  GT 40, BUT THIS IS EASILY MODIFIED.                                  
C                                                                       
C  THIS FUNCTION MUST BE INITIALIZED BY BEING CALLED WITH  N=0 BEFORE   
C  ANY COLLISION RATES ARE COMPUTED.                                    
C  THEORY: GEE, PERCIVAL, LODGE AND RICHARDS, MNRAS 175, 209-215 (1976) 
C                                                                       
C  RANGE OF VALIDITY OF GPLR RATES IS  10**6/N**2 LT T LLT 3*10**9      
C  OUTSIDE THIS RANGE, NUMERICAL INTEGRATION OF THE GPLR CROSS-SECTIONS 
C  IS RESORTED TO. THESE CROSS-SECTIONS ARE VALID DOWN TO ENERGIES OF   
C  4/N**2 RYDBERGS; THE CROSS-SECTION FORMULA CAN BE USED AT LOWER      
C  ENERGIES FOR BN CALCULATIONS, THE INACCURACY IN THE CROSS-SECTIONS   
C  HAVING LITTLE EFFECT ON THE BN'S.                                    
C                                                                       
      COMMON /COLINF/ AL18S4(506),S23TRM(506)                           
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION DRT,T,TE12                                       
      REAL L,J1,J2,J3,J4                                                
      IF (N.LE.0) GO TO 60                                              
      COLRAT=0.                                                         
      IS=NP-N                                                           
      IF (IS.GT.40) RETURN                                              
      S=IS                                                              
      EN2=N*N                                                           
      IF (SNGL(T).LT.1.E6/EN2) GO TO 50                                 
      EN=N                                                              
      ENP=NP                                                            
      IPOW=1+IS+IS                                                      
      POW=IPOW                                                          
      ENNP=N*NP                                                         
      BETA=1.58D5/T                                                     
      BETA1=1.4*SQRT(ENNP)                                              
      BETRT=BETA1/BETA                                                  
      BETSUM=BETA1+BETA                                                 
      F1=0.2*S/ENNP                                                     
      IF (F1.GT.0.02) GO TO 10                                          
      F1=1.-POW*F1                                                      
      GO TO 20                                                          
 10   F1=(1.-F1)**IPOW                                                  
 20   A=(2.666667/S)*(ENP/(S*EN))**3*S23TRM(IS)*F1                      
      L=0.85/BETA                                                       
      L=ALOG((1.+0.53*L*L*ENNP)/(1.+0.4*L))                             
      J1=1.333333*A*L*BETRT/BETSUM                                      
      DRT=DSQRT(2.D0-DFLOAT(N*N)/DFLOAT(NP*NP))                         
      F1=0.3*S/ENNP                                                     
      IF (F1.GT.0.02) GO TO 30                                          
      F1=1.-POW*F1                                                      
      GO TO 40                                                          
 30   F1=(1.-F1)**IPOW                                                  
 40   J2=0.                                                             
      ARG=BETA/BETA1                                                    
      IF (ARG.LE.150.) J2=1.777778*F1*(ENP*SNGL(DRT+1.D0)/((EN+ENP)*S))*
     1*3*EXP(-ARG)/(BETA/(1.-AL18S4(IS)))                               
      XI=2./(EN2*SNGL(DRT-1.D0))                                        
      Z=0.75*XI*BETSUM                                                  
      EXPZ=0.                                                           
      IF (Z.LE.150.) EXPZ=EXP(-Z)                                       
      J4=2./(Z*(2.+Z*(1.+EXPZ)))                                        
      J3=0.25*(EN2*XI/ENP)**3*J4*ALOG(1.+0.5*BETA*XI)/BETSUM            
      COLRAT=EN2*EN2*(J1+J2+J3)/SNGL(TE12)**3                           
      RETURN                                                            
C                                                                       
C  GPLR COLLISION RATE FORMULAS INVALID AT LOW TEMPERATURES. INTEGRATE  
C  CROSS-SECTIONS.                                                      
C                                                                       
 50   COLRAT=COLGL(N,NP,T,TE12)                                         
      RETURN                                                            
C                                                                       
C  INITIALIZE                                                           
C                                                                       
 60   DO 70 I=1,506                                                     
      Z=I                                                               
      AL18S4(I)=ALOG(18.*Z)/(4.*Z)                                      
      S23TRM(I)=(0.184-0.04*Z**(-.6666667))                             
 70   CONTINUE                                                          
      COLRAT=0.                                                         
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION COR (N,ISW)                                              
C                                                                       
C  CORRECTION TO RATE OF POPULATION OF LEVEL N DUE TO RADIATION FIELD.  
C  SEE WRITEUP. THIS SUBPROGRAM MAY BE REPLACED BY THE USER IF OTHER    
C  THAN A DILUTE BLACKBODY RADIATION FIELD IS REQUIRED.                 
C                                                                       
C  THE MAIN PROGRAM CALLS COR ONCE WITH ISW=0 BEFORE STARTING THE       
C  CALCULATIONS. ANY INITIALIZATION (INCLUDING THE READING OF DATA IF   
C  REQUIRED) SHOULD BE CARRIED OUT DURING THIS FIRST CALL.              
C                                                                       
C  COMPUTES TERMS OF                                                    
C  D(NU)RHO(NU)(N(N+1)B(N+1,N)+N(N-1)B(N-1,N)-N(N)(B(N,N-1)+B(N,N+1))   
C  WITH REMOVED FACTOR                                                  
C  (C/4PI)(H*H/2PI M KT)**3/2 N*N NE NI EXP(CHI1/N*N KT)                
C                                                                       
      COMMON /EXPDAT/ CXP(707),MAXN                                     
      DIMENSION DILT(508), DX(508)                                      
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,AM,AP,COR,CXP,EX,G                             
      LOGICAL NOFLD                                                     
      DATA DILT/508*0.5/,NOFLD/.FALSE./                                 
C     ADDED BY AAK TO TRY AND GET PROGRAM COMPILED.
C      DOUBLE PRECISION TBCK, EBCK
 
      IF (ISW.GT.0) GO TO 40                                            
      READ (3,*) TBCK,EBCK                                         
      IF (EBCK.EQ.0. .OR. TBCK.EQ.0.) NOFLD=.TRUE.                      
      IF (.NOT.NOFLD) GO TO 10                                          
      GO TO 80                                                          
 10   continue

	write(10,*)'tbck=',sngl(tbck),'  emission measure=',sngl(ebck)

      C15=15.778E4/TBCK                                                 
      MAXP=MAXN+1                                                       
      DO 20 I=1,MAXP                                                    
      AI=I*(I+1)                                                        
      ARG=C15*FLOAT(2*I+1)/AI**2                                        
      DX(I)=EXPM1(ARG)                                                  
 20   CONTINUE                                                          
      IF (EBCK.GE.0.9999E10) GO TO 80                                   
      TBCK32=TBCK*SQRT(TBCK)                                            
      TL=ALOG10(REAL(TBCK))                                                   
      DO 30 I=1,MAXP                                                    
      FREQG=6.58E6/FLOAT(I)**3                                          
      TOW=EBCK*CAPPA(TBCK,TL,TBCK32,FREQG)                              
      IF (TOW.LE.20.) DILT(I)=-0.5*EXPM1(-TOW)                          
 30   CONTINUE                                                          
      GO TO 80                                                          
 40   IF (NOFLD) GO TO 80                                               
      IF (N.GT.MAXN) WRITE (*,120) N                               
      GO TO (50,60,70), ISW                                             
C  N(N+1)B(N+1,N)                                                       
 50   CALL RAD (A,N,N+1,G)                                              
      EX=DX(N)                                                          
      COR=(DFLOAT((N+1)**2)/DFLOAT(N*N))*DILT(N)*A/(EX*(EX+1.D0))       
      RETURN                                                            
C  N(N-1)B(N-1,N)                                                       
 60   CALL RAD (A,N-1,N,G)                                              
      EX=DX(N-1)                                                        
      COR=DILT(N)*A/(EX/(EX+1.D0))                                      
      RETURN                                                            
C  N(N)(B(N,N-1) + B(N,N+1))                                            
 70   CALL RAD (AM,N-1,N,G)                                             
      CALL RAD (AP,N,N+1,EX)                                            
      EX=DX(N-1)                                                        
      G=DX(N)                                                           
      COR=-DILT(N)*(AM/EX+(DFLOAT((N+1)**2)/DFLOAT(N*N))*AP/G)          
      RETURN                                                            
 80   COR=0.D0                                                          
      RETURN                                                            
C                                                                       
 90   FORMAT (2G10.3)                                                   
 100  FORMAT (21H0NO BACKGROUND FIELD.//)                               
 110  FORMAT (32H0RADIATION FIELD - TEMPERATURE =,1PG12.5,21HK, EMISSION
     1 MEASURE =,G12.5//)                                               
 120  FORMAT (32H0*** COR CALLED WITH N TOO LARGE,I6/)                  
      END                                                               
C                                                                       
C                                                                       
      FUNCTION CROSS (N,NP,E)                                           
C                                                                       
C  COMPUTES CROSS-SECTION FOR TRANSITION FROM LEVEL  N  TO HIGHER LEVEL 
C  NP DUE TO COLLISION WITH ELECTRON OF ENERGY  E.                      
C                                                                       
C  THE FORMULA IS VALID FOR ENERGIES IN THE RANGE 4/N**2 LT E LLT 137**2
C  THIS SUBPROGRAM DOES NOT CHECK THAT  E  IS WITHIN THIS RANGE.        
C                                                                       
C  THEORY: GEE, PERCIVAL, LODGE AND RICHARDS, MNRAS 175, 209-215 (1976) 
C                                                                       
      COMMON /COLINF/ AL18S4(506),S23TRM(506)                           
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION DRT                                              
      REAL L                                                            
      C2(X,Y)=X*X*ALOG(1.+.6666667*X)/(Y+Y+1.5*X)                       
      EN=N                                                              
      ENP=NP                                                            
      IS=NP-N                                                           
      IPOW=1+IS+IS                                                      
      POW=IPOW                                                          
      S=IS                                                              
      ENNP=N*NP                                                         
      EENNP=E*E*ENNP                                                    
      EN2=N*N                                                           
      ENN=E*EN2                                                         
      D=0.2*S/ENNP                                                      
      IF (D.GT.0.02) GO TO 10                                           
      D=1.-POW*D                                                        
      GO TO 20                                                          
 10   D=(1.-D)**IPOW                                                    
 20   A=(2.666667/S)*(ENP/(S*EN))**3*S23TRM(IS)*D                       
      D=0.                                                              
      ARG=1./EENNP                                                      
      IF (ARG.LT.150.) D=EXP(-ARG)                                      
      L=ALOG((1.+0.53*EENNP)/(1.+0.4*E))                                
      F=(1.-0.3*S*D/ENNP)**IPOW                                         
      G=0.5*(ENN/ENP)**3                                                
      Y=1./(1.-D*AL18S4(IS))                                            
      DRT=DSQRT(2.D0-DFLOAT(N*N)/DFLOAT(NP*NP))                         
      XP=2./(ENN*SNGL(DRT+1.D0))                                        
      XM=2./(ENN*SNGL(DRT-1.D0))                                        
      H=C2(XM,Y)-C2(XP,Y)                                               
      CROSS=8.797016E-17*(EN2*EN2/E)*(A*D*L+F*G*H)                      
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION DMJ1 (K,KM)                                              
      COMMON /FITDAT/ AFIT(4,4),IVAL(4),NFIT                            
      COMMON /HIGHER/ STORE1(224,4),STORE2(224),VAL(24),B,STORE3(224,5),
     1RTVAL(24),LIMIT                                                   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION AFIT,B,DMJ1,RTVAL,STORE1,STORE2,STORE3,VAL       
      DMJ1=0.D0                                                         
      DO 10 J=1,NFIT                                                    
 10   DMJ1=DMJ1+AFIT(J,KM)*STORE1(K,J)                                  
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION DMJ2 (K)                                                 
      COMMON /FITDAT/ AFIT(4,4),IVAL(4),NFIT                            
      COMMON /HIGHER/ STORE1(224,4),STORE2(224),VAL(24),B,STORE3(224,5),
     1RTVAL(24),LIMIT                                                   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION AFIT,B,DMJ2,RTVAL,STORE1,STORE2,STORE3,SUM,TOT,VA
     1L                                                                 
      DMJ2=0.D0                                                         
      TOT=0.D0                                                          
      DO 20 J=1,NFIT                                                    
      SUM=0.D0                                                          
      DO 10 I=1,NFIT                                                    
 10   SUM=SUM+AFIT(J,I)                                                 
 20   TOT=TOT+SUM*STORE1(K,J)                                           
      DMJ2=1.D0-TOT                                                     
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION DPHI (IQ,LG,ITAU,N)                                      
C                                                                       
C  LAGRANGIAN INTERPOLATION                                             
C                                                                       
      DIMENSION IQ(8)                                                   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,DPHI                                           
      DPHI=0.D0                                                         
      DO 20 M=1,LG                                                      
      IF (M.EQ.ITAU) GO TO 20                                           
      A=1.D0                                                            
      DO 10 L=1,LG                                                      
      IF (L.EQ.ITAU.OR.L.EQ.M) GO TO 10                                 
      A=A*DFLOAT(N-IQ(L))                                               
 10   CONTINUE                                                          
      DPHI=DPHI+A                                                       
 20   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
      FUNCTION EXPM1 (X)                                                
C                                                                       
C  COMPUTES  EXP(X)-1  ACCURATELY, EVEN FOR SMALL  X                    
C                                                                       
      IF (ABS(X).LE.0.1) GO TO 10                                       
      EXPM1=EXP(X)-1.                                                   
      RETURN                                                            
 10   EXPM1=X*(.9999997+X*(0.5+X*(.1667708+X*.4166667E-1)))             
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE HELPME (AID,KM)                                        
      COMMON /HIGHER/ STORE1(224,4),STORE2(224),VAL(24),B,STORE3(224,5),
     1RTVAL(24),LIMIT                                                   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION AID,B,C,POL,RTVAL,STORE1,STORE2,STORE3,SUM,VAL,Y 
      SUM=0.D0                                                          
      DO 10 K=1,LIMIT                                                   
      C=STORE2(K)*STORE3(K,KM)                                          
 10   SUM=SUM+C                                                         
      SUM=SUM-0.5D0*C                                                   
      Y=  .6170614899993600D-2*(POL( 1,KM)+POL( 2,KM))                  
      Y=Y+.1426569431446683D-1*(POL( 3,KM)+POL( 4,KM))                  
      Y=Y+.2213871940870990D-1*(POL( 5,KM)+POL( 6,KM))                  
      Y=Y+.2964929245771839D-1*(POL( 7,KM)+POL( 8,KM))                  
      Y=Y+.3667324070554015D-1*(POL( 9,KM)+POL(10,KM))                  
      Y=Y+.4309508076597664D-1*(POL(11,KM)+POL(12,KM))                  
      Y=Y+.4880932605205694D-1*(POL(13,KM)+POL(14,KM))                  
      Y=Y+.5372213505798282D-1*(POL(15,KM)+POL(16,KM))                  
      Y=Y+.5775283402686280D-1*(POL(17,KM)+POL(18,KM))                  
      Y=Y+.6083523646390170D-1*(POL(19,KM)+POL(20,KM))                  
      Y=Y+.6291872817341415D-1*(POL(21,KM)+POL(22,KM))                  
      Y=Y+.6396909767337608D-1*(POL(23,KM)+POL(24,KM))                  
      AID=SUM+B*Y                                                       
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INTERP (M,CO,VAL,DVAL,IC,IR)                           
C                                                                       
C  COMPUTES SOLUTIONS AT ALL  N  FROM THOSE OBTAINED AT THE CONDENSED   
C  POINTS                                                               
C                                                                       
      DIMENSION CO(75), DVAL(507), IQ(8), M(75), VAL(507)               
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION CO,COT,DFL,DPHI,DVAL,FL,PHI,PHITAU,VAL           
C                                                                       
      LG=2*(IR+1)                                                       
      IB=IC-IR                                                          
      IBB=IB-1                                                          
      ICC=IC-1                                                          
C                                                                       
      K=M(IC)                                                           
      DO 10 I=1,K                                                       
      VAL(I)=0.D0                                                       
      DVAL(I)=0.D0                                                      
 10   CONTINUE                                                          
C                                                                       
      DO 20 IT=1,IC                                                     
      IA=IT                                                             
      MIT=M(IT)                                                         
      VAL(MIT)=CO(IT)                                                   
      IF ((M(IT+1)-MIT).GT.1) GO TO 30                                  
 20   CONTINUE                                                          
C                                                                       
 30   IF (IA.EQ.IC) GO TO 90                                            
      IF (IA.EQ.IB) GO TO 60                                            
C                                                                       
      DO 50 IT=IA,IBB                                                   
      N1=M(IT)+1                                                        
      N2=M(IT+1)                                                        
      ITR1=IT-IR-1                                                      
      DO 40 ITAU=1,LG                                                   
      IND=ITR1+ITAU                                                     
 40   IQ(ITAU)=M(IND)                                                   
      DO 50 ITAU=1,LG                                                   
      PHITAU=1.D0/PHI(IQ,LG,ITAU,IQ(ITAU))                              
      DO 50 N=N1,N2                                                     
      FL=PHI(IQ,LG,ITAU,N)*PHITAU                                       
      DFL=DPHI(IQ,LG,ITAU,N)*PHITAU                                     
      IND=ITR1+ITAU                                                     
      COT=CO(IND)                                                       
      VAL(N)=VAL(N)+FL*COT                                              
 50   DVAL(N)=DVAL(N)+DFL*COT                                           
 60   IF (IR.EQ.0) GO TO 90                                             
C                                                                       
      ICLG=IC-LG                                                        
      DO 70 ITAU=1,LG                                                   
      IND=ICLG+ITAU                                                     
 70   IQ(ITAU)=M(IND)                                                   
      DO 80 ITAU=1,LG                                                   
      PHITAU=1.D0/PHI(IQ,LG,ITAU,IQ(ITAU))                              
      DO 80 IT=IB,ICC                                                   
      N1=M(IT)+1                                                        
      N2=M(IT+1)                                                        
      DO 80 N=N1,N2                                                     
      FL=PHI(IQ,LG,ITAU,N)*PHITAU                                       
      DFL=DPHI(IQ,LG,ITAU,N)*PHITAU                                     
      IND=ICLG+ITAU                                                     
      COT=CO(IND)                                                       
      VAL(N)=VAL(N)+FL*COT                                              
 80   DVAL(N)=DVAL(N)+DFL*COT                                           
C                                                                       
 90   RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE JMD (SK,CO,MVAL,IC)                                    
      COMMON /FITDAT/ AFIT(4,4),IVAL(4),NFIT                            
      COMMON /GAUSS/ VALUE(12)                                          
      COMMON /HIGHER/ STORE1(224,4),STORE2(224),VAL(24),B,STORE3(224,5),
     1RTVAL(24),LIMIT                                                   
      COMMON /PARMS/ DENS,T,ITM                                         
      DIMENSION AZ(4), CO(75), IND(4,2), IPIV(4), MVAL(75), SK(75,75)   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,AC,AFIT,AID,AJ,AK,AKK,AZ,B,BK,CO,D,DENS,DMJ1,DM
     1J2,RTVAL,SK,SOS,STORE1,STORE2,STORE3,T,VAL,VALUE                  
      SOS(I,A)=DSQRT(-A)**(2*I+ITM)/DLOG(-A)                            
      LIMIT=200                                                         
      NG=24                                                             
      DO 10 J=1,NFIT                                                    
      K=IVAL(J)                                                         
      AJ=-1.D0/DFLOAT(MVAL(K))**2                                       
      DO 10 I=1,NFIT                                                    
 10   AFIT(J,I)=SOS(I,AJ)                                               
      CALL MATINV (AFIT,NFIT,AZ,0,D,IRROR,4,IPIV,IND)                   
      B=1.D0/DFLOAT(MVAL(IC)+LIMIT)**2                                  
      A=-0.5D0*B                                                        
      NH=NG/2                                                           
      DO 20 K=1,NH                                                      
      VAL(2*K-1)=A+VALUE(K)*B                                           
 20   VAL(2*K)=A-VALUE(K)*B                                             
      DO 30 K=1,LIMIT                                                   
      A=MVAL(IC)+K                                                      
      AC=-1.D0/A**2                                                     
      DO 30 J=1,NFIT                                                    
 30   STORE1(K,J)=SOS(J,AC)                                             
      DO 40 K=1,NG                                                      
      DO 40 J=1,NFIT                                                    
      INP=K+LIMIT                                                       
 40   STORE1(INP,J)=SOS(J,VAL(K))                                       
      KK=LIMIT+NG                                                       
      DO 60 K=1,KK                                                      
      DO 50 J=1,NFIT                                                    
 50   STORE3(K,J)=DMJ1(K,J)                                             
 60   STORE3(K,NFIT+1)=DMJ2(K)                                          
      DO 100 J=1,IC                                                     
      I=MVAL(J)                                                         
      DO 70 K=1,LIMIT                                                   
      KK=MVAL(IC)+K                                                     
      AKK=KK                                                            
 70   STORE2(K)=BK(I,KK,0)                                              
      DO 80 K=1,NG                                                      
      AK=DSQRT(-1.D0/VAL(K))                                            
      RTVAL(K)=AK**3*0.5D0                                              
      KK=AK                                                             
      INP=K+LIMIT                                                       
      STORE2(INP)=BK(I,KK,0)                                            
 80   CONTINUE                                                          
      CALL HELPME (AID,NFIT+1)                                          
      CO(J)=CO(J)-AID                                                   
      DO 90 KM=1,NFIT                                                   
      CALL HELPME (AID,KM)                                              
      L=IVAL(KM)                                                        
 90   SK(J,L)=SK(J,L)+AID                                               
 100  CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE MATINV (A,N,B,L,D,IRROR,NDA,IPIV,IND)                  
      DIMENSION A(NDA,NDA), B(NDA), IND(NDA,2), IPIV(NDA)               
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,AMAX,ATEMP,B,D                                 
C                                                                       
C     SOLVES SIMULTANEOUS EQUATIONS IF L=1                              
C                                                                       
C     INVERTS MATRIX A IF L=0                                           
C                                                                       
C     SOLUTIONS ARE RETURNED IN B                                       
C                                                                       
      M=IABS(L)                                                         
      D=1.D0                                                            
      DO 10 I=1,N                                                       
 10   IPIV(I)=0                                                         
      DO 190 I=1,N                                                      
      AMAX=0.D0                                                         
      DO 60 J=1,N                                                       
      IF (IPIV(J)) 70,20,60                                             
 20   DO 50 K=1,N                                                       
      IF (DABS(A(J,K)).LT.1.D-50) A(J,K)=0.D0                           
      IF (IPIV(K)-1) 30,50,70                                           
 30   IF (DABS(A(J,K))-AMAX) 50,50,40                                   
 40   IROW=J                                                            
      ICOL=K                                                            
      AMAX=DABS(A(J,K))                                                 
 50   CONTINUE                                                          
 60   CONTINUE                                                          
      IPIV(ICOL)=IPIV(ICOL)+1                                           
      IF (AMAX-1.D-50) 70,70,80                                         
 70   IRROR=1                                                           
      RETURN                                                            
 80   IF (IROW-ICOL) 90,120,90                                          
 90   D=-D                                                              
      DO 100 K=1,N                                                      
      AMAX=A(IROW,K)                                                    
      A(IROW,K)=A(ICOL,K)                                               
 100  A(ICOL,K)=AMAX                                                    
      IF (M) 120,120,110                                                
 110  AMAX=B(IROW)                                                      
      B(IROW)=B(ICOL)                                                   
      B(ICOL)=AMAX                                                      
 120  IND(I,1)=IROW                                                     
      IND(I,2)=ICOL                                                     
      AMAX=A(ICOL,ICOL)                                                 
      A(ICOL,ICOL)=1.D0                                                 
      DO 130 K=1,N                                                      
 130  A(ICOL,K)=A(ICOL,K)/AMAX                                          
      IF (M) 150,150,140                                                
 140  B(ICOL)=B(ICOL)/AMAX                                              
 150  DO 190 J=1,N                                                      
      IF (J-ICOL) 160,190,160                                           
 160  AMAX=A(J,ICOL)                                                    
      A(J,ICOL)=0.D0                                                    
      DO 170 K=1,N                                                      
      ATEMP=A(ICOL,K)*AMAX                                              
 170  A(J,K)=A(J,K)-ATEMP                                               
      IF (M) 190,190,180                                                
 180  B(J)=B(J)-B(ICOL)*AMAX                                            
 190  CONTINUE                                                          
      IF (L) 200,200,240                                                
 200  DO 230 I=1,N                                                      
      J=N+1-I                                                           
      IF (IND(J,1)-IND(J,2)) 210,230,210                                
 210  IROW=IND(J,1)                                                     
      ICOL=IND(J,2)                                                     
      DO 220 K=1,N                                                      
      AMAX=A(K,IROW)                                                    
      A(K,IROW)=A(K,ICOL)                                               
 220  A(K,ICOL)=AMAX                                                    
 230  CONTINUE                                                          
 240  IRROR=0                                                           
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION PHI (IQ,LG,ITAU,N)                                       
C                                                                       
C  LAGRANGIAN INTERPOLATION                                             
C                                                                       
      DIMENSION IQ(8)                                                   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION PHI                                              
      PHI=1.D0                                                          
      DO 10 L=1,LG                                                      
      IF (L.EQ.ITAU) GO TO 10                                           
      PHI=PHI*DFLOAT(N-IQ(L))                                           
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION POL (K,KM)                                               
      COMMON /HIGHER/ STORE1(224,4),STORE2(224),VAL(24),B,STORE3(224,5),
     1RTVAL(24),LIMIT                                                   
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION B,POL,RTVAL,STORE1,STORE2,STORE3,VAL             
      IND=K+LIMIT                                                       
      POL=RTVAL(K)*STORE3(IND,KM)*STORE2(IND)                           
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE RAD (A,N,NDASH,G)                                      
C                                                                       
C  COMPUTES EINSTEIN RADIATIVE RATES A AND GAUNT FACTORS G FOR          
C  TRANSITIONS FROM  NDASH  TO  N                                       
C                                                                       
      COMMON /GAUNTS/ GAUNT(50,22),IXV(12)                              
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,ALF2,EN,EN2,END2,G,TERM1,TERM2,TERM3           
C  ARRAY GAUNT IS NOT DOUBLE PRECISION                                  
      EN=N                                                              
      EN2=N*N                                                           
      END2=NDASH*NDASH                                                  
      NDMN=NDASH-N                                                      
      IF (NDMN.GT.50) GO TO 30                                          
      IF (N.LE.10) GO TO 20                                             
      DO 10 J=2,12                                                      
      IVJ=IXV(J)                                                        
      IVJ1=IXV(J-1)                                                     
      IF (N.GT.IVJ) GO TO 10                                            
      P=1./FLOAT(IVJ-IVJ1)                                              
      Q=FLOAT(IVJ-N)*P                                                  
      P=FLOAT(N-IVJ1)*P                                                 
      G=Q*GAUNT(NDMN,J+9)+GAUNT(NDMN,J+10)*P                            
      GO TO 40                                                          
 10   CONTINUE                                                          
 20   G=GAUNT(NDMN,N)                                                   
      GO TO 40                                                          
 30   ALF2=(EN2/END2)                                                   
      TERM1=(1.D0-ALF2)                                                 
      TERM1=(TERM1*EN)**.666666666666666D0                              
      TERM2=0.1728D0*(1.D0+ALF2)/TERM1                                  
      TERM3=0.0496D0*(1.D0-1.3333333333333D0*ALF2+ALF2*ALF2)/TERM1**2   
      G=1.D0-TERM2-TERM3                                                
 40   A=G*15.7457D9/(EN*(END2-EN2)*DFLOAT(NDASH)*END2)                  
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE RADCOL (T,MVAL,IC,NMIN)                                
      COMMON /EXPDAT/ CXP(707),MAXN                                     
      COMMON /RCRATS/ RADTOT(75),COLTOT(75)                             
      COMMON /TDEP/ TE32,TE12,CTE                                       
      DIMENSION MVAL(75)                                                
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION A,AKK,AN,C,COLTOT,CTE,CX,CXP,G,Q,RADTOT,T,TE12,TE
     132,TOT                                                            
C                                                                       
C     RADIATIVE CASCADE COEFFICIENTS                                    
C                                                                       
C     SUM FROM N TO ALL LEVELS DOWN TO NMIN (=1 OR 2)                   
C                                                                       
      DO 20 M=1,IC                                                      
      J=MVAL(M)                                                         
      TOT=0.D0                                                          
      IF (J.LE.NMIN) GO TO 20                                           
      K=J-1                                                             
      DO 10 I=NMIN,K                                                    
      CALL RAD (A,I,J,G)                                                
 10   TOT=TOT+A                                                         
 20   RADTOT(M)=TOT                                                     
C                                                                       
C     COLLISIONAL RATE TOTALS ON TO LEVEL N                             
C                                                                       
C     SUMMED FROM NMIN TO INFINITY                                      
C                                                                       
C       AAK: Modifications to go down to a quantum number of 20 are
C       below in lowercase.
c     AAK
      nlo=nmin
      DO 70 M=1,IC                                                      
      N=MVAL(M)                                                         
      AN=N                                                              
      TOT=0.D0                                                          
      L=N-1                                                             
c      DO 30 KK=20,L   AAK                                                     
      do 30 kk=nlo,l
      C=COLRAT(KK,N,T,TE12)                                             
 30   TOT=TOT+C                                                         
      L=N+1                                                             
      MAX=N+40                                                          
      DO 60 KK=L,MAX                                                    
      AKK=KK                                                            
      C=COLRAT(N,KK,T,TE12)                                             
      CX=CXP(KK)                                                        
      IF (CX.LT.1.D-30) GO TO 40                                        
      CXN=CXP(N)                                                        
      IF (CXN.LT.1.D-30) GO TO 40                                       
      Q=CXN/CX                                                          
      GO TO 50                                                          
 40   Q=DEXP(15.778D4*(1.D0/AKK**2-1.D0/AN**2)/T)                       
 50   C=C*(AKK/AN)**2*Q                                                 
 60   TOT=TOT+C                                                         
 70   COLTOT(M)=TOT                                                     
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE RECOMB (Z,ETEMP,N1,ALPHA)                              
C                                                                       
C  COMPUTES RECOMBINATION COEFFICIENT, ALPHA, ONTO LEVEL  N1  FOR IONS  
C  OF EFFECTIVE CHARGE  Z  AT ELECTRON TEMPERATURE  ETEMP               
C                                                                       
      COMMON /RCMB/ SV0(99),SV1(99),SV2(99)                             
C  SV0, SV1, SV2 NOT DOUBLE PRECISION                                   
      DIMENSION XV(99)                                                  
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION ALPHA,CONST,ETEMP,F,FL,P,Q,S0,S1,S2,TE,U,V,X,XV,X
     1VI,XVI1,XXX,Y,Z2,Z                                                
      Z2=Z*Z                                                            
      TE=ETEMP*1.0D-4                                                   
      CONST=15.778D0/TE                                                 
      FL=1.D0/(CONST*Z2)**.33333333333333D0                             
      CONST=5.197D-14*Z2*DSQRT(CONST)                                   
      XV(1)=.02D0                                                       
      DO 10 N=2,11                                                      
 10   XV(N)=XV(N-1)+.002D0                                              
      DO 20 N=12,23                                                     
 20   XV(N)=XV(N-1)+.005D0                                              
      DO 30 N=24,33                                                     
 30   XV(N)=XV(N-1)+.01D0                                               
      DO 40 N=34,65                                                     
 40   XV(N)=XV(N-32)*10.D0                                              
      DO 50 N=66,97                                                     
 50   XV(N)=XV(N-32)*10.D0                                              
      F=N1                                                              
      X=15.778D0*Z2/(F*F*TE)                                            
      IF (X.LT.0.02D0) GO TO 80                                         
      IF (X.GT.20.D0) GO TO 70                                          
      DO 60 I=2,99                                                      
      XVI=XV(I)                                                         
      IF (X.GT.XVI) GO TO 60                                            
      IM1=I-1                                                           
      XVI1=XV(IM1)                                                      
      P=1.D0/(XVI-XVI1)                                                 
      Q=(XVI-X)*P                                                       
      P=(X-XVI1)*P                                                      
      S0=P*SV0(I)+Q*SV0(IM1)                                            
      S1=P*SV1(I)+Q*SV1(IM1)                                            
      S2=P*SV2(I)+Q*SV2(IM1)                                            
      GO TO 90                                                          
 60   CONTINUE                                                          
 70   U=1.D0/X                                                          
      V=U/3.D0                                                          
      XXX=X**.333333333333333D0                                         
      S0=1.D0-U*(1.D0-2.D0*U*(1.D0-3.D0*U*(1.D0-4.D0*U)))               
      S1=-.1728D0*XXX*(1.D0-V*(8.D0-V*(70.D0-V*(800.D0-V*11440.D0))))   
      S2=-.0496D0*XXX**2*(1.D0-V*(3.D0-V*(32.D0-V*448.D0)))             
      GO TO 90                                                          
 80   S0=X*DEXP(X)*(-DLOG(X)-.5772D0+X)                                 
      XXX=X**.333333333333333D0                                         
      S1=.4629D0*X*(1.D0+4.D0*X)-1.0368D0*XXX**4*(1.D0+1.875D0*X)       
      S2=-.0672D0*X*(1.D0+3.D0*X)+.1488D0*XXX**5*(1.D0+1.8D0*X)         
 90   Y=(S0+FL*(S1+FL*S2))/F                                            
      ALPHA=CONST*Y                                                     
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE REDUCE (M,IC,IR,SK)                                    
C                                                                       
C   GIVEN A SET OF INTEGERS                                             
C         M(IT),IT=1,IC,SUCH THAT-                                      
C         1) M(IT+1)=M(IT)+1  FOR IT.LE.IA                              
C         WHERE IA.GE.1  AND                                            
C         2) (M(IT+1) - M(IT)).GT.1 FOR IT.GE.IA,                       
C   AND GIVEN A FUNCTION SUBPROGRAM  BK                                 
C         WHICH CALCULATES THE                                          
C         ELEMENTS OF A LARGE                                           
C         M(IC)*M(IC) MATRIX,                                           
C   THIS SUBROUTINE USES LAGRANGE                                       
C         INTERPOLATION OF ORDER                                        
C         2*(IR+1) TO CALCULATE A                                       
C         SMALLER IC*IC MATRIX SK                                       
C   REQUIRES A FUNCTION SUBPROGRAM                                      
C         PHI                                                           
C   IR MUST BE .LE. (IA-1)                                              
C                                                                       
      DIMENSION IQ(8), M(75), SK(75,75)                                 
      DIMENSION STORE1(8), STORE2(8)                                    
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION BK,DUCKIT,FL,PHI,PHITAU,SK,STORE1,STORE2         
C                                                                       
      LG=2*(IR+1)                                                       
      IB=IC-IR                                                          
      IBB=IB-1                                                          
      ICC=IC-1                                                          
C                                                                       
      DO 10 IS=1,IC                                                     
      DO 10 IT=1,IC                                                     
 10   SK(IS,IT)=BK(M(IS),M(IT),IS)                                      
C                                                                       
      DO 20 IT=1,IC                                                     
      IA=IT                                                             
      IF ((M(IT+1)-M(IT)).GT.1) GO TO 30                                
 20   CONTINUE                                                          
C                                                                       
 30   IF (IA.EQ.IC) GO TO 110                                           
      IF (IA.EQ.IB) GO TO 80                                            
C                                                                       
      DO 70 IT=IA,IBB                                                   
      N1=M(IT)+1                                                        
      N2=M(IT+1)-1                                                      
      DO 40 ITAU=1,LG                                                   
      IND=IT-IR-1+ITAU                                                  
 40   IQ(ITAU)=M(IND)                                                   
      DO 50 ITAU=1,LG                                                   
 50   STORE1(ITAU)=PHI(IQ,LG,ITAU,IQ(ITAU))                             
      DO 70 N=N1,N2                                                     
      DO 60 ITAU=1,LG                                                   
 60   STORE2(ITAU)=PHI(IQ,LG,ITAU,N)                                    
      DO 70 IS=1,IC                                                     
      DUCKIT=BK(M(IS),N,IS)                                             
      DO 70 ITAU=1,LG                                                   
      FL=STORE2(ITAU)/STORE1(ITAU)                                      
      IND=IT-IR-1+ITAU                                                  
 70   SK(IS,IND)=SK(IS,IND)+DUCKIT*FL                                   
C                                                                       
 80   IF (IR.EQ.0) GO TO 110                                            
C                                                                       
      DO 90 ITAU=1,LG                                                   
      IND=IC-LG+ITAU                                                    
 90   IQ(ITAU)=M(IND)                                                   
      DO 100 ITAU=1,LG                                                  
      PHITAU=1.D0/PHI(IQ,LG,ITAU,IQ(ITAU))                              
      DO 100 IT=IB,ICC                                                  
      N1=M(IT)+1                                                        
      N2=M(IT+1)-1                                                      
      DO 100 N=N1,N2                                                    
      FL=PHI(IQ,LG,ITAU,N)*PHITAU                                       
      DO 100 IS=1,IC                                                    
      IND=IC-LG+ITAU                                                    
 100  SK(IS,IND)=SK(IS,IND)+BK(M(IS),N,IS)*FL                           
C                                                                       
 110  RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE RHS (CO,MVAL,IC)                                       
C                                                                       
C  COMPUTES THE RIGHT HAND SIDE OF EQUATIONS (2.7) OF BROCKLEHURST,     
C  MNRAS 148, 417 (1970).                                               
C                                                                       
      COMMON /EXPDAT/ CXP(707),MAXN                                     
      COMMON /PARMS/ DENS,T,ITM                                         
      COMMON /TDEP/ TE32,TE12,CTE                                       
      DIMENSION MVAL(75), CO(75)                                        
c     DOUBLE PRECISION DABS,DSQRT,DBLE,DLOG,DLOG10,DFLOAT,DEXP          
      DOUBLE PRECISION ALFA,CO,CTE,CXP,DENS,RT,T,TE12,TE32              
      DO 10 I=1,IC                                                      
      J=MVAL(I)                                                         
      CALL COLION (J,1,T,RT)                                            
      CALL RECOMB (1.D0,T,J,ALFA)                                       
      CO(I)=-ALFA*CXP(J)*TE32*0.24146879D16/DFLOAT(J*J)-RT*DENS         
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C /*                                                                      
C //GO.FT07F001 DD UNIT=DISC,DSN=&&CARDS,DISP=(NEW,PASS),                 
C //  DCB=FB,SPACE=(CYL,(1,1))  CARD PUNCH OUTPUT SENT TO TEMPORARY FILE  
C //GO.SYSIN DD *  DATA INPUT CARDS                                       
C BN PROGRAM TEST                                                         
C 75  2  4                                                              
C 30 31 32 33 34 35 37 39 41 43 46 49 52 55 58 61 64 68 72 76 80 84 88  
C 92 97102107112117122127132138144150156162168174180187194201208215222  
C230238246254262270279288297306315325335345355365375386397408419430441  
C452463474485496507                                                     
C 75 72 69 66                                                           
C 10000.     10000000.                                                    
C 500.        1.          2                                               
C 500.        0.01        2    2   50  249                                
C                                                                       
C /*                                                                      
C //CARDS EXEC FILE  PRINT CARD O/P AFTER PRINTER O/P. NO CARDS PUNCHED   
C //FROM DD DSN=&&CARDS,DISP=(OLD,DELETE)                                 
C //TO DD SYSOUT=A                                                        
                                                                        
