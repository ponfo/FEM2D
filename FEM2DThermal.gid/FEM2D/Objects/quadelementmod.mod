  ã]  Ý   k820309    9          19.0        (áÑ]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DThermal.gid/FEM2D/Source/QuadElement.f90 QUADELEMENTMOD              QUADELEMENTTYPE                                                     
                            @                              
                            @                              
                            @                              
                            @                              
                         @                                '                     #ID    #X    #Y 	   #Z 
   #INIT    #SETID    #GETID    #SETX    #GETX    #SETY     #GETY $   #SETZ '   #GETZ +                 D                                                              D                                            
                 D                             	               
                 D                             
               
   1         À    $                                              #INIT    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                                     #POINTTYPE              
                                                      
                                      
                
                                      
                
                                      
      1         À    $                                              #SETID    #         @     @                                                #THIS    #ID              
                                                     #POINTTYPE              
                                            1         À    $                                             #GETID    %         @   @                                                      #THIS              
                                                     #POINTTYPE    1         À    $                                              #SETX    #         @     @                                                #THIS    #X              
                                                     #POINTTYPE              
                                      
      1         À    $                                        	     #GETX    %         @   @                                               
       #THIS              
                                                     #POINTTYPE    1         À    $                                          
     #SETY !   #         @     @                            !                    #THIS "   #Y #             
                                "                     #POINTTYPE              
                                 #     
      1         À    $                           $                  #GETY %   %         @   @                           %                    
       #THIS &             
                                &                     #POINTTYPE    1         À    $                            '                  #SETZ (   #         @     @                            (                    #THIS )   #Z *             
                                )                     #POINTTYPE              
                                 *     
      1         À    $                           +              	    #GETZ ,   %         @   @                           ,                    
       #THIS -             
                                -                     #POINTTYPE                   @  @                           .     '              
      #PTR /   #ALLOCATE 0   #SETID 4   #GETID 8   #SETX ;   #GETX ?   #SETY B   #GETY F   #SETZ I   #GETZ M                $                              /                            #POINTTYPE    1         À    $                            0                  #ALLOCATE 1   #         @     @                            1                    #THIS 2   #POINT 3             
                                2                    #POINTPTRTYPE .             
                                  3                    #POINTTYPE    1         À    $                            4                  #SETID 5   #         @     @                            5                    #THIS 6   #ID 7             
                                6                    #POINTPTRTYPE .             
                                 7           1         À    $                           8                  #GETID 9   %         @   @                           9                           #THIS :             
                                :                    #POINTPTRTYPE .   1         À    $                            ;                  #SETX <   #         @     @                            <                    #THIS =   #X >             
                                =                    #POINTPTRTYPE .             
                                 >     
      1         À    $                          ?                  #GETX @   %         @   @                          @                    
       #THIS A             
                                A                    #POINTPTRTYPE .   1         À    $                            B                  #SETY C   #         @     @                            C                    #THIS D   #Y E             
                                D                    #POINTPTRTYPE .             
                                 E     
      1         À    $                          F                  #GETY G   %         @   @                          G                    
       #THIS H             
                                H                    #POINTPTRTYPE .   1         À    $                            I             	     #SETZ J   #         @     @                            J                    #THIS K   #Z L             
                                K                    #POINTPTRTYPE .             
                                 L     
      1         À    $                           M             
 	    #GETZ N   %         @   @                           N                    
       #THIS O             
                                O                    #POINTPTRTYPE .                 @  @               D          P     'Ð                   #NPOINT Q   #NDOF R   #POINT S   #INTEGRATOR T   #MATERIAL t   #GEOMETRY    #GETNPOINT    #GETNDOF    #GETPOINTID    #GETINTEGRATOR    #GETMATERIAL    #GETGEOMETRY    #SETNPOINT    #SETNDOF    #SETPOINT £   #SETINTEGRATOR ¨   #GETONEPOINT ¬   #GETALLPOINTS °                 $                             Q                                 $                             R                              $                              S                                #POINTPTRTYPE .             &                                                         $                              T            P              #INTEGRATORPTRTYPE U                 @  @                         U     '                    #PTR V   #ALLOCATE p                $                             V                          #INTEGRATORTYPE W                  @  @               D           W     '                   #GAUSSORDER X   #INTEGTERMS Y   #WEIGHT Z   #GPOINT [   #SHAPEFUNC \   #DSHAPEFUNC ]   #INIT ^   #VALUEGPOINTS c   #GETG1D g   #GETGTRIANGLE j   #GETGSQUARE m                 $                             X                                 $                             Y                             $                             Z                             
            &                                                      $                             [            P                 
            &                   &                                                       $                             \            °                 
            &                   &                                                       $                             ]                            
            &                   &                   &                                           1         À    $                            ^                  #INIT _   #         @     @                            _                    #THIS `   #GAUSSORDER a   #TYPE b             
                                `                   #INTEGRATORTYPE W             
                                 a                     
                                b                    1 1         À    $                            c                  #VALUEGPOINTS d   #         @     @                            d                    #THIS e   #TYPE f             
                                e                   #INTEGRATORTYPE W             
                                f                    1 1         À    D                           g             	     #GETG1D h   #         @     @                            h                    #THIS i             
                                i                   #INTEGRATORTYPE W   1         À    D                           j             
     #GETGTRIANGLE k   #         @     @                            k                    #THIS l             
                                l                   #INTEGRATORTYPE W   1         À    D                           m                  #GETGSQUARE n   #         @     @                            n                    #THIS o             
                                o                   #INTEGRATORTYPE W   1         À    $                            p                  #ALLOCATE q   #         @     @                            q                    #THIS r   #INTEGRATOR s             
                                r                    #INTEGRATORPTRTYPE U             
                                  s                  #INTEGRATORTYPE W                 $                              t            Ð              #MATERIALPTRTYPE u                 @  @                         u     '                    #PTR v   #ALLOCATE                 $                             v                           #THMATERIALTYPE w                 @  @                         w     '                    #MATERIALTYPE x   #CONDUCTIVITY z   #INIT {                 $                              x                            #MATERIALTYPE y                 @  @                          y     '                                      $                             z                              
  p          p            p                          1         À    $                            {                  #INIT |   #         @     @                            |                    #THIS }   #KX ~   #KY              
                                }                    #THMATERIALTYPE w             
                                 ~     
                
                                      
      1         À    $                                              #ALLOCATE    #         @     @                                                #THIS    #MATERIAL                                                                  #MATERIALPTRTYPE u             
                                                     #THMATERIALTYPE w                 $                                          P             #GEOMETRYPTRTYPE                  @  @                              '                    #PTR                 $                                                         #GEOMETRYTYPE                   @  @                               '                        1         À    $                                             #GETNPOINT    %         @   @                                                      #THIS              
                                     Ð              #ELEMENTTYPE P   1         À    $                                             #GETNDOF    %         @   @                                                      #THIS              
                                     Ð              #ELEMENTTYPE P   1         À    $                                        	     #GETPOINTID    %         @   @                                                      #THIS    #I              
                                     Ð              #ELEMENTTYPE P             
                                            1         À    $                                        
     #GETINTEGRATOR    &         @   @                                                      #THIS    #INTEGRATORPTRTYPE U             
                                     Ð              #ELEMENTTYPE P   1         À    $                                             #GETMATERIAL    &         @   @                                                      #THIS    #MATERIALPTRTYPE u             
                                     Ð              #ELEMENTTYPE P   1         À    $                                             #GETGEOMETRY    &         @   @                                                      #THIS    #GEOMETRYPTRTYPE              
                                     Ð              #ELEMENTTYPE P   1         À    $                                              #SETNPOINT    #         @     @                                                #THIS    #NPOINT              
                                     Ð              #ELEMENTTYPE P             
                                            1         À    $                                              #SETNDOF     #         @     @                                                 #THIS ¡   #NDOF ¢             
                                ¡     Ð              #ELEMENTTYPE P             
                                 ¢           1         À    $                            £              	    #SETPOINT ¤   #         @     @                            ¤                    #THIS ¥   #I ¦   #POINT §             
                                ¥     Ð              #ELEMENTTYPE P             
                                 ¦                     
                                 §                    #POINTTYPE    1         À    $                            ¨              
    #SETINTEGRATOR ©   #         @     @                            ©                    #THIS ª   #INTEGRATOR «             
                                ª     Ð              #ELEMENTTYPE P             
                                  «                  #INTEGRATORTYPE W   1         À    $                           ¬                  #GETONEPOINT ­   &         @   @                           ­                           #THIS ®   #I ¯   #POINTPTRTYPE .             
                                ®     Ð              #ELEMENTTYPE P             
                                 ¯           1         À    $                           °                  #GETALLPOINTS ±   )        `   @                            ±                                         #THIS ²   #POINTPTRTYPE .   p          5 8 O#ELEMENTTYPE P    p        U  P   Q       5 8 O#ELEMENTTYPE P    p        U  P   Q                               
                                ²     Ð              #ELEMENTTYPE P                     @                        ³     'Ø                   #THELEMENT2DTYPE ´   #SETAREA Ù                 $                              ´     Ø                     #THELEMENT2DTYPE µ                 @  @                        µ     'Ø                   #ELEMENT2DTYPE ¶   #GETSTIFFNESS Ö                 $                              ¶     Ø                     #ELEMENT2DTYPE ·                  @  @                        ·     'Ø             	      #ELEMENTTYPE ¸   #AREA ¹   #GETAREA º   #GETSTIFFNESS ½   #SETAREA À   #SHAPEFUNC Ã   #DSHAPEFUNC È   #JACOBIAN Í   #JACOBIANDET Ò                 $                              ¸     Ð                     #ELEMENTTYPE P                 $                             ¹     Ð         
   1         À    $                           º                  #GETAREA »   %         @   @                           »                    
       #THIS ¼             
                                ¼     Ø              #ELEMENT2DTYPE ·   1         À    $                          ½                  #GETSTIFFNESSINTERF ¾   (        `   @                          ¾                                    
    #THIS ¿     p         5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R   p           5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R      5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R        5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R      5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R                               
                               ¿     Ø              #ELEMENT2DTYPE ·   1         À    $                           À                  #SETAREAINTERF Á   #         @     @                           Á     	               #THIS Â             
                               Â     Ø              #ELEMENT2DTYPE ·   1         À    $                          Ã                  #SHAPEFUNCINTERF Ä   (        `   @                          Ä                                    
    #THIS Å   #X Æ   #Y Ç   p           5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R        5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R                               
                               Å     Ø              #ELEMENT2DTYPE ·             
                                Æ     
                
                                Ç     
      1         À    $                          È                  #DSHAPEFUNCINTERF É   (        `   @                          É                                    
    #THIS Ê   #X Ë   #Y Ì   p          p           5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R       p           5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 O#ELEMENT2DTYPE ·    p        U #ELEMENTTYPE P    ·   ¸   U  P   R                               
                               Ê     Ø              #ELEMENT2DTYPE ·             
                                Ë     
                
                                Ì     
      1         À    $                           Í                  #JACOBIAN Î   (         `   @                           Î                                   
    #THIS Ï   #X Ð   #Y Ñ   p          p          p            p          p                                    
                                Ï     Ø              #ELEMENT2DTYPE ·             
                                 Ð     
                
                                 Ñ     
      1         À    $                           Ò             	     #JACOBIANDET Ó   %         @   @                           Ó                    
       #THIS Ô   #JACOBIAN Õ             
                                Ô     Ø              #ELEMENT2DTYPE ·             
                                 Õ                   
    p          p          p            p          p                          1         À    $                          Ö                  #GETSTIFFNESS ×   (        `   @                           ×                                    
    #THIS Ø     p         5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   R   p           5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   R      5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   R        5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   R      5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   Q   5 8 8 8 O#THELEMENT2DTYPE µ    p        U #ELEMENT2DTYPE ·    µ   ¶   U #ELEMENTTYPE P    ·   ¸   U  P   R                               
                                Ø     Ø              #THELEMENT2DTYPE µ   1         À    $                           Ù                  #SETAREA Ú   #         @     @                             Ú                    #THIS Û             
D @                              Û     Ø              #QUADELEMENTTYPE ³          u      fn#fn $         b   uapp(QUADELEMENTMOD    5  @   J  TOOLS    u  @   J  POINTMOD    µ  @   J  POINTPTRMOD    õ  @   J  ELEMENT2DMOD    5  @   J  THELEMENT2DMOD #   u  É       POINTTYPE+POINTMOD )   >  H   %   POINTTYPE%ID+POINTMOD=ID '     H   %   POINTTYPE%X+POINTMOD=X '   Î  H   %   POINTTYPE%Y+POINTMOD=Y '     H   %   POINTTYPE%Z+POINTMOD=Z (   ^  R   a   POINTTYPE%INIT+POINTMOD    °  o      INIT+POINTMOD #     W   a   INIT%THIS+POINTMOD !   v  @   a   INIT%ID+POINTMOD     ¶  @   a   INIT%X+POINTMOD     ö  @   a   INIT%Y+POINTMOD     6  @   a   INIT%Z+POINTMOD )   v  S   a   POINTTYPE%SETID+POINTMOD    É  Z      SETID+POINTMOD $   #  W   a   SETID%THIS+POINTMOD "   z  @   a   SETID%ID+POINTMOD )   º  S   a   POINTTYPE%GETID+POINTMOD      Z      GETID+POINTMOD $   g  W   a   GETID%THIS+POINTMOD (   ¾  R   a   POINTTYPE%SETX+POINTMOD    	  Y      SETX+POINTMOD #   i	  W   a   SETX%THIS+POINTMOD     À	  @   a   SETX%X+POINTMOD (    
  R   a   POINTTYPE%GETX+POINTMOD    R
  Z      GETX+POINTMOD #   ¬
  W   a   GETX%THIS+POINTMOD (     R   a   POINTTYPE%SETY+POINTMOD    U  Y      SETY+POINTMOD #   ®  W   a   SETY%THIS+POINTMOD       @   a   SETY%Y+POINTMOD (   E  R   a   POINTTYPE%GETY+POINTMOD      Z      GETY+POINTMOD #   ñ  W   a   GETY%THIS+POINTMOD (   H  R   a   POINTTYPE%SETZ+POINTMOD      Y      SETZ+POINTMOD #   ó  W   a   SETZ%THIS+POINTMOD     J  @   a   SETZ%Z+POINTMOD (     R   a   POINTTYPE%GETZ+POINTMOD    Ü  Z      GETZ+POINTMOD #   6  W   a   GETZ%THIS+POINTMOD )     ¹      POINTPTRTYPE+POINTPTRMOD -   F  _   a   POINTPTRTYPE%PTR+POINTPTRMOD 2   ¥  V   a   POINTPTRTYPE%ALLOCATE+POINTPTRMOD %   û  ]      ALLOCATE+POINTPTRMOD *   X  Z   a   ALLOCATE%THIS+POINTPTRMOD +   ²  W   a   ALLOCATE%POINT+POINTPTRMOD /   	  S   a   POINTPTRTYPE%SETID+POINTPTRMOD "   \  Z      SETID+POINTPTRMOD '   ¶  Z   a   SETID%THIS+POINTPTRMOD %     @   a   SETID%ID+POINTPTRMOD /   P  S   a   POINTPTRTYPE%GETID+POINTPTRMOD "   £  Z      GETID+POINTPTRMOD '   ý  Z   a   GETID%THIS+POINTPTRMOD .   W  R   a   POINTPTRTYPE%SETX+POINTPTRMOD !   ©  Y      SETX+POINTPTRMOD &     Z   a   SETX%THIS+POINTPTRMOD #   \  @   a   SETX%X+POINTPTRMOD .     R   a   POINTPTRTYPE%GETX+POINTPTRMOD !   î  Z      GETX+POINTPTRMOD &   H  Z   a   GETX%THIS+POINTPTRMOD .   ¢  R   a   POINTPTRTYPE%SETY+POINTPTRMOD !   ô  Y      SETY+POINTPTRMOD &   M  Z   a   SETY%THIS+POINTPTRMOD #   §  @   a   SETY%Y+POINTPTRMOD .   ç  R   a   POINTPTRTYPE%GETY+POINTPTRMOD !   9  Z      GETY+POINTPTRMOD &     Z   a   GETY%THIS+POINTPTRMOD .   í  R   a   POINTPTRTYPE%SETZ+POINTPTRMOD !   ?  Y      SETZ+POINTPTRMOD &     Z   a   SETZ%THIS+POINTPTRMOD #   ò  @   a   SETZ%Z+POINTPTRMOD .   2  R   a   POINTPTRTYPE%GETZ+POINTPTRMOD !     Z      GETZ+POINTPTRMOD &   Þ  Z   a   GETZ%THIS+POINTPTRMOD '   8  ^     ELEMENTTYPE+ELEMENTMOD .     H   a   ELEMENTTYPE%NPOINT+ELEMENTMOD ,   Þ  H   a   ELEMENTTYPE%NDOF+ELEMENTMOD -   &  ¦   a   ELEMENTTYPE%POINT+ELEMENTMOD 2   Ì  g   a   ELEMENTTYPE%INTEGRATOR+ELEMENTMOD 3   3  g      INTEGRATORPTRTYPE+INTEGRATORPTRMOD 7     d   a   INTEGRATORPTRTYPE%PTR+INTEGRATORPTRMOD -   þ  ñ      INTEGRATORTYPE+INTEGRATORMOD 8   ï  H   a   INTEGRATORTYPE%GAUSSORDER+INTEGRATORMOD 8   7   H   a   INTEGRATORTYPE%INTEGTERMS+INTEGRATORMOD 4         a   INTEGRATORTYPE%WEIGHT+INTEGRATORMOD 4   !  ¬   a   INTEGRATORTYPE%GPOINT+INTEGRATORMOD 7   ¿!  ¬   a   INTEGRATORTYPE%SHAPEFUNC+INTEGRATORMOD 8   k"  Ä   a   INTEGRATORTYPE%DSHAPEFUNC+INTEGRATORMOD 2   /#  R   a   INTEGRATORTYPE%INIT+INTEGRATORMOD #   #  l      INIT+INTEGRATORMOD (   í#  \   a   INIT%THIS+INTEGRATORMOD .   I$  @   a   INIT%GAUSSORDER+INTEGRATORMOD (   $  L   a   INIT%TYPE+INTEGRATORMOD :   Õ$  Z   a   INTEGRATORTYPE%VALUEGPOINTS+INTEGRATORMOD +   /%  \      VALUEGPOINTS+INTEGRATORMOD 0   %  \   a   VALUEGPOINTS%THIS+INTEGRATORMOD 0   ç%  L   a   VALUEGPOINTS%TYPE+INTEGRATORMOD ;   3&  T   %   INTEGRATORTYPE%GETG1D+INTEGRATORMOD=GETG1D %   &  R      GETG1D+INTEGRATORMOD *   Ù&  \   a   GETG1D%THIS+INTEGRATORMOD G   5'  Z   %   INTEGRATORTYPE%GETGTRIANGLE+INTEGRATORMOD=GETGTRIANGLE +   '  R      GETGTRIANGLE+INTEGRATORMOD 0   á'  \   a   GETGTRIANGLE%THIS+INTEGRATORMOD C   =(  X   %   INTEGRATORTYPE%GETGSQUARE+INTEGRATORMOD=GETGSQUARE )   (  R      GETGSQUARE+INTEGRATORMOD .   ç(  \   a   GETGSQUARE%THIS+INTEGRATORMOD <   C)  V   a   INTEGRATORPTRTYPE%ALLOCATE+INTEGRATORPTRMOD *   )  b      ALLOCATE+INTEGRATORPTRMOD /   û)  _   a   ALLOCATE%THIS+INTEGRATORPTRMOD 5   Z*  \   a   ALLOCATE%INTEGRATOR+INTEGRATORPTRMOD 0   ¶*  e   a   ELEMENTTYPE%MATERIAL+ELEMENTMOD /   +  g      MATERIALPTRTYPE+MATERIALPTRMOD 3   +  d   a   MATERIALPTRTYPE%PTR+MATERIALPTRMOD -   æ+  ~      THMATERIALTYPE+THMATERIALMOD :   d,  b   a   THMATERIALTYPE%MATERIALTYPE+THMATERIALMOD )   Æ,  P      MATERIALTYPE+MATERIALMOD :   -     a   THMATERIALTYPE%CONDUCTIVITY+THMATERIALMOD 2   ²-  R   a   THMATERIALTYPE%INIT+THMATERIALMOD #   .  b      INIT+THMATERIALMOD (   f.  \   a   INIT%THIS+THMATERIALMOD &   Â.  @   a   INIT%KX+THMATERIALMOD &   /  @   a   INIT%KY+THMATERIALMOD 8   B/  V   a   MATERIALPTRTYPE%ALLOCATE+MATERIALPTRMOD (   /  `      ALLOCATE+MATERIALPTRMOD -   ø/  ]   a   ALLOCATE%THIS+MATERIALPTRMOD 1   U0  \   a   ALLOCATE%MATERIAL+MATERIALPTRMOD 0   ±0  e   a   ELEMENTTYPE%GEOMETRY+ELEMENTMOD /   1  Y      GEOMETRYPTRTYPE+GEOMETRYPTRMOD 3   o1  b   a   GEOMETRYPTRTYPE%PTR+GEOMETRYPTRMOD )   Ñ1  P      GEOMETRYTYPE+GEOMETRYMOD 1   !2  W   a   ELEMENTTYPE%GETNPOINT+ELEMENTMOD %   x2  Z      GETNPOINT+ELEMENTMOD *   Ò2  Y   a   GETNPOINT%THIS+ELEMENTMOD /   +3  U   a   ELEMENTTYPE%GETNDOF+ELEMENTMOD #   3  Z      GETNDOF+ELEMENTMOD (   Ú3  Y   a   GETNDOF%THIS+ELEMENTMOD 2   34  X   a   ELEMENTTYPE%GETPOINTID+ELEMENTMOD &   4  a      GETPOINTID+ELEMENTMOD +   ì4  Y   a   GETPOINTID%THIS+ELEMENTMOD (   E5  @   a   GETPOINTID%I+ELEMENTMOD 5   5  [   a   ELEMENTTYPE%GETINTEGRATOR+ELEMENTMOD )   à5  q      GETINTEGRATOR+ELEMENTMOD .   Q6  Y   a   GETINTEGRATOR%THIS+ELEMENTMOD 3   ª6  Y   a   ELEMENTTYPE%GETMATERIAL+ELEMENTMOD '   7  o      GETMATERIAL+ELEMENTMOD ,   r7  Y   a   GETMATERIAL%THIS+ELEMENTMOD 3   Ë7  Y   a   ELEMENTTYPE%GETGEOMETRY+ELEMENTMOD '   $8  o      GETGEOMETRY+ELEMENTMOD ,   8  Y   a   GETGEOMETRY%THIS+ELEMENTMOD 1   ì8  W   a   ELEMENTTYPE%SETNPOINT+ELEMENTMOD %   C9  ^      SETNPOINT+ELEMENTMOD *   ¡9  Y   a   SETNPOINT%THIS+ELEMENTMOD ,   ú9  @   a   SETNPOINT%NPOINT+ELEMENTMOD /   ::  U   a   ELEMENTTYPE%SETNDOF+ELEMENTMOD #   :  \      SETNDOF+ELEMENTMOD (   ë:  Y   a   SETNDOF%THIS+ELEMENTMOD (   D;  @   a   SETNDOF%NDOF+ELEMENTMOD 0   ;  V   a   ELEMENTTYPE%SETPOINT+ELEMENTMOD $   Ú;  d      SETPOINT+ELEMENTMOD )   ><  Y   a   SETPOINT%THIS+ELEMENTMOD &   <  @   a   SETPOINT%I+ELEMENTMOD *   ×<  W   a   SETPOINT%POINT+ELEMENTMOD 5   .=  [   a   ELEMENTTYPE%SETINTEGRATOR+ELEMENTMOD )   =  b      SETINTEGRATOR+ELEMENTMOD .   ë=  Y   a   SETINTEGRATOR%THIS+ELEMENTMOD 4   D>  \   a   SETINTEGRATOR%INTEGRATOR+ELEMENTMOD 3    >  Y   a   ELEMENTTYPE%GETONEPOINT+ELEMENTMOD '   ù>  s      GETONEPOINT+ELEMENTMOD ,   l?  Y   a   GETONEPOINT%THIS+ELEMENTMOD )   Å?  @   a   GETONEPOINT%I+ELEMENTMOD 4   @  Z   a   ELEMENTTYPE%GETALLPOINTS+ELEMENTMOD (   _@       GETALLPOINTS+ELEMENTMOD -   qA  Y   a   GETALLPOINTS%THIS+ELEMENTMOD     ÊA  r       QUADELEMENTTYPE 0   <B  e   a   QUADELEMENTTYPE%THELEMENT2DTYPE /   ¡B  u      THELEMENT2DTYPE+THELEMENT2DMOD =   C  c   a   THELEMENT2DTYPE%ELEMENT2DTYPE+THELEMENT2DMOD +   yC  Õ      ELEMENT2DTYPE+ELEMENT2DMOD 7   ND  a   a   ELEMENT2DTYPE%ELEMENTTYPE+ELEMENT2DMOD 0   ¯D  H   a   ELEMENT2DTYPE%AREA+ELEMENT2DMOD 3   ÷D  U   a   ELEMENT2DTYPE%GETAREA+ELEMENT2DMOD %   LE  Z      GETAREA+ELEMENT2DMOD *   ¦E  [   a   GETAREA%THIS+ELEMENT2DMOD 8   F  `   a   ELEMENT2DTYPE%GETSTIFFNESS+ELEMENT2DMOD 0   aF  6     GETSTIFFNESSINTERF+ELEMENT2DMOD 5   J  [   a   GETSTIFFNESSINTERF%THIS+ELEMENT2DMOD 3   òJ  [   a   ELEMENT2DTYPE%SETAREA+ELEMENT2DMOD +   MK  R      SETAREAINTERF+ELEMENT2DMOD 0   K  [   a   SETAREAINTERF%THIS+ELEMENT2DMOD 5   úK  ]   a   ELEMENT2DTYPE%SHAPEFUNC+ELEMENT2DMOD -   WL       SHAPEFUNCINTERF+ELEMENT2DMOD 2   cN  [   a   SHAPEFUNCINTERF%THIS+ELEMENT2DMOD /   ¾N  @   a   SHAPEFUNCINTERF%X+ELEMENT2DMOD /   þN  @   a   SHAPEFUNCINTERF%Y+ELEMENT2DMOD 6   >O  ^   a   ELEMENT2DTYPE%DSHAPEFUNC+ELEMENT2DMOD .   O  ,     DSHAPEFUNCINTERF+ELEMENT2DMOD 3   ÈQ  [   a   DSHAPEFUNCINTERF%THIS+ELEMENT2DMOD 0   #R  @   a   DSHAPEFUNCINTERF%X+ELEMENT2DMOD 0   cR  @   a   DSHAPEFUNCINTERF%Y+ELEMENT2DMOD 4   £R  V   a   ELEMENT2DTYPE%JACOBIAN+ELEMENT2DMOD &   ùR  Ü      JACOBIAN+ELEMENT2DMOD +   ÕS  [   a   JACOBIAN%THIS+ELEMENT2DMOD (   0T  @   a   JACOBIAN%X+ELEMENT2DMOD (   pT  @   a   JACOBIAN%Y+ELEMENT2DMOD 7   °T  Y   a   ELEMENT2DTYPE%JACOBIANDET+ELEMENT2DMOD )   	U  h      JACOBIANDET+ELEMENT2DMOD .   qU  [   a   JACOBIANDET%THIS+ELEMENT2DMOD 2   ÌU  ´   a   JACOBIANDET%JACOBIAN+ELEMENT2DMOD <   V  Z   a   THELEMENT2DTYPE%GETSTIFFNESS+THELEMENT2DMOD ,   ÚV  ¨     GETSTIFFNESS+THELEMENT2DMOD 1   \  ]   a   GETSTIFFNESS%THIS+THELEMENT2DMOD (   ß\  U   a   QUADELEMENTTYPE%SETAREA    4]  R      SETAREA    ]  ]   a   SETAREA%THIS 