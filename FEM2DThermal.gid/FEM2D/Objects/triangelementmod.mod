  ë]  Ý   k820309    9          19.0        (áÑ]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DThermal.gid/FEM2D/Source/TriangElement.f90 TRIANGELEMENTMOD              TRIANGELEMENTTYPE                                                     
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
                                ²     Ð              #ELEMENTTYPE P                     @                        ³     'Ø                   #THELEMENT2DTYPE ´   #SETAREA Ù                 $                              ´     Ø                     #THELEMENT2DTYPE µ                 @  @                        µ     'Ø                   #ELEMENT2DTYPE ¶   #GETSTIFFNESS Ö                 $                              ¶     Ø                     #ELEMENT2DTYPE ·                  @  @                        ·     'Ø             	      #ELEMENTTYPE ¸   #AREA ¹   #GETAREA º   #GETSTIFFNESS ½   #SETAREA À   #SHAPEFUNC Ã   #DSHAPEFUNC È   #JACOBIAN Í   #JACOBIANDET Ò                 $                              ¸     Ð                     #ELEMENTTYPE P                 $                             ¹     Ð         
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
D @                              Û     Ø              #TRIANGELEMENTTYPE ³          y      fn#fn &     "   b   uapp(TRIANGELEMENTMOD    ;  @   J  TOOLS    {  @   J  POINTMOD    »  @   J  POINTPTRMOD    û  @   J  ELEMENT2DMOD    ;  @   J  THELEMENT2DMOD #   {  É       POINTTYPE+POINTMOD )   D  H   %   POINTTYPE%ID+POINTMOD=ID '     H   %   POINTTYPE%X+POINTMOD=X '   Ô  H   %   POINTTYPE%Y+POINTMOD=Y '     H   %   POINTTYPE%Z+POINTMOD=Z (   d  R   a   POINTTYPE%INIT+POINTMOD    ¶  o      INIT+POINTMOD #   %  W   a   INIT%THIS+POINTMOD !   |  @   a   INIT%ID+POINTMOD     ¼  @   a   INIT%X+POINTMOD     ü  @   a   INIT%Y+POINTMOD     <  @   a   INIT%Z+POINTMOD )   |  S   a   POINTTYPE%SETID+POINTMOD    Ï  Z      SETID+POINTMOD $   )  W   a   SETID%THIS+POINTMOD "     @   a   SETID%ID+POINTMOD )   À  S   a   POINTTYPE%GETID+POINTMOD      Z      GETID+POINTMOD $   m  W   a   GETID%THIS+POINTMOD (   Ä  R   a   POINTTYPE%SETX+POINTMOD    	  Y      SETX+POINTMOD #   o	  W   a   SETX%THIS+POINTMOD     Æ	  @   a   SETX%X+POINTMOD (   
  R   a   POINTTYPE%GETX+POINTMOD    X
  Z      GETX+POINTMOD #   ²
  W   a   GETX%THIS+POINTMOD (   	  R   a   POINTTYPE%SETY+POINTMOD    [  Y      SETY+POINTMOD #   ´  W   a   SETY%THIS+POINTMOD       @   a   SETY%Y+POINTMOD (   K  R   a   POINTTYPE%GETY+POINTMOD      Z      GETY+POINTMOD #   ÷  W   a   GETY%THIS+POINTMOD (   N  R   a   POINTTYPE%SETZ+POINTMOD       Y      SETZ+POINTMOD #   ù  W   a   SETZ%THIS+POINTMOD     P  @   a   SETZ%Z+POINTMOD (     R   a   POINTTYPE%GETZ+POINTMOD    â  Z      GETZ+POINTMOD #   <  W   a   GETZ%THIS+POINTMOD )     ¹      POINTPTRTYPE+POINTPTRMOD -   L  _   a   POINTPTRTYPE%PTR+POINTPTRMOD 2   «  V   a   POINTPTRTYPE%ALLOCATE+POINTPTRMOD %     ]      ALLOCATE+POINTPTRMOD *   ^  Z   a   ALLOCATE%THIS+POINTPTRMOD +   ¸  W   a   ALLOCATE%POINT+POINTPTRMOD /     S   a   POINTPTRTYPE%SETID+POINTPTRMOD "   b  Z      SETID+POINTPTRMOD '   ¼  Z   a   SETID%THIS+POINTPTRMOD %     @   a   SETID%ID+POINTPTRMOD /   V  S   a   POINTPTRTYPE%GETID+POINTPTRMOD "   ©  Z      GETID+POINTPTRMOD '     Z   a   GETID%THIS+POINTPTRMOD .   ]  R   a   POINTPTRTYPE%SETX+POINTPTRMOD !   ¯  Y      SETX+POINTPTRMOD &     Z   a   SETX%THIS+POINTPTRMOD #   b  @   a   SETX%X+POINTPTRMOD .   ¢  R   a   POINTPTRTYPE%GETX+POINTPTRMOD !   ô  Z      GETX+POINTPTRMOD &   N  Z   a   GETX%THIS+POINTPTRMOD .   ¨  R   a   POINTPTRTYPE%SETY+POINTPTRMOD !   ú  Y      SETY+POINTPTRMOD &   S  Z   a   SETY%THIS+POINTPTRMOD #   ­  @   a   SETY%Y+POINTPTRMOD .   í  R   a   POINTPTRTYPE%GETY+POINTPTRMOD !   ?  Z      GETY+POINTPTRMOD &     Z   a   GETY%THIS+POINTPTRMOD .   ó  R   a   POINTPTRTYPE%SETZ+POINTPTRMOD !   E  Y      SETZ+POINTPTRMOD &     Z   a   SETZ%THIS+POINTPTRMOD #   ø  @   a   SETZ%Z+POINTPTRMOD .   8  R   a   POINTPTRTYPE%GETZ+POINTPTRMOD !     Z      GETZ+POINTPTRMOD &   ä  Z   a   GETZ%THIS+POINTPTRMOD '   >  ^     ELEMENTTYPE+ELEMENTMOD .     H   a   ELEMENTTYPE%NPOINT+ELEMENTMOD ,   ä  H   a   ELEMENTTYPE%NDOF+ELEMENTMOD -   ,  ¦   a   ELEMENTTYPE%POINT+ELEMENTMOD 2   Ò  g   a   ELEMENTTYPE%INTEGRATOR+ELEMENTMOD 3   9  g      INTEGRATORPTRTYPE+INTEGRATORPTRMOD 7      d   a   INTEGRATORPTRTYPE%PTR+INTEGRATORPTRMOD -     ñ      INTEGRATORTYPE+INTEGRATORMOD 8   õ  H   a   INTEGRATORTYPE%GAUSSORDER+INTEGRATORMOD 8   =   H   a   INTEGRATORTYPE%INTEGTERMS+INTEGRATORMOD 4         a   INTEGRATORTYPE%WEIGHT+INTEGRATORMOD 4   !  ¬   a   INTEGRATORTYPE%GPOINT+INTEGRATORMOD 7   Å!  ¬   a   INTEGRATORTYPE%SHAPEFUNC+INTEGRATORMOD 8   q"  Ä   a   INTEGRATORTYPE%DSHAPEFUNC+INTEGRATORMOD 2   5#  R   a   INTEGRATORTYPE%INIT+INTEGRATORMOD #   #  l      INIT+INTEGRATORMOD (   ó#  \   a   INIT%THIS+INTEGRATORMOD .   O$  @   a   INIT%GAUSSORDER+INTEGRATORMOD (   $  L   a   INIT%TYPE+INTEGRATORMOD :   Û$  Z   a   INTEGRATORTYPE%VALUEGPOINTS+INTEGRATORMOD +   5%  \      VALUEGPOINTS+INTEGRATORMOD 0   %  \   a   VALUEGPOINTS%THIS+INTEGRATORMOD 0   í%  L   a   VALUEGPOINTS%TYPE+INTEGRATORMOD ;   9&  T   %   INTEGRATORTYPE%GETG1D+INTEGRATORMOD=GETG1D %   &  R      GETG1D+INTEGRATORMOD *   ß&  \   a   GETG1D%THIS+INTEGRATORMOD G   ;'  Z   %   INTEGRATORTYPE%GETGTRIANGLE+INTEGRATORMOD=GETGTRIANGLE +   '  R      GETGTRIANGLE+INTEGRATORMOD 0   ç'  \   a   GETGTRIANGLE%THIS+INTEGRATORMOD C   C(  X   %   INTEGRATORTYPE%GETGSQUARE+INTEGRATORMOD=GETGSQUARE )   (  R      GETGSQUARE+INTEGRATORMOD .   í(  \   a   GETGSQUARE%THIS+INTEGRATORMOD <   I)  V   a   INTEGRATORPTRTYPE%ALLOCATE+INTEGRATORPTRMOD *   )  b      ALLOCATE+INTEGRATORPTRMOD /   *  _   a   ALLOCATE%THIS+INTEGRATORPTRMOD 5   `*  \   a   ALLOCATE%INTEGRATOR+INTEGRATORPTRMOD 0   ¼*  e   a   ELEMENTTYPE%MATERIAL+ELEMENTMOD /   !+  g      MATERIALPTRTYPE+MATERIALPTRMOD 3   +  d   a   MATERIALPTRTYPE%PTR+MATERIALPTRMOD -   ì+  ~      THMATERIALTYPE+THMATERIALMOD :   j,  b   a   THMATERIALTYPE%MATERIALTYPE+THMATERIALMOD )   Ì,  P      MATERIALTYPE+MATERIALMOD :   -     a   THMATERIALTYPE%CONDUCTIVITY+THMATERIALMOD 2   ¸-  R   a   THMATERIALTYPE%INIT+THMATERIALMOD #   
.  b      INIT+THMATERIALMOD (   l.  \   a   INIT%THIS+THMATERIALMOD &   È.  @   a   INIT%KX+THMATERIALMOD &   /  @   a   INIT%KY+THMATERIALMOD 8   H/  V   a   MATERIALPTRTYPE%ALLOCATE+MATERIALPTRMOD (   /  `      ALLOCATE+MATERIALPTRMOD -   þ/  ]   a   ALLOCATE%THIS+MATERIALPTRMOD 1   [0  \   a   ALLOCATE%MATERIAL+MATERIALPTRMOD 0   ·0  e   a   ELEMENTTYPE%GEOMETRY+ELEMENTMOD /   1  Y      GEOMETRYPTRTYPE+GEOMETRYPTRMOD 3   u1  b   a   GEOMETRYPTRTYPE%PTR+GEOMETRYPTRMOD )   ×1  P      GEOMETRYTYPE+GEOMETRYMOD 1   '2  W   a   ELEMENTTYPE%GETNPOINT+ELEMENTMOD %   ~2  Z      GETNPOINT+ELEMENTMOD *   Ø2  Y   a   GETNPOINT%THIS+ELEMENTMOD /   13  U   a   ELEMENTTYPE%GETNDOF+ELEMENTMOD #   3  Z      GETNDOF+ELEMENTMOD (   à3  Y   a   GETNDOF%THIS+ELEMENTMOD 2   94  X   a   ELEMENTTYPE%GETPOINTID+ELEMENTMOD &   4  a      GETPOINTID+ELEMENTMOD +   ò4  Y   a   GETPOINTID%THIS+ELEMENTMOD (   K5  @   a   GETPOINTID%I+ELEMENTMOD 5   5  [   a   ELEMENTTYPE%GETINTEGRATOR+ELEMENTMOD )   æ5  q      GETINTEGRATOR+ELEMENTMOD .   W6  Y   a   GETINTEGRATOR%THIS+ELEMENTMOD 3   °6  Y   a   ELEMENTTYPE%GETMATERIAL+ELEMENTMOD '   	7  o      GETMATERIAL+ELEMENTMOD ,   x7  Y   a   GETMATERIAL%THIS+ELEMENTMOD 3   Ñ7  Y   a   ELEMENTTYPE%GETGEOMETRY+ELEMENTMOD '   *8  o      GETGEOMETRY+ELEMENTMOD ,   8  Y   a   GETGEOMETRY%THIS+ELEMENTMOD 1   ò8  W   a   ELEMENTTYPE%SETNPOINT+ELEMENTMOD %   I9  ^      SETNPOINT+ELEMENTMOD *   §9  Y   a   SETNPOINT%THIS+ELEMENTMOD ,    :  @   a   SETNPOINT%NPOINT+ELEMENTMOD /   @:  U   a   ELEMENTTYPE%SETNDOF+ELEMENTMOD #   :  \      SETNDOF+ELEMENTMOD (   ñ:  Y   a   SETNDOF%THIS+ELEMENTMOD (   J;  @   a   SETNDOF%NDOF+ELEMENTMOD 0   ;  V   a   ELEMENTTYPE%SETPOINT+ELEMENTMOD $   à;  d      SETPOINT+ELEMENTMOD )   D<  Y   a   SETPOINT%THIS+ELEMENTMOD &   <  @   a   SETPOINT%I+ELEMENTMOD *   Ý<  W   a   SETPOINT%POINT+ELEMENTMOD 5   4=  [   a   ELEMENTTYPE%SETINTEGRATOR+ELEMENTMOD )   =  b      SETINTEGRATOR+ELEMENTMOD .   ñ=  Y   a   SETINTEGRATOR%THIS+ELEMENTMOD 4   J>  \   a   SETINTEGRATOR%INTEGRATOR+ELEMENTMOD 3   ¦>  Y   a   ELEMENTTYPE%GETONEPOINT+ELEMENTMOD '   ÿ>  s      GETONEPOINT+ELEMENTMOD ,   r?  Y   a   GETONEPOINT%THIS+ELEMENTMOD )   Ë?  @   a   GETONEPOINT%I+ELEMENTMOD 4   @  Z   a   ELEMENTTYPE%GETALLPOINTS+ELEMENTMOD (   e@       GETALLPOINTS+ELEMENTMOD -   wA  Y   a   GETALLPOINTS%THIS+ELEMENTMOD "   ÐA  r       TRIANGELEMENTTYPE 2   BB  e   a   TRIANGELEMENTTYPE%THELEMENT2DTYPE /   §B  u      THELEMENT2DTYPE+THELEMENT2DMOD =   C  c   a   THELEMENT2DTYPE%ELEMENT2DTYPE+THELEMENT2DMOD +   C  Õ      ELEMENT2DTYPE+ELEMENT2DMOD 7   TD  a   a   ELEMENT2DTYPE%ELEMENTTYPE+ELEMENT2DMOD 0   µD  H   a   ELEMENT2DTYPE%AREA+ELEMENT2DMOD 3   ýD  U   a   ELEMENT2DTYPE%GETAREA+ELEMENT2DMOD %   RE  Z      GETAREA+ELEMENT2DMOD *   ¬E  [   a   GETAREA%THIS+ELEMENT2DMOD 8   F  `   a   ELEMENT2DTYPE%GETSTIFFNESS+ELEMENT2DMOD 0   gF  6     GETSTIFFNESSINTERF+ELEMENT2DMOD 5   J  [   a   GETSTIFFNESSINTERF%THIS+ELEMENT2DMOD 3   øJ  [   a   ELEMENT2DTYPE%SETAREA+ELEMENT2DMOD +   SK  R      SETAREAINTERF+ELEMENT2DMOD 0   ¥K  [   a   SETAREAINTERF%THIS+ELEMENT2DMOD 5    L  ]   a   ELEMENT2DTYPE%SHAPEFUNC+ELEMENT2DMOD -   ]L       SHAPEFUNCINTERF+ELEMENT2DMOD 2   iN  [   a   SHAPEFUNCINTERF%THIS+ELEMENT2DMOD /   ÄN  @   a   SHAPEFUNCINTERF%X+ELEMENT2DMOD /   O  @   a   SHAPEFUNCINTERF%Y+ELEMENT2DMOD 6   DO  ^   a   ELEMENT2DTYPE%DSHAPEFUNC+ELEMENT2DMOD .   ¢O  ,     DSHAPEFUNCINTERF+ELEMENT2DMOD 3   ÎQ  [   a   DSHAPEFUNCINTERF%THIS+ELEMENT2DMOD 0   )R  @   a   DSHAPEFUNCINTERF%X+ELEMENT2DMOD 0   iR  @   a   DSHAPEFUNCINTERF%Y+ELEMENT2DMOD 4   ©R  V   a   ELEMENT2DTYPE%JACOBIAN+ELEMENT2DMOD &   ÿR  Ü      JACOBIAN+ELEMENT2DMOD +   ÛS  [   a   JACOBIAN%THIS+ELEMENT2DMOD (   6T  @   a   JACOBIAN%X+ELEMENT2DMOD (   vT  @   a   JACOBIAN%Y+ELEMENT2DMOD 7   ¶T  Y   a   ELEMENT2DTYPE%JACOBIANDET+ELEMENT2DMOD )   U  h      JACOBIANDET+ELEMENT2DMOD .   wU  [   a   JACOBIANDET%THIS+ELEMENT2DMOD 2   ÒU  ´   a   JACOBIANDET%JACOBIAN+ELEMENT2DMOD <   V  Z   a   THELEMENT2DTYPE%GETSTIFFNESS+THELEMENT2DMOD ,   àV  ¨     GETSTIFFNESS+THELEMENT2DMOD 1   \  ]   a   GETSTIFFNESS%THIS+THELEMENT2DMOD *   å\  U   a   TRIANGELEMENTTYPE%SETAREA    :]  R      SETAREA    ]  _   a   SETAREA%THIS 