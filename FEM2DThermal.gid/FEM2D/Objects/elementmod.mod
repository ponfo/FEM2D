  řD  ž   k820309    9          19.0        m9Ü]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DThermal.gid/FEM2D/Source/Element/Element.f90 ELEMENTMOD              ELEMENTTYPE                                                     
                            @                              
                            @                              
                            @                              
                            @                              
                            @                              
                            @                              
                      @  @                               '              
      #PTR 	   #ALLOCATE 2   #SETID 6   #GETID :   #SETX =   #GETX A   #SETY D   #GETY H   #SETZ K   #GETZ O                $                              	                            #POINTTYPE 
                     @                           
     '                     #ID    #X    #Y    #Z    #INIT    #SETID    #GETID    #SETX    #GETX !   #SETY $   #GETY (   #SETZ +   #GETZ /                 D                                                              D                                            
                 D                                            
                 D                                            
   1         Ŕ    $                                              #INIT    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                                     #POINTTYPE 
             
                                                      
                                      
                
                                      
                
                                      
      1         Ŕ    $                                              #SETID    #         @     @                                                #THIS    #ID              
                                                     #POINTTYPE 
             
                                            1         Ŕ    $                                             #GETID    %         @   @                                                      #THIS              
                                                     #POINTTYPE 
   1         Ŕ    $                                              #SETX    #         @     @                                                #THIS    #X               
                                                     #POINTTYPE 
             
                                       
      1         Ŕ    $                           !             	     #GETX "   %         @   @                           "                    
       #THIS #             
                                #                     #POINTTYPE 
   1         Ŕ    $                            $             
     #SETY %   #         @     @                            %                    #THIS &   #Y '             
                                &                     #POINTTYPE 
             
                                 '     
      1         Ŕ    $                           (                  #GETY )   %         @   @                           )                    
       #THIS *             
                                *                     #POINTTYPE 
   1         Ŕ    $                            +                  #SETZ ,   #         @     @                            ,                    #THIS -   #Z .             
                                -                     #POINTTYPE 
             
                                 .     
      1         Ŕ    $                           /              	    #GETZ 0   %         @   @                           0                    
       #THIS 1             
                                1                     #POINTTYPE 
   1         Ŕ    $                            2                  #ALLOCATE 3   #         @     @                            3                    #THIS 4   #POINT 5             
                                4                    #POINTPTRTYPE              
                                  5                    #POINTTYPE 
   1         Ŕ    $                            6                  #SETID 7   #         @     @                            7                    #THIS 8   #ID 9             
                                8                    #POINTPTRTYPE              
                                 9           1         Ŕ    $                          :                  #GETID ;   %         @   @                          ;                           #THIS <             
                                <                    #POINTPTRTYPE    1         Ŕ    $                            =                  #SETX >   #         @     @                            >                    #THIS ?   #X @             
                                ?                    #POINTPTRTYPE              
                                 @     
      1         Ŕ    $                           A                  #GETX B   %         @   @                           B                    
       #THIS C             
                                C                    #POINTPTRTYPE    1         Ŕ    $                            D                  #SETY E   #         @     @                            E                    #THIS F   #Y G             
                                F                    #POINTPTRTYPE              
                                 G     
      1         Ŕ    $                           H                  #GETY I   %         @   @                           I                    
       #THIS J             
                                J                    #POINTPTRTYPE    1         Ŕ    $                            K             	     #SETZ L   #         @     @                            L                    #THIS M   #Z N             
                                M                    #POINTPTRTYPE              
                                 N     
      1         Ŕ    $                           O             
 	    #GETZ P   %         @   @                           P                    
       #THIS Q             
                                Q                    #POINTPTRTYPE                  @  @                          R     '                    #PTR S   #ALLOCATE m                $                             S                          #INTEGRATORTYPE T                     @               @           T     '                   #GAUSSORDER U   #INTEGTERMS V   #WEIGHT W   #GPOINT X   #SHAPEFUNC Y   #DSHAPEFUNC Z   #INIT [   #VALUEGPOINTS `   #GETG1D d   #GETGTRIANGLE g   #GETGSQUARE j                 $                             U                                 $                             V                             $                             W                             
            &                                                      $                             X            P                 
            &                   &                                                       $                             Y            °                 
            &                   &                                                       $                             Z                            
            &                   &                   &                                           1         Ŕ    $                            [                  #INIT \   #         @     @                            \                    #THIS ]   #GAUSSORDER ^   #TYPE _             
                                ]                   #INTEGRATORTYPE T             
                                 ^                     
                                _                    1 1         Ŕ    $                            `                  #VALUEGPOINTS a   #         @     @                            a                    #THIS b   #TYPE c             
                                b                   #INTEGRATORTYPE T             
                                c                    1 1         Ŕ    D                            d             	     #GETG1D e   #         @     @                            e                    #THIS f             
                                f                   #INTEGRATORTYPE T   1         Ŕ    D                            g             
     #GETGTRIANGLE h   #         @     @                            h                    #THIS i             
                                i                   #INTEGRATORTYPE T   1         Ŕ    D                            j                  #GETGSQUARE k   #         @     @                            k                    #THIS l             
                                l                   #INTEGRATORTYPE T   1         Ŕ    $                            m                  #ALLOCATE n   #         @     @                            n                    #THIS o   #INTEGRATOR p             
                                o                    #INTEGRATORPTRTYPE R             
                                  p                  #INTEGRATORTYPE T                 @  @                          q     '                    #PTR r   #ALLOCATE |                $                             r                           #THMATERIALTYPE s                 @  @                         s     '                    #MATERIALTYPE t   #CONDUCTIVITY v   #INIT w                 $                              t                            #MATERIALTYPE u                 @  @                          u     '                                      $                             v                              
  p          p            p                          1         Ŕ    $                            w                  #INIT x   #         @     @                            x                    #THIS y   #KX z   #KY {             
                                y                    #THMATERIALTYPE s             
                                 z     
                
                                 {     
      1         Ŕ    $                            |                  #ALLOCATE }   #         @     @                            }                    #THIS ~   #MATERIAL                                              ~                    #MATERIALPTRTYPE q             
                                                     #THMATERIALTYPE s                 @  @                               '                    #PTR                 $                                                         #GEOMETRYTYPE                   @  @                               '                                          @               @               'Ř                   #ID    #NPOINT    #NDOF    #POINT    #INTEGRATOR    #MATERIAL    #GEOMETRY    #GETID    #GETNPOINT    #GETNDOF    #GETPOINTID    #GETINTEGRATOR    #GETMATERIAL    #GETGEOMETRY    #SETID Ą   #SETNPOINT Ľ   #SETNDOF Š   #SETPOINT ­   #SETINTEGRATOR ˛   #GETONEPOINT ś   #GETALLPOINTS ş                 $                                                             $                                                             $                                                           $                                                              #POINTPTRTYPE              &                                                         $                                          X              #INTEGRATORPTRTYPE R                 $                                          Ř              #MATERIALPTRTYPE q                 $                                          X             #GEOMETRYPTRTYPE    1         Ŕ    $                                             #GETID    %         @   @                                                       #THIS              
                                     Ř              #ELEMENTTYPE    1         Ŕ    $                                        	     #GETNPOINT    %         @   @                                                       #THIS              
                                     Ř              #ELEMENTTYPE    1         Ŕ    $                                        
     #GETNDOF    %         @   @                                                       #THIS              
                                     Ř              #ELEMENTTYPE    1         Ŕ    $                                             #GETPOINTID    %         @   @                                                       #THIS    #I              
D @                                   Ř              #ELEMENTTYPE              
                                            1         Ŕ    $                                             #GETINTEGRATOR    &         @   @                                                       #THIS    #INTEGRATORPTRTYPE R             
                                     Ř              #ELEMENTTYPE    1         Ŕ    $                                             #GETMATERIAL    &         @   @                                                       #THIS    #MATERIALPTRTYPE q             
                                     Ř              #ELEMENTTYPE    1         Ŕ    $                                             #GETGEOMETRY    &         @   @                                                       #THIS     #GEOMETRYPTRTYPE              
                                      Ř              #ELEMENTTYPE    1         Ŕ    $                            Ą                  #SETID ˘   #         @     @                             ˘                    #THIS Ł   #ID ¤             
D                                Ł     Ř              #ELEMENTTYPE              
                                 ¤           1         Ŕ    $                            Ľ              	    #SETNPOINT Ś   #         @     @                             Ś                    #THIS §   #NPOINT ¨             
D                                §     Ř              #ELEMENTTYPE              
                                 ¨           1         Ŕ    $                            Š              
    #SETNDOF Ş   #         @     @                             Ş                    #THIS Ť   #NDOF Ź             
D                                Ť     Ř              #ELEMENTTYPE              
                                 Ź           1         Ŕ    $                            ­                  #SETPOINT Ž   #         @     @                             Ž                    #THIS Ż   #I °   #POINT ą             
D                                Ż     Ř              #ELEMENTTYPE              
                                 °                     
                                 ą                    #POINTTYPE 
   1         Ŕ    $                            ˛                  #SETINTEGRATOR ł   #         @     @                             ł                    #THIS ´   #INTEGRATOR ľ             
D                                ´     Ř              #ELEMENTTYPE              
                                  ľ                  #INTEGRATORTYPE T   1         Ŕ    $                           ś                  #GETONEPOINT ˇ   &         @   @                            ˇ                           #THIS ¸   #I š   #POINTPTRTYPE              
                                ¸     Ř              #ELEMENTTYPE              
                                 š           1         Ŕ    $                           ş                  #GETALLPOINTS ť   )        `   @                             ť                                         #THIS ź   #POINTPTRTYPE    p          5 8 O#ELEMENTTYPE     p        U            5 8 O#ELEMENTTYPE     p        U                                   
                                ź     Ř              #ELEMENTTYPE           u      fn#fn          b   uapp(ELEMENTMOD    1  @   J  TOOLS    q  @   J  POINTMOD    ą  @   J  POINTPTRMOD    ń  @   J  INTEGRATORMOD !   1  @   J  INTEGRATORPTRMOD    q  @   J  MATERIALPTRMOD    ą  @   J  GEOMETRYPTRMOD )   ń  š      POINTPTRTYPE+POINTPTRMOD -   Ş  _   a   POINTPTRTYPE%PTR+POINTPTRMOD #   	  É       POINTTYPE+POINTMOD )   Ň  H   %   POINTTYPE%ID+POINTMOD=ID '     H   %   POINTTYPE%X+POINTMOD=X '   b  H   %   POINTTYPE%Y+POINTMOD=Y '   Ş  H   %   POINTTYPE%Z+POINTMOD=Z (   ň  R   a   POINTTYPE%INIT+POINTMOD    D  o      INIT+POINTMOD #   ł  W   a   INIT%THIS+POINTMOD !   
  @   a   INIT%ID+POINTMOD     J  @   a   INIT%X+POINTMOD       @   a   INIT%Y+POINTMOD     Ę  @   a   INIT%Z+POINTMOD )   
  S   a   POINTTYPE%SETID+POINTMOD    ]  Z      SETID+POINTMOD $   ˇ  W   a   SETID%THIS+POINTMOD "   	  @   a   SETID%ID+POINTMOD )   N	  S   a   POINTTYPE%GETID+POINTMOD    Ą	  Z      GETID+POINTMOD $   ű	  W   a   GETID%THIS+POINTMOD (   R
  R   a   POINTTYPE%SETX+POINTMOD    ¤
  Y      SETX+POINTMOD #   ý
  W   a   SETX%THIS+POINTMOD     T  @   a   SETX%X+POINTMOD (     R   a   POINTTYPE%GETX+POINTMOD    ć  Z      GETX+POINTMOD #   @  W   a   GETX%THIS+POINTMOD (     R   a   POINTTYPE%SETY+POINTMOD    é  Y      SETY+POINTMOD #   B  W   a   SETY%THIS+POINTMOD       @   a   SETY%Y+POINTMOD (   Ů  R   a   POINTTYPE%GETY+POINTMOD    +  Z      GETY+POINTMOD #     W   a   GETY%THIS+POINTMOD (   Ü  R   a   POINTTYPE%SETZ+POINTMOD    .  Y      SETZ+POINTMOD #     W   a   SETZ%THIS+POINTMOD     Ţ  @   a   SETZ%Z+POINTMOD (     R   a   POINTTYPE%GETZ+POINTMOD    p  Z      GETZ+POINTMOD #   Ę  W   a   GETZ%THIS+POINTMOD 2   !  V   a   POINTPTRTYPE%ALLOCATE+POINTPTRMOD %   w  ]      ALLOCATE+POINTPTRMOD *   Ô  Z   a   ALLOCATE%THIS+POINTPTRMOD +   .  W   a   ALLOCATE%POINT+POINTPTRMOD /     S   a   POINTPTRTYPE%SETID+POINTPTRMOD "   Ř  Z      SETID+POINTPTRMOD '   2  Z   a   SETID%THIS+POINTPTRMOD %     @   a   SETID%ID+POINTPTRMOD /   Ě  S   a   POINTPTRTYPE%GETID+POINTPTRMOD "     Z      GETID+POINTPTRMOD '   y  Z   a   GETID%THIS+POINTPTRMOD .   Ó  R   a   POINTPTRTYPE%SETX+POINTPTRMOD !   %  Y      SETX+POINTPTRMOD &   ~  Z   a   SETX%THIS+POINTPTRMOD #   Ř  @   a   SETX%X+POINTPTRMOD .     R   a   POINTPTRTYPE%GETX+POINTPTRMOD !   j  Z      GETX+POINTPTRMOD &   Ä  Z   a   GETX%THIS+POINTPTRMOD .     R   a   POINTPTRTYPE%SETY+POINTPTRMOD !   p  Y      SETY+POINTPTRMOD &   É  Z   a   SETY%THIS+POINTPTRMOD #   #  @   a   SETY%Y+POINTPTRMOD .   c  R   a   POINTPTRTYPE%GETY+POINTPTRMOD !   ľ  Z      GETY+POINTPTRMOD &     Z   a   GETY%THIS+POINTPTRMOD .   i  R   a   POINTPTRTYPE%SETZ+POINTPTRMOD !   ť  Y      SETZ+POINTPTRMOD &     Z   a   SETZ%THIS+POINTPTRMOD #   n  @   a   SETZ%Z+POINTPTRMOD .   Ž  R   a   POINTPTRTYPE%GETZ+POINTPTRMOD !      Z      GETZ+POINTPTRMOD &   Z  Z   a   GETZ%THIS+POINTPTRMOD 3   ´  g      INTEGRATORPTRTYPE+INTEGRATORPTRMOD 7     d   a   INTEGRATORPTRTYPE%PTR+INTEGRATORPTRMOD -     ń       INTEGRATORTYPE+INTEGRATORMOD 8   p  H   a   INTEGRATORTYPE%GAUSSORDER+INTEGRATORMOD 8   ¸  H   a   INTEGRATORTYPE%INTEGTERMS+INTEGRATORMOD 4         a   INTEGRATORTYPE%WEIGHT+INTEGRATORMOD 4     Ź   a   INTEGRATORTYPE%GPOINT+INTEGRATORMOD 7   @  Ź   a   INTEGRATORTYPE%SHAPEFUNC+INTEGRATORMOD 8   ě  Ä   a   INTEGRATORTYPE%DSHAPEFUNC+INTEGRATORMOD 2   °   R   a   INTEGRATORTYPE%INIT+INTEGRATORMOD #   !  l      INIT+INTEGRATORMOD (   n!  \   a   INIT%THIS+INTEGRATORMOD .   Ę!  @   a   INIT%GAUSSORDER+INTEGRATORMOD (   
"  L   a   INIT%TYPE+INTEGRATORMOD :   V"  Z   a   INTEGRATORTYPE%VALUEGPOINTS+INTEGRATORMOD +   °"  \      VALUEGPOINTS+INTEGRATORMOD 0   #  \   a   VALUEGPOINTS%THIS+INTEGRATORMOD 0   h#  L   a   VALUEGPOINTS%TYPE+INTEGRATORMOD ;   ´#  T   %   INTEGRATORTYPE%GETG1D+INTEGRATORMOD=GETG1D %   $  R      GETG1D+INTEGRATORMOD *   Z$  \   a   GETG1D%THIS+INTEGRATORMOD G   ś$  Z   %   INTEGRATORTYPE%GETGTRIANGLE+INTEGRATORMOD=GETGTRIANGLE +   %  R      GETGTRIANGLE+INTEGRATORMOD 0   b%  \   a   GETGTRIANGLE%THIS+INTEGRATORMOD C   ž%  X   %   INTEGRATORTYPE%GETGSQUARE+INTEGRATORMOD=GETGSQUARE )   &  R      GETGSQUARE+INTEGRATORMOD .   h&  \   a   GETGSQUARE%THIS+INTEGRATORMOD <   Ä&  V   a   INTEGRATORPTRTYPE%ALLOCATE+INTEGRATORPTRMOD *   '  b      ALLOCATE+INTEGRATORPTRMOD /   |'  _   a   ALLOCATE%THIS+INTEGRATORPTRMOD 5   Ű'  \   a   ALLOCATE%INTEGRATOR+INTEGRATORPTRMOD /   7(  g      MATERIALPTRTYPE+MATERIALPTRMOD 3   (  d   a   MATERIALPTRTYPE%PTR+MATERIALPTRMOD -   )  ~      THMATERIALTYPE+THMATERIALMOD :   )  b   a   THMATERIALTYPE%MATERIALTYPE+THMATERIALMOD )   â)  P      MATERIALTYPE+MATERIALMOD :   2*     a   THMATERIALTYPE%CONDUCTIVITY+THMATERIALMOD 2   Î*  R   a   THMATERIALTYPE%INIT+THMATERIALMOD #    +  b      INIT+THMATERIALMOD (   +  \   a   INIT%THIS+THMATERIALMOD &   Ţ+  @   a   INIT%KX+THMATERIALMOD &   ,  @   a   INIT%KY+THMATERIALMOD 8   ^,  V   a   MATERIALPTRTYPE%ALLOCATE+MATERIALPTRMOD (   ´,  `      ALLOCATE+MATERIALPTRMOD -   -  ]   a   ALLOCATE%THIS+MATERIALPTRMOD 1   q-  \   a   ALLOCATE%MATERIAL+MATERIALPTRMOD /   Í-  Y      GEOMETRYPTRTYPE+GEOMETRYPTRMOD 3   &.  b   a   GEOMETRYPTRTYPE%PTR+GEOMETRYPTRMOD )   .  P      GEOMETRYTYPE+GEOMETRYMOD    Ř.  |      ELEMENTTYPE    T0  H   a   ELEMENTTYPE%ID #   0  H   a   ELEMENTTYPE%NPOINT !   ä0  H   a   ELEMENTTYPE%NDOF "   ,1  Ś   a   ELEMENTTYPE%POINT '   Ň1  g   a   ELEMENTTYPE%INTEGRATOR %   92  e   a   ELEMENTTYPE%MATERIAL %   2  e   a   ELEMENTTYPE%GEOMETRY "   3  S   a   ELEMENTTYPE%GETID    V3  Z      GETID    °3  Y   a   GETID%THIS &   	4  W   a   ELEMENTTYPE%GETNPOINT    `4  Z      GETNPOINT    ş4  Y   a   GETNPOINT%THIS $   5  U   a   ELEMENTTYPE%GETNDOF    h5  Z      GETNDOF    Â5  Y   a   GETNDOF%THIS '   6  X   a   ELEMENTTYPE%GETPOINTID    s6  a      GETPOINTID     Ô6  Y   a   GETPOINTID%THIS    -7  @   a   GETPOINTID%I *   m7  [   a   ELEMENTTYPE%GETINTEGRATOR    Č7  q      GETINTEGRATOR #   98  Y   a   GETINTEGRATOR%THIS (   8  Y   a   ELEMENTTYPE%GETMATERIAL    ë8  o      GETMATERIAL !   Z9  Y   a   GETMATERIAL%THIS (   ł9  Y   a   ELEMENTTYPE%GETGEOMETRY    :  o      GETGEOMETRY !   {:  Y   a   GETGEOMETRY%THIS "   Ô:  S   a   ELEMENTTYPE%SETID    ';  Z      SETID    ;  Y   a   SETID%THIS    Ú;  @   a   SETID%ID &   <  W   a   ELEMENTTYPE%SETNPOINT    q<  ^      SETNPOINT    Ď<  Y   a   SETNPOINT%THIS !   (=  @   a   SETNPOINT%NPOINT $   h=  U   a   ELEMENTTYPE%SETNDOF    ˝=  \      SETNDOF    >  Y   a   SETNDOF%THIS    r>  @   a   SETNDOF%NDOF %   ˛>  V   a   ELEMENTTYPE%SETPOINT    ?  d      SETPOINT    l?  Y   a   SETPOINT%THIS    Ĺ?  @   a   SETPOINT%I    @  W   a   SETPOINT%POINT *   \@  [   a   ELEMENTTYPE%SETINTEGRATOR    ˇ@  b      SETINTEGRATOR #   A  Y   a   SETINTEGRATOR%THIS )   rA  \   a   SETINTEGRATOR%INTEGRATOR (   ÎA  Y   a   ELEMENTTYPE%GETONEPOINT    'B  s      GETONEPOINT !   B  Y   a   GETONEPOINT%THIS    óB  @   a   GETONEPOINT%I )   3C  Z   a   ELEMENTTYPE%GETALLPOINTS    C       GETALLPOINTS "   D  Y   a   GETALLPOINTS%THIS 