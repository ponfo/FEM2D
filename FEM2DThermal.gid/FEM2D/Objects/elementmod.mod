  7B  ¶   k820309    9          19.0        Ï]                                                                                                          
       /home/facundo/Documents/FEM2DIUA_Thermal/FEM2DThermal.gid/FEM2D/Source/Element.f90 ELEMENTMOD              ELEMENTTYPE                                                     
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
   1         À    $                                              #INIT    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                                     #POINTTYPE 
             
                                                      
                                      
                
                                      
                
                                      
      1         À    $                                              #SETID    #         @     @                                                #THIS    #ID              
                                                     #POINTTYPE 
             
                                            1         À    $                                             #GETID    %         @   @                                                      #THIS              
                                                     #POINTTYPE 
   1         À    $                                              #SETX    #         @     @                                                #THIS    #X               
                                                     #POINTTYPE 
             
                                       
      1         À    $                           !             	     #GETX "   %         @   @                           "                    
       #THIS #             
                                #                     #POINTTYPE 
   1         À    $                            $             
     #SETY %   #         @     @                            %                    #THIS &   #Y '             
                                &                     #POINTTYPE 
             
                                 '     
      1         À    $                           (                  #GETY )   %         @   @                           )                    
       #THIS *             
                                *                     #POINTTYPE 
   1         À    $                            +                  #SETZ ,   #         @     @                            ,                    #THIS -   #Z .             
                                -                     #POINTTYPE 
             
                                 .     
      1         À    $                           /              	    #GETZ 0   %         @   @                           0                    
       #THIS 1             
                                1                     #POINTTYPE 
   1         À    $                            2                  #ALLOCATE 3   #         @     @                            3                    #THIS 4   #POINT 5             
                                4                    #POINTPTRTYPE              
                                  5                    #POINTTYPE 
   1         À    $                            6                  #SETID 7   #         @     @                            7                    #THIS 8   #ID 9             
                                8                    #POINTPTRTYPE              
                                 9           1         À    $                          :                  #GETID ;   %         @   @                          ;                           #THIS <             
                                <                    #POINTPTRTYPE    1         À    $                            =                  #SETX >   #         @     @                            >                    #THIS ?   #X @             
                                ?                    #POINTPTRTYPE              
                                 @     
      1         À    $                           A                  #GETX B   %         @   @                           B                    
       #THIS C             
                                C                    #POINTPTRTYPE    1         À    $                            D                  #SETY E   #         @     @                            E                    #THIS F   #Y G             
                                F                    #POINTPTRTYPE              
                                 G     
      1         À    $                           H                  #GETY I   %         @   @                           I                    
       #THIS J             
                                J                    #POINTPTRTYPE    1         À    $                            K             	     #SETZ L   #         @     @                            L                    #THIS M   #Z N             
                                M                    #POINTPTRTYPE              
                                 N     
      1         À    $                           O             
 	    #GETZ P   %         @   @                           P                    
       #THIS Q             
                                Q                    #POINTPTRTYPE                  @  @                          R     '                    #PTR S   #ALLOCATE m                $                             S                          #INTEGRATORTYPE T                     @               @           T     '                   #GAUSSORDER U   #INTEGTERMS V   #WEIGHT W   #GPOINT X   #SHAPEFUNC Y   #DSHAPEFUNC Z   #INIT [   #VALUEGPOINTS `   #GETG1D d   #GETGTRIANGLE g   #GETGSQUARE j                 $                             U                                 $                             V                             $                             W                             
            &                                                      $                             X            P                 
            &                   &                                                       $                             Y            °                 
            &                   &                                                       $                             Z                            
            &                   &                   &                                           1         À    $                            [                  #INIT \   #         @     @                            \                    #THIS ]   #GAUSSORDER ^   #TYPE _             
                                ]                   #INTEGRATORTYPE T             
                                 ^                     
                                _                    1 1         À    $                            `                  #VALUEGPOINTS a   #         @     @                            a                    #THIS b   #TYPE c             
                                b                   #INTEGRATORTYPE T             
                                c                    1 1         À    D                            d             	     #GETG1D e   #         @     @                            e                    #THIS f             
                                f                   #INTEGRATORTYPE T   1         À    D                            g             
     #GETGTRIANGLE h   #         @     @                            h                    #THIS i             
                                i                   #INTEGRATORTYPE T   1         À    D                            j                  #GETGSQUARE k   #         @     @                            k                    #THIS l             
                                l                   #INTEGRATORTYPE T   1         À    $                            m                  #ALLOCATE n   #         @     @                            n                    #THIS o   #INTEGRATOR p             
                                o                    #INTEGRATORPTRTYPE R             
                                  p                  #INTEGRATORTYPE T                 @  @                          q     '                    #PTR r   #ALLOCATE |                $                             r                           #THMATERIALTYPE s                 @  @                         s     '                    #MATERIALTYPE t   #CONDUCTIVITY v   #INIT w                 $                              t                            #MATERIALTYPE u                 @  @                          u     '                                      $                             v                              
  p          p            p                          1         À    $                            w                  #INIT x   #         @     @                            x                    #THIS y   #KX z   #KY {             
                                y                    #THMATERIALTYPE s             
                                 z     
                
                                 {     
      1         À    $                            |                  #ALLOCATE }   #         @     @                            }                    #THIS ~   #MATERIAL                                              ~                    #MATERIALPTRTYPE q             
                                                     #THMATERIALTYPE s                 @  @                               '                    #PTR                 $                                                         #GEOMETRYTYPE                   @  @                               '                                          @               @               'Ğ                   #NPOINT    #NDOF    #POINT    #INTEGRATOR    #MATERIAL    #GEOMETRY    #GETNPOINT    #GETNDOF    #GETPOINTID    #GETINTEGRATOR    #GETMATERIAL    #GETGEOMETRY    #SETNPOINT    #SETNDOF ¡   #SETPOINT ¥   #SETINTEGRATOR ª   #GETONEPOINT ®   #GETALLPOINTS ²                $                                                              $                                                           $                                                              #POINTPTRTYPE              &                                                         $                                          P              #INTEGRATORPTRTYPE R                 $                                          Ğ              #MATERIALPTRTYPE q                 $                                          P             #GEOMETRYPTRTYPE    1         À    $                                             #GETNPOINT    %         @   @                                                       #THIS              
                                     Ğ              #ELEMENTTYPE    1         À    $                                             #GETNDOF    %         @   @                                                       #THIS              
                                     Ğ              #ELEMENTTYPE    1         À    $                                        	     #GETPOINTID    %         @   @                                                       #THIS    #I              
D @                                   Ğ              #ELEMENTTYPE              
                                            1         À    $                                        
     #GETINTEGRATOR    &         @   @                                                       #THIS    #INTEGRATORPTRTYPE R             
                                     Ğ              #ELEMENTTYPE    1         À    $                                             #GETMATERIAL    &         @   @                                                       #THIS    #MATERIALPTRTYPE q             
                                     Ğ              #ELEMENTTYPE    1         À    $                                             #GETGEOMETRY    &         @   @                                                       #THIS    #GEOMETRYPTRTYPE              
                                     Ğ              #ELEMENTTYPE    1         À    $                                              #SETNPOINT    #         @     @                                                 #THIS    #NPOINT               
D                                     Ğ              #ELEMENTTYPE              
                                             1         À    $                            ¡                  #SETNDOF ¢   #         @     @                             ¢                    #THIS £   #NDOF ¤             
D                                £     Ğ              #ELEMENTTYPE              
                                 ¤           1         À    $                            ¥              	    #SETPOINT ¦   #         @     @                             ¦                    #THIS §   #I ¨   #POINT ©             
D                                §     Ğ              #ELEMENTTYPE              
                                 ¨                     
                                 ©                    #POINTTYPE 
   1         À    $                            ª              
    #SETINTEGRATOR «   #         @     @                             «                    #THIS ¬   #INTEGRATOR ­             
D                                ¬     Ğ              #ELEMENTTYPE              
                                  ­                  #INTEGRATORTYPE T   1         À    $                           ®                  #GETONEPOINT ¯   &         @   @                            ¯                           #THIS °   #I ±   #POINTPTRTYPE              
                                °     Ğ              #ELEMENTTYPE              
                                 ±           1         À    $                           ²                  #GETALLPOINTS ³   )        `   @                             ³                                         #THIS ´   #POINTPTRTYPE    p          5 8 O#ELEMENTTYPE     p        U            5 8 O#ELEMENTTYPE     p        U                                   
                                ´     Ğ              #ELEMENTTYPE           f      fn#fn          b   uapp(ELEMENTMOD    "  @   J  TOOLS    b  @   J  POINTMOD    ¢  @   J  POINTPTRMOD    â  @   J  INTEGRATORMOD !   "  @   J  INTEGRATORPTRMOD    b  @   J  MATERIALPTRMOD    ¢  @   J  GEOMETRYPTRMOD )   â  ¹      POINTPTRTYPE+POINTPTRMOD -     _   a   POINTPTRTYPE%PTR+POINTPTRMOD #   ú  É       POINTTYPE+POINTMOD )   Ã  H   %   POINTTYPE%ID+POINTMOD=ID '     H   %   POINTTYPE%X+POINTMOD=X '   S  H   %   POINTTYPE%Y+POINTMOD=Y '     H   %   POINTTYPE%Z+POINTMOD=Z (   ã  R   a   POINTTYPE%INIT+POINTMOD    5  o      INIT+POINTMOD #   ¤  W   a   INIT%THIS+POINTMOD !   û  @   a   INIT%ID+POINTMOD     ;  @   a   INIT%X+POINTMOD     {  @   a   INIT%Y+POINTMOD     »  @   a   INIT%Z+POINTMOD )   û  S   a   POINTTYPE%SETID+POINTMOD    N  Z      SETID+POINTMOD $   ¨  W   a   SETID%THIS+POINTMOD "   ÿ  @   a   SETID%ID+POINTMOD )   ?	  S   a   POINTTYPE%GETID+POINTMOD    	  Z      GETID+POINTMOD $   ì	  W   a   GETID%THIS+POINTMOD (   C
  R   a   POINTTYPE%SETX+POINTMOD    
  Y      SETX+POINTMOD #   î
  W   a   SETX%THIS+POINTMOD     E  @   a   SETX%X+POINTMOD (     R   a   POINTTYPE%GETX+POINTMOD    ×  Z      GETX+POINTMOD #   1  W   a   GETX%THIS+POINTMOD (     R   a   POINTTYPE%SETY+POINTMOD    Ú  Y      SETY+POINTMOD #   3  W   a   SETY%THIS+POINTMOD       @   a   SETY%Y+POINTMOD (   Ê  R   a   POINTTYPE%GETY+POINTMOD      Z      GETY+POINTMOD #   v  W   a   GETY%THIS+POINTMOD (   Í  R   a   POINTTYPE%SETZ+POINTMOD      Y      SETZ+POINTMOD #   x  W   a   SETZ%THIS+POINTMOD     Ï  @   a   SETZ%Z+POINTMOD (     R   a   POINTTYPE%GETZ+POINTMOD    a  Z      GETZ+POINTMOD #   »  W   a   GETZ%THIS+POINTMOD 2     V   a   POINTPTRTYPE%ALLOCATE+POINTPTRMOD %   h  ]      ALLOCATE+POINTPTRMOD *   Å  Z   a   ALLOCATE%THIS+POINTPTRMOD +     W   a   ALLOCATE%POINT+POINTPTRMOD /   v  S   a   POINTPTRTYPE%SETID+POINTPTRMOD "   É  Z      SETID+POINTPTRMOD '   #  Z   a   SETID%THIS+POINTPTRMOD %   }  @   a   SETID%ID+POINTPTRMOD /   ½  S   a   POINTPTRTYPE%GETID+POINTPTRMOD "     Z      GETID+POINTPTRMOD '   j  Z   a   GETID%THIS+POINTPTRMOD .   Ä  R   a   POINTPTRTYPE%SETX+POINTPTRMOD !     Y      SETX+POINTPTRMOD &   o  Z   a   SETX%THIS+POINTPTRMOD #   É  @   a   SETX%X+POINTPTRMOD .   	  R   a   POINTPTRTYPE%GETX+POINTPTRMOD !   [  Z      GETX+POINTPTRMOD &   µ  Z   a   GETX%THIS+POINTPTRMOD .     R   a   POINTPTRTYPE%SETY+POINTPTRMOD !   a  Y      SETY+POINTPTRMOD &   º  Z   a   SETY%THIS+POINTPTRMOD #     @   a   SETY%Y+POINTPTRMOD .   T  R   a   POINTPTRTYPE%GETY+POINTPTRMOD !   ¦  Z      GETY+POINTPTRMOD &      Z   a   GETY%THIS+POINTPTRMOD .   Z  R   a   POINTPTRTYPE%SETZ+POINTPTRMOD !   ¬  Y      SETZ+POINTPTRMOD &     Z   a   SETZ%THIS+POINTPTRMOD #   _  @   a   SETZ%Z+POINTPTRMOD .     R   a   POINTPTRTYPE%GETZ+POINTPTRMOD !   ñ  Z      GETZ+POINTPTRMOD &   K  Z   a   GETZ%THIS+POINTPTRMOD 3   ¥  g      INTEGRATORPTRTYPE+INTEGRATORPTRMOD 7     d   a   INTEGRATORPTRTYPE%PTR+INTEGRATORPTRMOD -   p  ñ       INTEGRATORTYPE+INTEGRATORMOD 8   a  H   a   INTEGRATORTYPE%GAUSSORDER+INTEGRATORMOD 8   ©  H   a   INTEGRATORTYPE%INTEGTERMS+INTEGRATORMOD 4   ñ     a   INTEGRATORTYPE%WEIGHT+INTEGRATORMOD 4     ¬   a   INTEGRATORTYPE%GPOINT+INTEGRATORMOD 7   1  ¬   a   INTEGRATORTYPE%SHAPEFUNC+INTEGRATORMOD 8   İ  Ä   a   INTEGRATORTYPE%DSHAPEFUNC+INTEGRATORMOD 2   ¡   R   a   INTEGRATORTYPE%INIT+INTEGRATORMOD #   ó   l      INIT+INTEGRATORMOD (   _!  \   a   INIT%THIS+INTEGRATORMOD .   »!  @   a   INIT%GAUSSORDER+INTEGRATORMOD (   û!  L   a   INIT%TYPE+INTEGRATORMOD :   G"  Z   a   INTEGRATORTYPE%VALUEGPOINTS+INTEGRATORMOD +   ¡"  \      VALUEGPOINTS+INTEGRATORMOD 0   ı"  \   a   VALUEGPOINTS%THIS+INTEGRATORMOD 0   Y#  L   a   VALUEGPOINTS%TYPE+INTEGRATORMOD ;   ¥#  T   %   INTEGRATORTYPE%GETG1D+INTEGRATORMOD=GETG1D %   ù#  R      GETG1D+INTEGRATORMOD *   K$  \   a   GETG1D%THIS+INTEGRATORMOD G   §$  Z   %   INTEGRATORTYPE%GETGTRIANGLE+INTEGRATORMOD=GETGTRIANGLE +   %  R      GETGTRIANGLE+INTEGRATORMOD 0   S%  \   a   GETGTRIANGLE%THIS+INTEGRATORMOD C   ¯%  X   %   INTEGRATORTYPE%GETGSQUARE+INTEGRATORMOD=GETGSQUARE )   &  R      GETGSQUARE+INTEGRATORMOD .   Y&  \   a   GETGSQUARE%THIS+INTEGRATORMOD <   µ&  V   a   INTEGRATORPTRTYPE%ALLOCATE+INTEGRATORPTRMOD *   '  b      ALLOCATE+INTEGRATORPTRMOD /   m'  _   a   ALLOCATE%THIS+INTEGRATORPTRMOD 5   Ì'  \   a   ALLOCATE%INTEGRATOR+INTEGRATORPTRMOD /   ((  g      MATERIALPTRTYPE+MATERIALPTRMOD 3   (  d   a   MATERIALPTRTYPE%PTR+MATERIALPTRMOD -   ó(  ~      THMATERIALTYPE+THMATERIALMOD :   q)  b   a   THMATERIALTYPE%MATERIALTYPE+THMATERIALMOD )   Ó)  P      MATERIALTYPE+MATERIALMOD :   #*     a   THMATERIALTYPE%CONDUCTIVITY+THMATERIALMOD 2   ¿*  R   a   THMATERIALTYPE%INIT+THMATERIALMOD #   +  b      INIT+THMATERIALMOD (   s+  \   a   INIT%THIS+THMATERIALMOD &   Ï+  @   a   INIT%KX+THMATERIALMOD &   ,  @   a   INIT%KY+THMATERIALMOD 8   O,  V   a   MATERIALPTRTYPE%ALLOCATE+MATERIALPTRMOD (   ¥,  `      ALLOCATE+MATERIALPTRMOD -   -  ]   a   ALLOCATE%THIS+MATERIALPTRMOD 1   b-  \   a   ALLOCATE%MATERIAL+MATERIALPTRMOD /   ¾-  Y      GEOMETRYPTRTYPE+GEOMETRYPTRMOD 3   .  b   a   GEOMETRYPTRTYPE%PTR+GEOMETRYPTRMOD )   y.  P      GEOMETRYTYPE+GEOMETRYMOD    É.  ^      ELEMENTTYPE #   '0  H   a   ELEMENTTYPE%NPOINT !   o0  H   a   ELEMENTTYPE%NDOF "   ·0  ¦   a   ELEMENTTYPE%POINT '   ]1  g   a   ELEMENTTYPE%INTEGRATOR %   Ä1  e   a   ELEMENTTYPE%MATERIAL %   )2  e   a   ELEMENTTYPE%GEOMETRY &   2  W   a   ELEMENTTYPE%GETNPOINT    å2  Z      GETNPOINT    ?3  Y   a   GETNPOINT%THIS $   3  U   a   ELEMENTTYPE%GETNDOF    í3  Z      GETNDOF    G4  Y   a   GETNDOF%THIS '    4  X   a   ELEMENTTYPE%GETPOINTID    ø4  a      GETPOINTID     Y5  Y   a   GETPOINTID%THIS    ²5  @   a   GETPOINTID%I *   ò5  [   a   ELEMENTTYPE%GETINTEGRATOR    M6  q      GETINTEGRATOR #   ¾6  Y   a   GETINTEGRATOR%THIS (   7  Y   a   ELEMENTTYPE%GETMATERIAL    p7  o      GETMATERIAL !   ß7  Y   a   GETMATERIAL%THIS (   88  Y   a   ELEMENTTYPE%GETGEOMETRY    8  o      GETGEOMETRY !    9  Y   a   GETGEOMETRY%THIS &   Y9  W   a   ELEMENTTYPE%SETNPOINT    °9  ^      SETNPOINT    :  Y   a   SETNPOINT%THIS !   g:  @   a   SETNPOINT%NPOINT $   §:  U   a   ELEMENTTYPE%SETNDOF    ü:  \      SETNDOF    X;  Y   a   SETNDOF%THIS    ±;  @   a   SETNDOF%NDOF %   ñ;  V   a   ELEMENTTYPE%SETPOINT    G<  d      SETPOINT    «<  Y   a   SETPOINT%THIS    =  @   a   SETPOINT%I    D=  W   a   SETPOINT%POINT *   =  [   a   ELEMENTTYPE%SETINTEGRATOR    ö=  b      SETINTEGRATOR #   X>  Y   a   SETINTEGRATOR%THIS )   ±>  \   a   SETINTEGRATOR%INTEGRATOR (   ?  Y   a   ELEMENTTYPE%GETONEPOINT    f?  s      GETONEPOINT !   Ù?  Y   a   GETONEPOINT%THIS    2@  @   a   GETONEPOINT%I )   r@  Z   a   ELEMENTTYPE%GETALLPOINTS    Ì@       GETALLPOINTS "   ŞA  Y   a   GETALLPOINTS%THIS 