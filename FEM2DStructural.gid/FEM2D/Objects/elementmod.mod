  ¢E  Â   k820309    9          19.0        {¯Ş]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DStructural.gid/FEM2D/Source/Element/Element.f90 ELEMENTMOD              ELEMENTTYPE                                                     
                            @                              
                            @                              
                            @                              
                            @                              
                            @                              
                      @  @                               '              
      #PTR    #ALLOCATE 1   #SETID 5   #GETID 9   #SETX <   #GETX @   #SETY C   #GETY G   #SETZ J   #GETZ N                $                                                          #POINTTYPE 	                     @                           	     '                     #ID 
   #X    #Y    #Z    #INIT    #SETID    #GETID    #SETX    #GETX     #SETY #   #GETY '   #SETZ *   #GETZ .                 D                             
                                 D                                            
                 D                                            
                 D                                            
   1         À    $                                              #INIT    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                                     #POINTTYPE 	             
                                                      
                                      
                
                                      
                
                                      
      1         À    $                                              #SETID    #         @     @                                                #THIS    #ID              
                                                     #POINTTYPE 	             
                                            1         À    $                                             #GETID    %         @   @                                                      #THIS              
                                                     #POINTTYPE 	   1         À    $                                              #SETX    #         @     @                                                #THIS    #X              
                                                     #POINTTYPE 	             
                                      
      1         À    $                                         	     #GETX !   %         @   @                           !                    
       #THIS "             
                                "                     #POINTTYPE 	   1         À    $                            #             
     #SETY $   #         @     @                            $                    #THIS %   #Y &             
                                %                     #POINTTYPE 	             
                                 &     
      1         À    $                           '                  #GETY (   %         @   @                           (                    
       #THIS )             
                                )                     #POINTTYPE 	   1         À    $                            *                  #SETZ +   #         @     @                            +                    #THIS ,   #Z -             
                                ,                     #POINTTYPE 	             
                                 -     
      1         À    $                           .              	    #GETZ /   %         @   @                           /                    
       #THIS 0             
                                0                     #POINTTYPE 	   1         À    $                            1                  #ALLOCATE 2   #         @     @                            2                    #THIS 3   #POINT 4             
                                3                    #POINTPTRTYPE              
                                  4                    #POINTTYPE 	   1         À    $                            5                  #SETID 6   #         @     @                            6                    #THIS 7   #ID 8             
                                7                    #POINTPTRTYPE              
                                 8           1         À    $                          9                  #GETID :   %         @   @                          :                           #THIS ;             
                                ;                    #POINTPTRTYPE    1         À    $                            <                  #SETX =   #         @     @                            =                    #THIS >   #X ?             
                                >                    #POINTPTRTYPE              
                                 ?     
      1         À    $                           @                  #GETX A   %         @   @                           A                    
       #THIS B             
                                B                    #POINTPTRTYPE    1         À    $                            C                  #SETY D   #         @     @                            D                    #THIS E   #Y F             
                                E                    #POINTPTRTYPE              
                                 F     
      1         À    $                           G                  #GETY H   %         @   @                           H                    
       #THIS I             
                                I                    #POINTPTRTYPE    1         À    $                            J             	     #SETZ K   #         @     @                            K                    #THIS L   #Z M             
                                L                    #POINTPTRTYPE              
                                 M     
      1         À    $                           N             
 	    #GETZ O   %         @   @                           O                    
       #THIS P             
                                P                    #POINTPTRTYPE                  @  @                          Q     '                    #PTR R   #ALLOCATE l                $                             R                          #INTEGRATORTYPE S                     @               @           S     '                   #GAUSSORDER T   #INTEGTERMS U   #WEIGHT V   #GPOINT W   #SHAPEFUNC X   #DSHAPEFUNC Y   #INIT Z   #VALUEGPOINTS _   #GETG1D c   #GETGTRIANGLE f   #GETGSQUARE i                 $                             T                                 $                             U                             $                             V                             
            &                                                      $                             W            P                 
            &                   &                                                       $                             X            °                 
            &                   &                                                       $                             Y                            
            &                   &                   &                                           1         À    $                            Z                  #INIT [   #         @     @                            [                    #THIS \   #GAUSSORDER ]   #TYPE ^             
                                \                   #INTEGRATORTYPE S             
                                 ]                     
                                ^                    1 1         À    $                            _                  #VALUEGPOINTS `   #         @     @                            `                    #THIS a   #TYPE b             
                                a                   #INTEGRATORTYPE S             
                                b                    1 1         À    D                            c             	     #GETG1D d   #         @     @                            d                    #THIS e             
                                e                   #INTEGRATORTYPE S   1         À    D                            f             
     #GETGTRIANGLE g   #         @     @                            g                    #THIS h             
                                h                   #INTEGRATORTYPE S   1         À    D                            i                  #GETGSQUARE j   #         @     @                            j                    #THIS k             
                                k                   #INTEGRATORTYPE S   1         À    $                            l                  #ALLOCATE m   #         @     @                            m                    #THIS n   #INTEGRATOR o             
                                n                    #INTEGRATORPTRTYPE Q             
                                  o                  #INTEGRATORTYPE S                 @  @                          p     '                    #PTR q   #ALLOCATE                 $                             q     P                      #STRUCTMATERIALTYPE r                 @  @                         r     'P                    #MATERIALTYPE s   #YOUNG u   #POISSONCOEF v   #THERMALCOEF w   #AREA x   #THICKNESS y   #D11 z   #D12 {   #D21 |   #D22 }   #D33 ~   #INIT                  $                              s                            #MATERIALTYPE t                 @  @                          t     '                                      $                             u                
                 $                             v               
                 $                             w               
                 $                             x               
                 $                             y                
                 $                             z     (          
                 $                             {     0          
                 $                             |     8       	   
                 $                             }     @       
   
                 $                             ~     H          
   1         À    $                                              #INIT    #         @     @                                                #THIS    #YOUNG    #POISSONCOEF    #THERMALCOEF    #AREA    #THICKNESS              
                                     P               #STRUCTMATERIALTYPE r             
                                      
                
                                      
                
                                      
                
                                      
                
                                      
      1         À    $                                              #ALLOCATE    #         @     @                                                #THIS    #MATERIAL                                                                  #MATERIALPTRTYPE p             
                                       P              #STRUCTMATERIALTYPE r                     @               @               'X                   #ID    #NPOINT    #NDOF    #POINT    #INTEGRATOR    #MATERIAL    #GETID    #GETNPOINT    #GETNDOF    #GETPOINTID    #GETINTEGRATOR    #GETMATERIAL ¢   #SETID ¥   #SETNPOINT ©   #SETNDOF ­   #SETPOINT ±   #SETINTEGRATOR ¶   #GETONEPOINT º   #GETALLPOINTS ¾                 $                                                             $                                                             $                                                           $                                                              #POINTPTRTYPE              &                                                         $                                          X              #INTEGRATORPTRTYPE Q                 $                                          Ø              #MATERIALPTRTYPE p   1         À    $                                             #GETID    %         @   @                                                       #THIS              
                                     X              #ELEMENTTYPE    1         À    $                                             #GETNPOINT    %         @   @                                                       #THIS              
                                     X              #ELEMENTTYPE    1         À    $                                        	     #GETNDOF    %         @   @                                                       #THIS              
                                     X              #ELEMENTTYPE    1         À    $                                        
     #GETPOINTID    %         @   @                                                       #THIS    #I              
D @                                   X              #ELEMENTTYPE              
                                            1         À    $                                             #GETINTEGRATOR     &         @   @                                                        #THIS ¡   #INTEGRATORPTRTYPE Q             
                                ¡     X              #ELEMENTTYPE    1         À    $                           ¢                  #GETMATERIAL £   &         @   @                            £                           #THIS ¤   #MATERIALPTRTYPE p             
                                ¤     X              #ELEMENTTYPE    1         À    $                            ¥                  #SETID ¦   #         @     @                             ¦                    #THIS §   #ID ¨             
D                                §     X              #ELEMENTTYPE              
                                 ¨           1         À    $                            ©                  #SETNPOINT ª   #         @     @                             ª                    #THIS «   #NPOINT ¬             
D                                «     X              #ELEMENTTYPE              
                                 ¬           1         À    $                            ­              	    #SETNDOF ®   #         @     @                             ®                    #THIS ¯   #NDOF °             
D                                ¯     X              #ELEMENTTYPE              
                                 °           1         À    $                            ±              
    #SETPOINT ²   #         @     @                             ²                    #THIS ³   #I ´   #POINT µ             
D                                ³     X              #ELEMENTTYPE              
                                 ´                     
                                 µ                    #POINTTYPE 	   1         À    $                            ¶                  #SETINTEGRATOR ·   #         @     @                             ·                    #THIS ¸   #INTEGRATOR ¹             
D                                ¸     X              #ELEMENTTYPE              
                                  ¹                  #INTEGRATORTYPE S   1         À    $                           º                  #GETONEPOINT »   &         @   @                            »                           #THIS ¼   #I ½   #POINTPTRTYPE              
                                ¼     X              #ELEMENTTYPE              
                                 ½           1         À    $                           ¾                  #GETALLPOINTS ¿   )        `   @                             ¿                                         #THIS À   #POINTPTRTYPE    p          5 8 O#ELEMENTTYPE     p        U            5 8 O#ELEMENTTYPE     p        U                                   
                                À     X              #ELEMENTTYPE           x      fn#fn          b   uapp(ELEMENTMOD    4  @   J  TOOLS    t  @   J  POINTMOD    ´  @   J  POINTPTRMOD    ô  @   J  INTEGRATORMOD !   4  @   J  INTEGRATORPTRMOD    t  @   J  MATERIALPTRMOD )   ´  ¹      POINTPTRTYPE+POINTPTRMOD -   m  _   a   POINTPTRTYPE%PTR+POINTPTRMOD #   Ì  É       POINTTYPE+POINTMOD )     H   %   POINTTYPE%ID+POINTMOD=ID '   İ  H   %   POINTTYPE%X+POINTMOD=X '   %  H   %   POINTTYPE%Y+POINTMOD=Y '   m  H   %   POINTTYPE%Z+POINTMOD=Z (   µ  R   a   POINTTYPE%INIT+POINTMOD      o      INIT+POINTMOD #   v  W   a   INIT%THIS+POINTMOD !   Í  @   a   INIT%ID+POINTMOD       @   a   INIT%X+POINTMOD     M  @   a   INIT%Y+POINTMOD       @   a   INIT%Z+POINTMOD )   Í  S   a   POINTTYPE%SETID+POINTMOD       Z      SETID+POINTMOD $   z  W   a   SETID%THIS+POINTMOD "   Ñ  @   a   SETID%ID+POINTMOD )   	  S   a   POINTTYPE%GETID+POINTMOD    d	  Z      GETID+POINTMOD $   ¾	  W   a   GETID%THIS+POINTMOD (   
  R   a   POINTTYPE%SETX+POINTMOD    g
  Y      SETX+POINTMOD #   À
  W   a   SETX%THIS+POINTMOD       @   a   SETX%X+POINTMOD (   W  R   a   POINTTYPE%GETX+POINTMOD    ©  Z      GETX+POINTMOD #     W   a   GETX%THIS+POINTMOD (   Z  R   a   POINTTYPE%SETY+POINTMOD    ¬  Y      SETY+POINTMOD #     W   a   SETY%THIS+POINTMOD     \  @   a   SETY%Y+POINTMOD (     R   a   POINTTYPE%GETY+POINTMOD    î  Z      GETY+POINTMOD #   H  W   a   GETY%THIS+POINTMOD (     R   a   POINTTYPE%SETZ+POINTMOD    ñ  Y      SETZ+POINTMOD #   J  W   a   SETZ%THIS+POINTMOD     ¡  @   a   SETZ%Z+POINTMOD (   á  R   a   POINTTYPE%GETZ+POINTMOD    3  Z      GETZ+POINTMOD #     W   a   GETZ%THIS+POINTMOD 2   ä  V   a   POINTPTRTYPE%ALLOCATE+POINTPTRMOD %   :  ]      ALLOCATE+POINTPTRMOD *     Z   a   ALLOCATE%THIS+POINTPTRMOD +   ñ  W   a   ALLOCATE%POINT+POINTPTRMOD /   H  S   a   POINTPTRTYPE%SETID+POINTPTRMOD "     Z      SETID+POINTPTRMOD '   õ  Z   a   SETID%THIS+POINTPTRMOD %   O  @   a   SETID%ID+POINTPTRMOD /     S   a   POINTPTRTYPE%GETID+POINTPTRMOD "   â  Z      GETID+POINTPTRMOD '   <  Z   a   GETID%THIS+POINTPTRMOD .     R   a   POINTPTRTYPE%SETX+POINTPTRMOD !   è  Y      SETX+POINTPTRMOD &   A  Z   a   SETX%THIS+POINTPTRMOD #     @   a   SETX%X+POINTPTRMOD .   Û  R   a   POINTPTRTYPE%GETX+POINTPTRMOD !   -  Z      GETX+POINTPTRMOD &     Z   a   GETX%THIS+POINTPTRMOD .   á  R   a   POINTPTRTYPE%SETY+POINTPTRMOD !   3  Y      SETY+POINTPTRMOD &     Z   a   SETY%THIS+POINTPTRMOD #   æ  @   a   SETY%Y+POINTPTRMOD .   &  R   a   POINTPTRTYPE%GETY+POINTPTRMOD !   x  Z      GETY+POINTPTRMOD &   Ò  Z   a   GETY%THIS+POINTPTRMOD .   ,  R   a   POINTPTRTYPE%SETZ+POINTPTRMOD !   ~  Y      SETZ+POINTPTRMOD &   ×  Z   a   SETZ%THIS+POINTPTRMOD #   1  @   a   SETZ%Z+POINTPTRMOD .   q  R   a   POINTPTRTYPE%GETZ+POINTPTRMOD !   Ã  Z      GETZ+POINTPTRMOD &     Z   a   GETZ%THIS+POINTPTRMOD 3   w  g      INTEGRATORPTRTYPE+INTEGRATORPTRMOD 7   Ş  d   a   INTEGRATORPTRTYPE%PTR+INTEGRATORPTRMOD -   B  ñ       INTEGRATORTYPE+INTEGRATORMOD 8   3  H   a   INTEGRATORTYPE%GAUSSORDER+INTEGRATORMOD 8   {  H   a   INTEGRATORTYPE%INTEGTERMS+INTEGRATORMOD 4   Ã     a   INTEGRATORTYPE%WEIGHT+INTEGRATORMOD 4   W  ¬   a   INTEGRATORTYPE%GPOINT+INTEGRATORMOD 7     ¬   a   INTEGRATORTYPE%SHAPEFUNC+INTEGRATORMOD 8   ¯  Ä   a   INTEGRATORTYPE%DSHAPEFUNC+INTEGRATORMOD 2   s   R   a   INTEGRATORTYPE%INIT+INTEGRATORMOD #   Å   l      INIT+INTEGRATORMOD (   1!  \   a   INIT%THIS+INTEGRATORMOD .   !  @   a   INIT%GAUSSORDER+INTEGRATORMOD (   Í!  L   a   INIT%TYPE+INTEGRATORMOD :   "  Z   a   INTEGRATORTYPE%VALUEGPOINTS+INTEGRATORMOD +   s"  \      VALUEGPOINTS+INTEGRATORMOD 0   Ï"  \   a   VALUEGPOINTS%THIS+INTEGRATORMOD 0   +#  L   a   VALUEGPOINTS%TYPE+INTEGRATORMOD ;   w#  T   %   INTEGRATORTYPE%GETG1D+INTEGRATORMOD=GETG1D %   Ë#  R      GETG1D+INTEGRATORMOD *   $  \   a   GETG1D%THIS+INTEGRATORMOD G   y$  Z   %   INTEGRATORTYPE%GETGTRIANGLE+INTEGRATORMOD=GETGTRIANGLE +   Ó$  R      GETGTRIANGLE+INTEGRATORMOD 0   %%  \   a   GETGTRIANGLE%THIS+INTEGRATORMOD C   %  X   %   INTEGRATORTYPE%GETGSQUARE+INTEGRATORMOD=GETGSQUARE )   Ù%  R      GETGSQUARE+INTEGRATORMOD .   +&  \   a   GETGSQUARE%THIS+INTEGRATORMOD <   &  V   a   INTEGRATORPTRTYPE%ALLOCATE+INTEGRATORPTRMOD *   İ&  b      ALLOCATE+INTEGRATORPTRMOD /   ?'  _   a   ALLOCATE%THIS+INTEGRATORPTRMOD 5   '  \   a   ALLOCATE%INTEGRATOR+INTEGRATORPTRMOD /   ú'  g      MATERIALPTRTYPE+MATERIALPTRMOD 3   a(  h   a   MATERIALPTRTYPE%PTR+MATERIALPTRMOD 5   É(  ß      STRUCTMATERIALTYPE+STRUCTMATERIALMOD B   ¨)  b   a   STRUCTMATERIALTYPE%MATERIALTYPE+STRUCTMATERIALMOD )   
*  P      MATERIALTYPE+MATERIALMOD ;   Z*  H   a   STRUCTMATERIALTYPE%YOUNG+STRUCTMATERIALMOD A   ¢*  H   a   STRUCTMATERIALTYPE%POISSONCOEF+STRUCTMATERIALMOD A   ê*  H   a   STRUCTMATERIALTYPE%THERMALCOEF+STRUCTMATERIALMOD :   2+  H   a   STRUCTMATERIALTYPE%AREA+STRUCTMATERIALMOD ?   z+  H   a   STRUCTMATERIALTYPE%THICKNESS+STRUCTMATERIALMOD 9   Â+  H   a   STRUCTMATERIALTYPE%D11+STRUCTMATERIALMOD 9   
,  H   a   STRUCTMATERIALTYPE%D12+STRUCTMATERIALMOD 9   R,  H   a   STRUCTMATERIALTYPE%D21+STRUCTMATERIALMOD 9   ,  H   a   STRUCTMATERIALTYPE%D22+STRUCTMATERIALMOD 9   â,  H   a   STRUCTMATERIALTYPE%D33+STRUCTMATERIALMOD :   *-  R   a   STRUCTMATERIALTYPE%INIT+STRUCTMATERIALMOD '   |-        INIT+STRUCTMATERIALMOD ,   .  `   a   INIT%THIS+STRUCTMATERIALMOD -   t.  @   a   INIT%YOUNG+STRUCTMATERIALMOD 3   ´.  @   a   INIT%POISSONCOEF+STRUCTMATERIALMOD 3   ô.  @   a   INIT%THERMALCOEF+STRUCTMATERIALMOD ,   4/  @   a   INIT%AREA+STRUCTMATERIALMOD 1   t/  @   a   INIT%THICKNESS+STRUCTMATERIALMOD 8   ´/  V   a   MATERIALPTRTYPE%ALLOCATE+MATERIALPTRMOD (   
0  `      ALLOCATE+MATERIALPTRMOD -   j0  ]   a   ALLOCATE%THIS+MATERIALPTRMOD 1   Ç0  `   a   ALLOCATE%MATERIAL+MATERIALPTRMOD    '1  ]      ELEMENTTYPE    2  H   a   ELEMENTTYPE%ID #   Ì2  H   a   ELEMENTTYPE%NPOINT !   3  H   a   ELEMENTTYPE%NDOF "   \3  ¦   a   ELEMENTTYPE%POINT '   4  g   a   ELEMENTTYPE%INTEGRATOR %   i4  e   a   ELEMENTTYPE%MATERIAL "   Î4  S   a   ELEMENTTYPE%GETID    !5  Z      GETID    {5  Y   a   GETID%THIS &   Ô5  W   a   ELEMENTTYPE%GETNPOINT    +6  Z      GETNPOINT    6  Y   a   GETNPOINT%THIS $   Ş6  U   a   ELEMENTTYPE%GETNDOF    37  Z      GETNDOF    7  Y   a   GETNDOF%THIS '   æ7  X   a   ELEMENTTYPE%GETPOINTID    >8  a      GETPOINTID     8  Y   a   GETPOINTID%THIS    ø8  @   a   GETPOINTID%I *   89  [   a   ELEMENTTYPE%GETINTEGRATOR    9  q      GETINTEGRATOR #   :  Y   a   GETINTEGRATOR%THIS (   ]:  Y   a   ELEMENTTYPE%GETMATERIAL    ¶:  o      GETMATERIAL !   %;  Y   a   GETMATERIAL%THIS "   ~;  S   a   ELEMENTTYPE%SETID    Ñ;  Z      SETID    +<  Y   a   SETID%THIS    <  @   a   SETID%ID &   Ä<  W   a   ELEMENTTYPE%SETNPOINT    =  ^      SETNPOINT    y=  Y   a   SETNPOINT%THIS !   Ò=  @   a   SETNPOINT%NPOINT $   >  U   a   ELEMENTTYPE%SETNDOF    g>  \      SETNDOF    Ã>  Y   a   SETNDOF%THIS    ?  @   a   SETNDOF%NDOF %   \?  V   a   ELEMENTTYPE%SETPOINT    ²?  d      SETPOINT    @  Y   a   SETPOINT%THIS    o@  @   a   SETPOINT%I    ¯@  W   a   SETPOINT%POINT *   A  [   a   ELEMENTTYPE%SETINTEGRATOR    aA  b      SETINTEGRATOR #   ÃA  Y   a   SETINTEGRATOR%THIS )   B  \   a   SETINTEGRATOR%INTEGRATOR (   xB  Y   a   ELEMENTTYPE%GETONEPOINT    ÑB  s      GETONEPOINT !   DC  Y   a   GETONEPOINT%THIS    C  @   a   GETONEPOINT%I )   İC  Z   a   ELEMENTTYPE%GETALLPOINTS    7D       GETALLPOINTS "   IE  Y   a   GETALLPOINTS%THIS 