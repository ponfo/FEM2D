  I  Ä   k820309    9          19.0        }$Ü]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DThermal.gid/FEM2D/Source/BoundaryCondition/ThBoundaryCondition1D.f90 THBOUNDARYCONDITION1DMOD              THBOUNDARYCONDITION1DTYPE gen@THERMALBOUNDARYCONDITION1D                                                     
                       @                                  
                            @                              
                            @                              
                            @                              
                @          @                              
                                                              u #CONSTRUCTOR    &         @   @                                Ø                      #NDIRICHLET    #NNORMALFLUX 	   #NCONVECTION 
   #THBOUNDARYCONDITION1DTYPE              
  @                                                   
  @                              	                     
  @                              
                             @                               '                    #ID    #VALUE    #INIT    #APPLY                  $                                                              $                                            
   1         À    $                                              #INIT    #         @     @                                                #THIS    #ID    #VALUE              
                                                    #DIRICHLETPOINTTYPE              
                                                      
                                      
      1         À    $                                             #APPLY    #         @     @                                               #THIS    #STIFFNESS    #RHS              
                                                    #DIRICHLETPOINTTYPE              
                                                   #SPARSE              
                                                   
               &                                                             @                               '                    #POINTID    #VALUE    #INIT    #APPLY "                 $                                                              $                                            
   1         À    $                                              #INIT    #         @     @                                                #THIS    #POINTID     #VALUE !             
                                                    #NORMALFLUXPOINTTYPE              
                                                       
                                 !     
      1         À    $                           "                  #APPLY #   #         @     @                           #                    #THIS $   #RHS %             
                                $                    #NORMALFLUXPOINTTYPE              
                                %                   
               &                                                             @                          &     '                    #POINTID '   #COEF (   #TEMP )   #INIT *   #APPLY 0                 $                             '                                 $                             (               
                 $                             )               
   1         À    $                            *                  #INIT +   #         @     @                            +                    #THIS ,   #POINTID -   #COEF .   #TEMP /             
                                ,                    #CONVECTIONPOINTTYPE &             
                                 -                     
                                 .     
                
                                 /     
      1         À    $                           0                  #APPLY 1   #         @     @                           1                    #THIS 2   #STIFFNESS 3   #RHS 4             
                                2                    #CONVECTIONPOINTTYPE &             
                                 3                   #SPARSE              
                                4                   
               &                                                          @  @               Ä                '                   #A 5   #AI 6   #AJ 7   #ROWCOUNTER 8   #N 9   #NNZ :   #COUNTER ;   #TRIPLET <   #INIT A   #UPDATE F   #APPEND K   #MAKECRS Q   #GET U   #GETNNZ Z   #GETN ]   #PRINTVALUE `   #PRINTNONZEROS f   #PRINTALL j   #DELETEROWANDCOL n   #FREE s   #HANDLEDUPLICATES v               $                             5                              
            &                                                      $                             6            H                             &                                                      $                             7                                         &                                                       $                             8            Ø                             &                                                         $                             9                                $                             :     $                          $                             ;     (                          $                              <     Ø       0             #TRIPLET =                 @  @              D           =     'Ø                    #A >   #ROW ?   #COL @                                            >                              
            &                                                                                    ?            H                             &                                                                                    @                                         &                                           1         À    $                            A             	     #INIT B   #         @     @                            B                    #THIS C   #NNZ D   #ROWS E             
                                C                   #SPARSE              
                                 D                     
                                 E           1         À    $                            F             
     #UPDATE G   #         @     @                            G                    #THIS H   #NNZ I   #ROWS J             
                                H                   #SPARSE              
                                 I                     
                                 J           1         À    $                            K                  #APPEND L   #         @     @                            L                    #THIS M   #VAL N   #ROW O   #COL P             
                                M                   #SPARSE              
                                 N     
                
                                 O                     
                                 P           1         À    $                            Q                  #MAKECRS R   #         @     @                            R                    #THIS S   #SORTROWS T             
                                S                   #SPARSE              
                                 T           1         À    $                           U                  #GET V   %         @   @                           V                    
       #THIS W   #I X   #J Y             
                                W                   #SPARSE              
                                 X                     
                                 Y           1         À    $                           Z                  #GETNNZ [   %         @   @                           [                           #THIS \             
                                \                   #SPARSE    1         À    $                           ]                  #GETN ^   %         @   @                           ^                           #THIS _             
                                _                   #SPARSE    1         À    $                            `                  #PRINTVALUE a   #         @     @                            a                    #THIS b   #I c   #J d   #FILENAME e             
                                b                   #SPARSE              
                                 c                     
                                 d                     
                               e                    1 1         À    $                            f              	    #PRINTNONZEROS g   #         @     @                            g                    #THIS h   #FILENAME i             
                                h                   #SPARSE              
                               i                    1 1         À    $                            j              
    #PRINTALL k   #         @     @                            k                    #THIS l   #FILENAME m             
                                l                   #SPARSE              
                               m                    1 1         À    $                            n                  #DELETEROWANDCOL o   #         @     @                            o                    #THIS p   #ROW q   #COL r             
                                p                   #SPARSE              
                                 q                     
                                 r           1         À    $                            s                  #FREE t   #         @     @                            t                    #THIS u             
                                u                   #SPARSE    1         À    D                           v                  #HANDLEDUPLICATES w   #         @     @                            w                    #THIS x             
                                x                   #SPARSE                      @               À           y     '                   #A z   #AI {   #AJ |   #ROWCOUNTER }   #N ~   #NNZ    #COUNTER    #TRIPLET    #INIT    #UPDATE    #APPEND    #MAKECRS    #GET    #GETNNZ    #GETN    #PRINTVALUE    #PRINTNONZEROS    #PRINTALL    #DELETEROWANDCOL    #FREE    #HANDLEDUPLICATES                $                             z                              
            &                                                      $                             {            H                             &                                                      $                             |                                         &                                                       $                             }            Ø                             &                                                         $                             ~                                $                                  $                          $                                  (                          $                                   Ø       0             #TRIPLET =   1         À    $                                         	     #INIT B   1         À    $                                         
     #UPDATE G   1         À    $                                              #APPEND L   1         À    $                                              #MAKECRS R   1         À    $                                             #GET V   1         À    $                                             #GETNNZ [   1         À    $                                             #GETN ^   1         À    $                                              #PRINTVALUE a   1         À    $                                          	    #PRINTNONZEROS g   1         À    $                                          
    #PRINTALL k   1         À    $                                              #DELETEROWANDCOL o   1         À    $                                              #FREE t   1         À    D                                              #HANDLEDUPLICATES w                     @               @                'Ø                    #DIRICHLETPOINT    #NORMALFLUXPOINT    #CONVECTIONPOINT    #INIT    #ADDDIRICHLETPOINT    #GETDIRICHLETPOINT    #GETNDIRICHLET ¡   #ADDNORMALFLUXPOINT ¤   #GETNORMALFLUXPOINT ©   #GETNNORMALFLUX ­   #ADDCONVECTIONPOINT °   #GETCONVECTIONPOINT ¶   #GETNCONVECTION º   #APPLY ½              $                                                               #DIRICHLETPOINTTYPE              &                                                      $                                          H                    #NORMALFLUXPOINTTYPE              &                                                      $                                                              #CONVECTIONPOINTTYPE &             &                                           1         À    $                                             #INIT    #         @     @                                                #THIS    #NDIRICHLET    #NNORMALFLUX    #NCONVECTION              
D                                     Ø               #THBOUNDARYCONDITION1DTYPE              
                                                      
                                                      
                                            1         À    $                                              #ADDDIRICHLETPOINT    #         @     @                                                 #THIS    #ID    #VALUE              
D                                     Ø               #THBOUNDARYCONDITION1DTYPE              
  @                                                   
  @                                   
      1         À    $                                             #GETDIRICHLETPOINT    &         @   @                                                       #THIS    #I     #DIRICHLETPOINTTYPE              
                                     Ø               #THBOUNDARYCONDITION1DTYPE              
                                             1         À    $                          ¡                  #GETNDIRICHLET ¢   %         @   @                           ¢                           #THIS £             
                                £     Ø               #THBOUNDARYCONDITION1DTYPE    1         À    $                            ¤                  #ADDNORMALFLUXPOINT ¥   #         @     @                             ¥                    #THIS ¦   #POINTID §   #VALUE ¨             
D                                ¦     Ø               #THBOUNDARYCONDITION1DTYPE              
  @                              §                     
  @                              ¨     
      1         À    $                           ©             	     #GETNORMALFLUXPOINT ª   &         @   @                            ª                           #THIS «   #I ¬   #NORMALFLUXPOINTTYPE              
                                «     Ø               #THBOUNDARYCONDITION1DTYPE              
                                 ¬           1         À    $                          ­             
     #GETNNORMALFLUX ®   %         @   @                           ®                           #THIS ¯             
                                ¯     Ø               #THBOUNDARYCONDITION1DTYPE    1         À    $                            °                  #ADDCONVECTIONPOINT ±   #         @     @                             ±                    #THIS ²   #POINTID ³   #COEF ´   #TEMP µ             
D                                ²     Ø               #THBOUNDARYCONDITION1DTYPE              
  @                              ³                     
  @                              ´     
                
  @                              µ     
      1         À    $                           ¶              	    #GETCONVECTIONPOINT ·   &         @   @                            ·                           #THIS ¸   #I ¹   #CONVECTIONPOINTTYPE &             
                                ¸     Ø               #THBOUNDARYCONDITION1DTYPE              
                                 ¹           1         À    $                          º              
    #GETNCONVECTION »   %         @   @                           »                           #THIS ¼             
                                ¼     Ø               #THBOUNDARYCONDITION1DTYPE    1         À    $                            ½                  #APPLY ¾   #         @     @                             ¾                    #THIS ¿   #STIFFNESS À   #RHS Á             
D @                              ¿     Ø               #THBOUNDARYCONDITION1DTYPE              
D @                              À                   #SPARSE y             
D @                              Á                   
               &                                                        fn#fn .   ;  I   b   uapp(THBOUNDARYCONDITION1DMOD      @   J  TOOLS    Ä  @   J  DEBUGGERMOD "     @   J  DIRICHLETPOINTMOD #   D  @   J  NORMALFLUXPOINTMOD #     @   J  CONVECTIONPOINTMOD    Ä  @   J  SPARSEKIT /     Q       gen@THERMALBOUNDARYCONDITION1D    U  ¡      CONSTRUCTOR '   ö  @   a   CONSTRUCTOR%NDIRICHLET (   6  @   a   CONSTRUCTOR%NNORMALFLUX (   v  @   a   CONSTRUCTOR%NCONVECTION 5   ¶  x       DIRICHLETPOINTTYPE+DIRICHLETPOINTMOD 8   .  H   a   DIRICHLETPOINTTYPE%ID+DIRICHLETPOINTMOD ;   v  H   a   DIRICHLETPOINTTYPE%VALUE+DIRICHLETPOINTMOD :   ¾  R   a   DIRICHLETPOINTTYPE%INIT+DIRICHLETPOINTMOD '     e      INIT+DIRICHLETPOINTMOD ,   u  `   a   INIT%THIS+DIRICHLETPOINTMOD *   Õ  @   a   INIT%ID+DIRICHLETPOINTMOD -     @   a   INIT%VALUE+DIRICHLETPOINTMOD ;   U  S   a   DIRICHLETPOINTTYPE%APPLY+DIRICHLETPOINTMOD (   ¨  j      APPLY+DIRICHLETPOINTMOD -     `   a   APPLY%THIS+DIRICHLETPOINTMOD 2   r  T   a   APPLY%STIFFNESS+DIRICHLETPOINTMOD ,   Æ     a   APPLY%RHS+DIRICHLETPOINTMOD 7   R	  }       NORMALFLUXPOINTTYPE+NORMALFLUXPOINTMOD ?   Ï	  H   a   NORMALFLUXPOINTTYPE%POINTID+NORMALFLUXPOINTMOD =   
  H   a   NORMALFLUXPOINTTYPE%VALUE+NORMALFLUXPOINTMOD <   _
  R   a   NORMALFLUXPOINTTYPE%INIT+NORMALFLUXPOINTMOD (   ±
  j      INIT+NORMALFLUXPOINTMOD -     a   a   INIT%THIS+NORMALFLUXPOINTMOD 0   |  @   a   INIT%POINTID+NORMALFLUXPOINTMOD .   ¼  @   a   INIT%VALUE+NORMALFLUXPOINTMOD =   ü  S   a   NORMALFLUXPOINTTYPE%APPLY+NORMALFLUXPOINTMOD )   O  [      APPLY+NORMALFLUXPOINTMOD .   ª  a   a   APPLY%THIS+NORMALFLUXPOINTMOD -        a   APPLY%RHS+NORMALFLUXPOINTMOD 7            CONVECTIONPOINTTYPE+CONVECTIONPOINTMOD ?     H   a   CONVECTIONPOINTTYPE%POINTID+CONVECTIONPOINTMOD <   e  H   a   CONVECTIONPOINTTYPE%COEF+CONVECTIONPOINTMOD <   ­  H   a   CONVECTIONPOINTTYPE%TEMP+CONVECTIONPOINTMOD <   õ  R   a   CONVECTIONPOINTTYPE%INIT+CONVECTIONPOINTMOD (   G  s      INIT+CONVECTIONPOINTMOD -   º  a   a   INIT%THIS+CONVECTIONPOINTMOD 0     @   a   INIT%POINTID+CONVECTIONPOINTMOD -   [  @   a   INIT%COEF+CONVECTIONPOINTMOD -     @   a   INIT%TEMP+CONVECTIONPOINTMOD =   Û  S   a   CONVECTIONPOINTTYPE%APPLY+CONVECTIONPOINTMOD )   .  j      APPLY+CONVECTIONPOINTMOD .     a   a   APPLY%THIS+CONVECTIONPOINTMOD 3   ù  T   a   APPLY%STIFFNESS+CONVECTIONPOINTMOD -   M     a   APPLY%RHS+CONVECTIONPOINTMOD !   Ù  U     SPARSE+SPARSEKIT #   .     a   SPARSE%A+SPARSEKIT $   Â     a   SPARSE%AI+SPARSEKIT $   V     a   SPARSE%AJ+SPARSEKIT ,   ê     a   SPARSE%ROWCOUNTER+SPARSEKIT #   ~  H   a   SPARSE%N+SPARSEKIT %   Æ  H   a   SPARSE%NNZ+SPARSEKIT )     H   a   SPARSE%COUNTER+SPARSEKIT )   V  ]   a   SPARSE%TRIPLET+SPARSEKIT "   ³  i      TRIPLET+SPARSEKIT $        a   TRIPLET%A+SPARSEKIT &   °     a   TRIPLET%ROW+SPARSEKIT &   D     a   TRIPLET%COL+SPARSEKIT &   Ø  R   a   SPARSE%INIT+SPARSEKIT    *  e      INIT+SPARSEKIT $     T   a   INIT%THIS+SPARSEKIT #   ã  @   a   INIT%NNZ+SPARSEKIT $   #  @   a   INIT%ROWS+SPARSEKIT (   c  T   a   SPARSE%UPDATE+SPARSEKIT !   ·  e      UPDATE+SPARSEKIT &     T   a   UPDATE%THIS+SPARSEKIT %   p  @   a   UPDATE%NNZ+SPARSEKIT &   °  @   a   UPDATE%ROWS+SPARSEKIT (   ð  T   a   SPARSE%APPEND+SPARSEKIT !   D  m      APPEND+SPARSEKIT &   ±  T   a   APPEND%THIS+SPARSEKIT %     @   a   APPEND%VAL+SPARSEKIT %   E  @   a   APPEND%ROW+SPARSEKIT %     @   a   APPEND%COL+SPARSEKIT )   Å  U   a   SPARSE%MAKECRS+SPARSEKIT "     `      MAKECRS+SPARSEKIT '   z  T   a   MAKECRS%THIS+SPARSEKIT +   Î  @   a   MAKECRS%SORTROWS+SPARSEKIT %      Q   a   SPARSE%GET+SPARSEKIT    _   h      GET+SPARSEKIT #   Ç   T   a   GET%THIS+SPARSEKIT     !  @   a   GET%I+SPARSEKIT     [!  @   a   GET%J+SPARSEKIT (   !  T   a   SPARSE%GETNNZ+SPARSEKIT !   ï!  Z      GETNNZ+SPARSEKIT &   I"  T   a   GETNNZ%THIS+SPARSEKIT &   "  R   a   SPARSE%GETN+SPARSEKIT    ï"  Z      GETN+SPARSEKIT $   I#  T   a   GETN%THIS+SPARSEKIT ,   #  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   õ#  n      PRINTVALUE+SPARSEKIT *   c$  T   a   PRINTVALUE%THIS+SPARSEKIT '   ·$  @   a   PRINTVALUE%I+SPARSEKIT '   ÷$  @   a   PRINTVALUE%J+SPARSEKIT .   7%  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   %  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   Þ%  `      PRINTNONZEROS+SPARSEKIT -   >&  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   &  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   Þ&  V   a   SPARSE%PRINTALL+SPARSEKIT #   4'  `      PRINTALL+SPARSEKIT (   '  T   a   PRINTALL%THIS+SPARSEKIT ,   è'  L   a   PRINTALL%FILENAME+SPARSEKIT 1   4(  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   (  d      DELETEROWANDCOL+SPARSEKIT /   õ(  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   I)  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   )  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   É)  R   a   SPARSE%FREE+SPARSEKIT    *  R      FREE+SPARSEKIT $   m*  T   a   FREE%THIS+SPARSEKIT C   Á*  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   +  R      HANDLEDUPLICATES+SPARSEKIT 0   q+  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT !   Å+  U      SPARSE+SPARSEKIT #   -     a   SPARSE%A+SPARSEKIT $   ®-     a   SPARSE%AI+SPARSEKIT $   B.     a   SPARSE%AJ+SPARSEKIT ,   Ö.     a   SPARSE%ROWCOUNTER+SPARSEKIT #   j/  H   a   SPARSE%N+SPARSEKIT %   ²/  H   a   SPARSE%NNZ+SPARSEKIT )   ú/  H   a   SPARSE%COUNTER+SPARSEKIT )   B0  ]   a   SPARSE%TRIPLET+SPARSEKIT &   0  R   a   SPARSE%INIT+SPARSEKIT (   ñ0  T   a   SPARSE%UPDATE+SPARSEKIT (   E1  T   a   SPARSE%APPEND+SPARSEKIT )   1  U   a   SPARSE%MAKECRS+SPARSEKIT %   î1  Q   a   SPARSE%GET+SPARSEKIT (   ?2  T   a   SPARSE%GETNNZ+SPARSEKIT &   2  R   a   SPARSE%GETN+SPARSEKIT ,   å2  X   a   SPARSE%PRINTVALUE+SPARSEKIT /   =3  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT *   3  V   a   SPARSE%PRINTALL+SPARSEKIT 1   î3  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT &   K4  R   a   SPARSE%FREE+SPARSEKIT C   4  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES *   û4  l      THBOUNDARYCONDITION1DTYPE 9   g6  ¬   a   THBOUNDARYCONDITION1DTYPE%DIRICHLETPOINT :   7  ­   a   THBOUNDARYCONDITION1DTYPE%NORMALFLUXPOINT :   À7  ­   a   THBOUNDARYCONDITION1DTYPE%CONVECTIONPOINT /   m8  R   a   THBOUNDARYCONDITION1DTYPE%INIT    ¿8        INIT    C9  g   a   INIT%THIS     ª9  @   a   INIT%NDIRICHLET !   ê9  @   a   INIT%NNORMALFLUX !   *:  @   a   INIT%NCONVECTION <   j:  _   a   THBOUNDARYCONDITION1DTYPE%ADDDIRICHLETPOINT "   É:  e      ADDDIRICHLETPOINT '   .;  g   a   ADDDIRICHLETPOINT%THIS %   ;  @   a   ADDDIRICHLETPOINT%ID (   Õ;  @   a   ADDDIRICHLETPOINT%VALUE <   <  _   a   THBOUNDARYCONDITION1DTYPE%GETDIRICHLETPOINT "   t<  y      GETDIRICHLETPOINT '   í<  g   a   GETDIRICHLETPOINT%THIS $   T=  @   a   GETDIRICHLETPOINT%I 8   =  [   a   THBOUNDARYCONDITION1DTYPE%GETNDIRICHLET    ï=  Z      GETNDIRICHLET #   I>  g   a   GETNDIRICHLET%THIS =   °>  `   a   THBOUNDARYCONDITION1DTYPE%ADDNORMALFLUXPOINT #   ?  j      ADDNORMALFLUXPOINT (   z?  g   a   ADDNORMALFLUXPOINT%THIS +   á?  @   a   ADDNORMALFLUXPOINT%POINTID )   !@  @   a   ADDNORMALFLUXPOINT%VALUE =   a@  `   a   THBOUNDARYCONDITION1DTYPE%GETNORMALFLUXPOINT #   Á@  z      GETNORMALFLUXPOINT (   ;A  g   a   GETNORMALFLUXPOINT%THIS %   ¢A  @   a   GETNORMALFLUXPOINT%I 9   âA  \   a   THBOUNDARYCONDITION1DTYPE%GETNNORMALFLUX    >B  Z      GETNNORMALFLUX $   B  g   a   GETNNORMALFLUX%THIS =   ÿB  `   a   THBOUNDARYCONDITION1DTYPE%ADDCONVECTIONPOINT #   _C  s      ADDCONVECTIONPOINT (   ÒC  g   a   ADDCONVECTIONPOINT%THIS +   9D  @   a   ADDCONVECTIONPOINT%POINTID (   yD  @   a   ADDCONVECTIONPOINT%COEF (   ¹D  @   a   ADDCONVECTIONPOINT%TEMP =   ùD  `   a   THBOUNDARYCONDITION1DTYPE%GETCONVECTIONPOINT #   YE  z      GETCONVECTIONPOINT (   ÓE  g   a   GETCONVECTIONPOINT%THIS %   :F  @   a   GETCONVECTIONPOINT%I 9   zF  \   a   THBOUNDARYCONDITION1DTYPE%GETNCONVECTION    ÖF  Z      GETNCONVECTION $   0G  g   a   GETNCONVECTION%THIS 0   G  S   a   THBOUNDARYCONDITION1DTYPE%APPLY    êG  j      APPLY    TH  g   a   APPLY%THIS     »H  T   a   APPLY%STIFFNESS    I     a   APPLY%RHS 