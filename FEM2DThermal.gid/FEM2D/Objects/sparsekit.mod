  �:  �   k820309    9          19.0        ��]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DThermal.gid/FEM2D/Source/Lib/SparseKit.f90 SPARSEKIT              SPARSE INVERSEGMRESD ID i@ i@ i@ gen@SPARSE gen@TRANSPOSE gen@NORM gen@TRACE gen@DET gen@LCHOLESKY gen@GMRES gen@INVERSE gen@JACOBIEIGEN gen@BICGOMP gen@BICGRAD gen@CGOMP                                                     
                                                           
                                                                  #SPARSE_SPARSE_PROD    #SPARSE_VECT_PROD    #COEF_SPARSE_PROD    &         @    @X                                                     #A    #B    #SPARSE              
  @                                                #SPARSE              
  @                                                #SPARSE    (        `   @X                                                               
    #MAT    #VECT 	   p          H r 
     7
S
O
 p        j            j                                      H r 
     7
S
O
 p        j            j                                                              
                                                   #SPARSE           0  
 @                              	                   
              &                                           &         @    @X                                                     #COEF    #MAT    #SPARSE              
                                      
                
                                                    #SPARSE                                                               #SPARSE_SPARSE_ADD    &         @    @X                                                     #A    #B    #SPARSE              
  @                                                #SPARSE              
                                                   #SPARSE                                                               #SPARSE_SPARSE_SUB    &         @    @X                                                     #A    #B    #SPARSE              
  @                                                #SPARSE              
                                                   #SPARSE                                                           u #CONSTRUCTOR    &         @   @                                                    #NNZ    #ROWS    #SPARSE              
  @                                                   
  @                                                  @    X                                    u #TRANSPOSE             @   X                                    u #NORM             @   X                                    u #TRACE             @   X                                     u #DET             @    X                                    u #LCHOLESKY            `   X                                     u #GMRES             @    X                                     u #INVERSE             @     X                                      u #JACOBIEIGEN            D    X                                     u #BICGOMP            `   X                                     u #BICGRAD             `   X                                    u #CGOMP !                     @               �                '                   #A "   #AI #   #AJ $   #ROWCOUNTER %   #N &   #NNZ '   #COUNTER (   #TRIPLET )   #INIT .   #UPDATE 3   #APPEND 8   #MAKECRS >   #GET B   #GETNNZ G   #GETN J   #PRINTVALUE M   #PRINTNONZEROS S   #PRINTALL W   #DELETEROWANDCOL [   #FREE `   #HANDLEDUPLICATES c              � $                             "                              
            &                                                     � $                             #            H                             &                                                     � $                             $            �                             &                                                      � $                             %            �                             &                                                       � $                             &                              � $                             '     $                         � $                             (     (                         � $                              )     �       0             #TRIPLET *                 @  @              @           *     '�                    #A +   #ROW ,   #COL -              �                              +                              
            &                                                      �                              ,            H                             &                                                      �                              -            �                             &                                           1         �   � $                     �      .             	     #INIT /   #         @     @                            /                    #THIS 0   #NNZ 1   #ROWS 2             
D                                0                   #SPARSE              
                                 1                     
                                 2           1         �   � $                      �      3             
     #UPDATE 4   #         @     @                             4                    #THIS 5   #NNZ 6   #ROWS 7             
D                                5                   #SPARSE              
                                 6                     
                                 7           1         �   � $                     �      8                  #APPEND 9   #         @     @                            9                    #THIS :   #VAL ;   #ROW <   #COL =             
D                                :                   #SPARSE              
                                 ;     
                
                                 <                     
                                 =           1         �   � $                     �      >                  #MAKECRS ?   #         @     @                            ?                    #THIS @   #SORTROWS A             
D @                              @                   #SPARSE              
 @                               A           1         �   � $                    �      B                  #GET C   %         @   @                           C                    
       #THIS D   #I E   #J F             
                                D                   #SPARSE              
                                 E                     
                                 F           1         �   � $                     �      G                  #GETNNZ H   %         @   @                            H                           #THIS I             
                                I                   #SPARSE    1         �   � $                     �      J                  #GETN K   %         @   @                            K                           #THIS L             
                                L                   #SPARSE    1         �   � $                      �      M                  #PRINTVALUE N   #         @     @                             N                    #THIS O   #I P   #J Q   #FILENAME R             
D @                              O                   #SPARSE              
  @                              P                     
  @                              Q                     
 @                             R                    1 1         �   � $                      �      S              	    #PRINTNONZEROS T   #         @     @                             T                    #THIS U   #FILENAME V             
                                U                   #SPARSE              
 @                             V                    1 1         �   � $                      �      W              
    #PRINTALL X   #         @     @                             X                    #THIS Y   #FILENAME Z             
D @                              Y                   #SPARSE              
 @                             Z                    1 1         �   � $                      �      [                  #DELETEROWANDCOL \   #         @     @                             \                    #THIS ]   #ROW ^   #COL _             
D                                ]                   #SPARSE              
                                 ^                     
                                 _           1         �   � $                     �      `                  #FREE a   #         @     @                            a                    #THIS b             
D                                b                   #SPARSE    1         �   � D                     �      c                  #HANDLEDUPLICATES d   #         @     @                            d                    #THIS e             
D @                              e                   #SPARSE    &         @    X                                                     #A f   #SPARSE             
  @                              f                  #SPARSE    %         @   X                                               
       #A g             
                                 g                  #SPARSE    (        `   X                                                                
    #A h   #RHS i   #ITRMAXIN j   #MRIN k   #TOLABSIN l   #TOLRELIN m   p          5 8 O#SPARSE     p        U     &       5 8 O#SPARSE     p        U     &                              
@ @                              h                  #SPARSE             
                                 i                    
    p          5 8 �#SPARSE     p        r#SPARSE     h   U     &       5 8 �#SPARSE     p        r#SPARSE     h   U     &                               
 @   �                          j                      
 @   �                          k                      
 @   �                          l     
                 
 @   �                          m     
       &         @    X                                                      #A n   #SPARSE             
  @                              n                  #SPARSE    #         @     X                                                 #A_INPUT o   #EIGENVAL p   #EIGENVEC q            
                                 o                  #SPARSE             D                                p                    
 "    p          5 8 �#SPARSE     p        r#SPARSE     o   U     &       5 8 �#SPARSE     p        r#SPARSE     o   U     &                              D                                q                    
 #      p        5 8 �#SPARSE     p        r#SPARSE     o   U     &   p          5 8 �#SPARSE     p        r#SPARSE     o   U     &     5 8 �#SPARSE     p        r#SPARSE     o   U     &       5 8 �#SPARSE     p        r#SPARSE     o   U     &     5 8 �#SPARSE     p        r#SPARSE     o   U     &                     %         @   X                                               
       #A r             
                                 r                  #SPARSE    &         @                                 s                          #A t   #SPARSE              
  @                              t                  #SPARSE    &         @                                u                          #N v   #SPARSE              
  @                              v           &         @    X                                                     #A w   #SPARSE             
D @                              w                   #SPARSE    %         @   X                                                
       #A x             
D @                              x                   #SPARSE    (        D    X                                               B                
    #A y   #RHS z             &                                                    
@ @                              y                  #SPARSE              
@ @                              z                   
 A             &                                           (        `   X                                                 l                
    #A {   #RHS |   p          5 8 O#SPARSE     p        U     &       5 8 O#SPARSE     p        U     &                              
@ @                              {                  #SPARSE             
@ @                              |                    
 h   p          5 8 �#SPARSE     p        r#SPARSE     {   U     &       5 8 �#SPARSE     p        r#SPARSE     {   U     &                     (        `   X                           !                    �                
    #M }   #B ~   p          5 8 O#SPARSE     p        U     &       5 8 O#SPARSE     p        U     &                              
@ @                              }                  #SPARSE             
@ @                              ~                    
 �   p          5 8 �#SPARSE     p        r#SPARSE     }   U     &       5 8 �#SPARSE     p        r#SPARSE     }   U     &                                   @                           
     SIZE    �   r      fn#fn      �   b   uapp(SPARSEKIT    �  @   J  TOOLS      @   J  QUICKSORTMOD    S  �      i@ #   �  j      SPARSE_SPARSE_PROD %   A  T   a   SPARSE_SPARSE_PROD%A %   �  T   a   SPARSE_SPARSE_PROD%B !   �  w     SPARSE_VECT_PROD %   `  T   a   SPARSE_VECT_PROD%MAT &   �  �   a   SPARSE_VECT_PROD%VECT !   @  o      COEF_SPARSE_PROD &   �  @   a   COEF_SPARSE_PROD%COEF %   �  T   a   COEF_SPARSE_PROD%MAT    C  W      i@ "   �  j      SPARSE_SPARSE_ADD $     T   a   SPARSE_SPARSE_ADD%A $   X  T   a   SPARSE_SPARSE_ADD%B    �  W      i@ "   	  j      SPARSE_SPARSE_SUB $   m	  T   a   SPARSE_SPARSE_SUB%A $   �	  T   a   SPARSE_SPARSE_SUB%B    
  Q       gen@SPARSE    f
  o      CONSTRUCTOR     �
  @   a   CONSTRUCTOR%NNZ !     @   a   CONSTRUCTOR%ROWS    U  O       gen@TRANSPOSE    �  J       gen@NORM    �  K       gen@TRACE    9  I       gen@DET    �  O       gen@LCHOLESKY    �  K       gen@GMRES      M       gen@INVERSE     i  Q       gen@JACOBIEIGEN    �  M       gen@BICGOMP      M       gen@BICGRAD    T  K       gen@CGOMP    �  U      SPARSE    �  �   a   SPARSE%A    �  �   a   SPARSE%AI      �   a   SPARSE%AJ "   �  �   a   SPARSE%ROWCOUNTER    D  H   a   SPARSE%N    �  H   a   SPARSE%NNZ    �  H   a   SPARSE%COUNTER      ]   a   SPARSE%TRIPLET    y  i       TRIPLET    �  �   a   TRIPLET%A    v  �   a   TRIPLET%ROW    
  �   a   TRIPLET%COL    �  R   a   SPARSE%INIT    �  e      INIT    U  T   a   INIT%THIS    �  @   a   INIT%NNZ    �  @   a   INIT%ROWS    )  T   a   SPARSE%UPDATE    }  e      UPDATE    �  T   a   UPDATE%THIS    6  @   a   UPDATE%NNZ    v  @   a   UPDATE%ROWS    �  T   a   SPARSE%APPEND    
  m      APPEND    w  T   a   APPEND%THIS    �  @   a   APPEND%VAL      @   a   APPEND%ROW    K  @   a   APPEND%COL    �  U   a   SPARSE%MAKECRS    �  `      MAKECRS    @  T   a   MAKECRS%THIS !   �  @   a   MAKECRS%SORTROWS    �  Q   a   SPARSE%GET    %  h      GET    �  T   a   GET%THIS    �  @   a   GET%I    !  @   a   GET%J    a  T   a   SPARSE%GETNNZ    �  Z      GETNNZ      T   a   GETNNZ%THIS    c  R   a   SPARSE%GETN    �  Z      GETN      T   a   GETN%THIS "   c  X   a   SPARSE%PRINTVALUE    �  n      PRINTVALUE     )   T   a   PRINTVALUE%THIS    }   @   a   PRINTVALUE%I    �   @   a   PRINTVALUE%J $   �   L   a   PRINTVALUE%FILENAME %   I!  [   a   SPARSE%PRINTNONZEROS    �!  `      PRINTNONZEROS #   "  T   a   PRINTNONZEROS%THIS '   X"  L   a   PRINTNONZEROS%FILENAME     �"  V   a   SPARSE%PRINTALL    �"  `      PRINTALL    Z#  T   a   PRINTALL%THIS "   �#  L   a   PRINTALL%FILENAME '   �#  ]   a   SPARSE%DELETEROWANDCOL     W$  d      DELETEROWANDCOL %   �$  T   a   DELETEROWANDCOL%THIS $   %  @   a   DELETEROWANDCOL%ROW $   O%  @   a   DELETEROWANDCOL%COL    �%  R   a   SPARSE%FREE    �%  R      FREE    3&  T   a   FREE%THIS (   �&  ^   !   SPARSE%HANDLEDUPLICATES !   �&  R      HANDLEDUPLICATES &   7'  T   a   HANDLEDUPLICATES%THIS    �'  c       TRANSPOSE    �'  T   a   TRANSPOSE%A    B(  W       NORM    �(  T   a   NORM%A    �(  0      GMRES    *  T   a   GMRES%A    q*  
  a   GMRES%RHS    {+  @   a   GMRES%ITRMAXIN    �+  @   a   GMRES%MRIN    �+  @   a   GMRES%TOLABSIN    ;,  @   a   GMRES%TOLRELIN    {,  c       INVERSE    �,  T   a   INVERSE%A    2-  q       JACOBIEIGEN $   �-  T   a   JACOBIEIGEN%A_INPUT %   �-  
  a   JACOBIEIGEN%EIGENVAL %   /  �  a   JACOBIEIGEN%EIGENVEC    �0  W       TRACE    S1  T   a   TRACE%A    �1  c       INVERSEGMRESD     
2  T   a   INVERSEGMRESD%A    ^2  c       ID    �2  @   a   ID%N    3  c       LCHOLESKY    d3  T   a   LCHOLESKY%A    �3  W       DET    4  T   a   DET%A    c4  �       BICGOMP    5  T   a   BICGOMP%A    c5  �   a   BICGOMP%RHS    �5  �       BICGRAD    �6  T   a   BICGRAD%A    ?7  
  a   BICGRAD%RHS    I8  �       CGOMP    C9  T   a   CGOMP%M    �9  
  a   CGOMP%B &   �:  =      SPARSE_VECT_PROD%SIZE 