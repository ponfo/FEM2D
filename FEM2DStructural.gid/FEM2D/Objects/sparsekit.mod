  �:  �   k820309    9          19.0        XIC^                                                                                                          
       /home/marco/Gid/problemtypes/FEM2D/FEM2DStructural.gid/FEM2D/Source/Lib/SparseKit.f90 SPARSEKIT              SPARSE INVERSEGMRESD ID i@ i@ i@ gen@SPARSE gen@TRANSPOSE gen@NORM gen@TRACE gen@DET gen@LCHOLESKY gen@GMRES gen@INVERSE gen@JACOBIEIGEN gen@BICGOMP gen@BICGRAD gen@CGOMP                                                     
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
     SIZE    �   h      fn#fn      �   b   uapp(SPARSEKIT    �  @   J  TOOLS    	  @   J  QUICKSORTMOD    I  �      i@ #   �  j      SPARSE_SPARSE_PROD %   7  T   a   SPARSE_SPARSE_PROD%A %   �  T   a   SPARSE_SPARSE_PROD%B !   �  w     SPARSE_VECT_PROD %   V  T   a   SPARSE_VECT_PROD%MAT &   �  �   a   SPARSE_VECT_PROD%VECT !   6  o      COEF_SPARSE_PROD &   �  @   a   COEF_SPARSE_PROD%COEF %   �  T   a   COEF_SPARSE_PROD%MAT    9  W      i@ "   �  j      SPARSE_SPARSE_ADD $   �  T   a   SPARSE_SPARSE_ADD%A $   N  T   a   SPARSE_SPARSE_ADD%B    �  W      i@ "   �  j      SPARSE_SPARSE_SUB $   c	  T   a   SPARSE_SPARSE_SUB%A $   �	  T   a   SPARSE_SPARSE_SUB%B    
  Q       gen@SPARSE    \
  o      CONSTRUCTOR     �
  @   a   CONSTRUCTOR%NNZ !     @   a   CONSTRUCTOR%ROWS    K  O       gen@TRANSPOSE    �  J       gen@NORM    �  K       gen@TRACE    /  I       gen@DET    x  O       gen@LCHOLESKY    �  K       gen@GMRES      M       gen@INVERSE     _  Q       gen@JACOBIEIGEN    �  M       gen@BICGOMP    �  M       gen@BICGRAD    J  K       gen@CGOMP    �  U      SPARSE    �  �   a   SPARSE%A    ~  �   a   SPARSE%AI      �   a   SPARSE%AJ "   �  �   a   SPARSE%ROWCOUNTER    :  H   a   SPARSE%N    �  H   a   SPARSE%NNZ    �  H   a   SPARSE%COUNTER      ]   a   SPARSE%TRIPLET    o  i       TRIPLET    �  �   a   TRIPLET%A    l  �   a   TRIPLET%ROW       �   a   TRIPLET%COL    �  R   a   SPARSE%INIT    �  e      INIT    K  T   a   INIT%THIS    �  @   a   INIT%NNZ    �  @   a   INIT%ROWS      T   a   SPARSE%UPDATE    s  e      UPDATE    �  T   a   UPDATE%THIS    ,  @   a   UPDATE%NNZ    l  @   a   UPDATE%ROWS    �  T   a   SPARSE%APPEND       m      APPEND    m  T   a   APPEND%THIS    �  @   a   APPEND%VAL      @   a   APPEND%ROW    A  @   a   APPEND%COL    �  U   a   SPARSE%MAKECRS    �  `      MAKECRS    6  T   a   MAKECRS%THIS !   �  @   a   MAKECRS%SORTROWS    �  Q   a   SPARSE%GET      h      GET    �  T   a   GET%THIS    �  @   a   GET%I      @   a   GET%J    W  T   a   SPARSE%GETNNZ    �  Z      GETNNZ      T   a   GETNNZ%THIS    Y  R   a   SPARSE%GETN    �  Z      GETN      T   a   GETN%THIS "   Y  X   a   SPARSE%PRINTVALUE    �  n      PRINTVALUE        T   a   PRINTVALUE%THIS    s   @   a   PRINTVALUE%I    �   @   a   PRINTVALUE%J $   �   L   a   PRINTVALUE%FILENAME %   ?!  [   a   SPARSE%PRINTNONZEROS    �!  `      PRINTNONZEROS #   �!  T   a   PRINTNONZEROS%THIS '   N"  L   a   PRINTNONZEROS%FILENAME     �"  V   a   SPARSE%PRINTALL    �"  `      PRINTALL    P#  T   a   PRINTALL%THIS "   �#  L   a   PRINTALL%FILENAME '   �#  ]   a   SPARSE%DELETEROWANDCOL     M$  d      DELETEROWANDCOL %   �$  T   a   DELETEROWANDCOL%THIS $   %  @   a   DELETEROWANDCOL%ROW $   E%  @   a   DELETEROWANDCOL%COL    �%  R   a   SPARSE%FREE    �%  R      FREE    )&  T   a   FREE%THIS (   }&  ^   !   SPARSE%HANDLEDUPLICATES !   �&  R      HANDLEDUPLICATES &   -'  T   a   HANDLEDUPLICATES%THIS    �'  c       TRANSPOSE    �'  T   a   TRANSPOSE%A    8(  W       NORM    �(  T   a   NORM%A    �(  0      GMRES    *  T   a   GMRES%A    g*  
  a   GMRES%RHS    q+  @   a   GMRES%ITRMAXIN    �+  @   a   GMRES%MRIN    �+  @   a   GMRES%TOLABSIN    1,  @   a   GMRES%TOLRELIN    q,  c       INVERSE    �,  T   a   INVERSE%A    (-  q       JACOBIEIGEN $   �-  T   a   JACOBIEIGEN%A_INPUT %   �-  
  a   JACOBIEIGEN%EIGENVAL %   �.  �  a   JACOBIEIGEN%EIGENVEC    �0  W       TRACE    I1  T   a   TRACE%A    �1  c       INVERSEGMRESD      2  T   a   INVERSEGMRESD%A    T2  c       ID    �2  @   a   ID%N    �2  c       LCHOLESKY    Z3  T   a   LCHOLESKY%A    �3  W       DET    4  T   a   DET%A    Y4  �       BICGOMP    5  T   a   BICGOMP%A    Y5  �   a   BICGOMP%RHS    �5  �       BICGRAD    �6  T   a   BICGRAD%A    57  
  a   BICGRAD%RHS    ?8  �       CGOMP    99  T   a   CGOMP%M    �9  
  a   CGOMP%B &   �:  =      SPARSE_VECT_PROD%SIZE 