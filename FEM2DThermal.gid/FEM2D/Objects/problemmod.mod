  �  M   k820309    9          19.0        !�]                                                                                                          
       /home/facundo/Documents/FEM2DIUA_Thermal/FEM2DThermal.gid/FEM2D/Source/Problem.f90 PROBLEMMOD              PROBLEMTYPE                                                     
                            @                              
                         @              �                '                   #A    #AI    #AJ    #ROWCOUNTER    #N    #NNZ 	   #COUNTER 
   #TRIPLET    #INIT    #UPDATE    #APPEND    #MAKECRS     #GET $   #GETNNZ )   #GETN ,   #PRINTVALUE /   #PRINTNONZEROS 5   #PRINTALL 9   #DELETEROWANDCOL =   #FREE B   #HANDLEDUPLICATES E              � $                                                           
            &                                                     � $                                         H                             &                                                     � $                                         �                             &                                                      � $                                         �                             &                                                        � $                                                            � $                             	     $                         � $                             
     (                         � $                                   �       0             #TRIPLET                  @  @              D                '�                    #A    #ROW    #COL               �                                                            
            &                                                      �                                          H                             &                                                      �                                          �                             &                                           1         �   � $                      �                   	     #INIT    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                   
     #UPDATE    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                        #APPEND    #         @     @                                                #THIS    #VAL    #ROW    #COL              
                                                   #SPARSE              
                                      
                
                                                      
                                            1         �   � $                      �                         #MAKECRS !   #         @     @                            !                    #THIS "   #SORTROWS #             
                                "                   #SPARSE              
                                 #           1         �   � $                     �      $                  #GET %   %         @   @                           %                    
       #THIS &   #I '   #J (             
                                &                   #SPARSE              
                                 '                     
                                 (           1         �   � $                     �      )                  #GETNNZ *   %         @   @                           *                           #THIS +             
                                +                   #SPARSE    1         �   � $                     �      ,                  #GETN -   %         @   @                           -                           #THIS .             
                                .                   #SPARSE    1         �   � $                      �      /                  #PRINTVALUE 0   #         @     @                            0                    #THIS 1   #I 2   #J 3   #FILENAME 4             
                                1                   #SPARSE              
                                 2                     
                                 3                     
                               4                    1 1         �   � $                      �      5              	    #PRINTNONZEROS 6   #         @     @                            6                    #THIS 7   #FILENAME 8             
                                7                   #SPARSE              
                               8                    1 1         �   � $                      �      9              
    #PRINTALL :   #         @     @                            :                    #THIS ;   #FILENAME <             
                                ;                   #SPARSE              
                               <                    1 1         �   � $                      �      =                  #DELETEROWANDCOL >   #         @     @                            >                    #THIS ?   #ROW @   #COL A             
                                ?                   #SPARSE              
                                 @                     
                                 A           1         �   � $                      �      B                  #FREE C   #         @     @                            C                    #THIS D             
                                D                   #SPARSE    1         �   � D                      �      E                  #HANDLEDUPLICATES F   #         @     @                            F                    #THIS G             
                                G                   #SPARSE                      @               �           H     '�                   #STIFFNESS I   #RHS J   #DOF K                � $                              I                          #SPARSE               � $                             J                            
            &                                                      � $                             K            P                
            &                                              �   f      fn#fn          b   uapp(PROBLEMMOD    "  @   J  TOOLS    b  @   J  SPARSEKIT !   �  U      SPARSE+SPARSEKIT #   �  �   a   SPARSE%A+SPARSEKIT $   �  �   a   SPARSE%AI+SPARSEKIT $     �   a   SPARSE%AJ+SPARSEKIT ,   �  �   a   SPARSE%ROWCOUNTER+SPARSEKIT #   G  H   a   SPARSE%N+SPARSEKIT %   �  H   a   SPARSE%NNZ+SPARSEKIT )   �  H   a   SPARSE%COUNTER+SPARSEKIT )     ]   a   SPARSE%TRIPLET+SPARSEKIT "   |  i      TRIPLET+SPARSEKIT $   �  �   a   TRIPLET%A+SPARSEKIT &   y  �   a   TRIPLET%ROW+SPARSEKIT &     �   a   TRIPLET%COL+SPARSEKIT &   �  R   a   SPARSE%INIT+SPARSEKIT    �  e      INIT+SPARSEKIT $   X	  T   a   INIT%THIS+SPARSEKIT #   �	  @   a   INIT%NNZ+SPARSEKIT $   �	  @   a   INIT%ROWS+SPARSEKIT (   ,
  T   a   SPARSE%UPDATE+SPARSEKIT !   �
  e      UPDATE+SPARSEKIT &   �
  T   a   UPDATE%THIS+SPARSEKIT %   9  @   a   UPDATE%NNZ+SPARSEKIT &   y  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !     m      APPEND+SPARSEKIT &   z  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %     @   a   APPEND%ROW+SPARSEKIT %   N  @   a   APPEND%COL+SPARSEKIT )   �  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   C  T   a   MAKECRS%THIS+SPARSEKIT +   �  @   a   MAKECRS%SORTROWS+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT    (  h      GET+SPARSEKIT #   �  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT     $  @   a   GET%J+SPARSEKIT (   d  T   a   SPARSE%GETNNZ+SPARSEKIT !   �  Z      GETNNZ+SPARSEKIT &     T   a   GETNNZ%THIS+SPARSEKIT &   f  R   a   SPARSE%GETN+SPARSEKIT    �  Z      GETN+SPARSEKIT $     T   a   GETN%THIS+SPARSEKIT ,   f  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *   ,  T   a   PRINTVALUE%THIS+SPARSEKIT '   �  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .      L   a   PRINTVALUE%FILENAME+SPARSEKIT /   L  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �  `      PRINTNONZEROS+SPARSEKIT -     T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   [  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �  V   a   SPARSE%PRINTALL+SPARSEKIT #   �  `      PRINTALL+SPARSEKIT (   ]  T   a   PRINTALL%THIS+SPARSEKIT ,   �  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   Z  d      DELETEROWANDCOL+SPARSEKIT /   �  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .     @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   R  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   �  R   a   SPARSE%FREE+SPARSEKIT    �  R      FREE+SPARSEKIT $   6  T   a   FREE%THIS+SPARSEKIT C   �  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �  R      HANDLEDUPLICATES+SPARSEKIT 0   :  T   a   HANDLEDUPLICATES%THIS+SPARSEKIT    �  q       PROBLEMTYPE &   �  \   a   PROBLEMTYPE%STIFFNESS     [  �   a   PROBLEMTYPE%RHS     �  �   a   PROBLEMTYPE%DOF 