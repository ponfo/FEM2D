  �!  ]   k820309    9          19.0        �]                                                                                                          
       /home/facundo/Documents/FEM2DIUA_Thermal/FEM2DThermal.gid/FEM2D/Source/ConvectionPoint.f90 CONVECTIONPOINTMOD              CONVECTIONPOINTTYPE gen@CONVECTIONPOINT                                                     
                            @                              
                                                              u #CONSTRUCTOR    &         @   @                                                      #POINTID    #COEF    #TEMP    #CONVECTIONPOINTTYPE              
  @                                                   
  @                                   
                
  @                                   
                        @               �                '                   #A 	   #AI 
   #AJ    #ROWCOUNTER    #N    #NNZ    #COUNTER    #TRIPLET    #INIT    #UPDATE    #APPEND    #MAKECRS %   #GET )   #GETNNZ .   #GETN 1   #PRINTVALUE 4   #PRINTNONZEROS :   #PRINTALL >   #DELETEROWANDCOL B   #FREE G   #HANDLEDUPLICATES J              � $                             	                              
            &                                                     � $                             
            H                             &                                                     � $                                         �                             &                                                      � $                                         �                             &                                                        � $                                                            � $                                  $                         � $                                  (                         � $                                   �       0             #TRIPLET                  @  @              D                '�                    #A    #ROW    #COL               �                                                            
            &                                                      �                                          H                             &                                                      �                                          �                             &                                           1         �   � $                      �                   	     #INIT    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                   
     #UPDATE    #         @     @                                                #THIS    #NNZ    #ROWS              
                                                   #SPARSE              
                                                      
                                            1         �   � $                      �                        #APPEND     #         @     @                                                 #THIS !   #VAL "   #ROW #   #COL $             
                                !                   #SPARSE              
                                 "     
                
                                 #                     
                                 $           1         �   � $                      �      %                  #MAKECRS &   #         @     @                            &                    #THIS '   #SORTROWS (             
                                '                   #SPARSE              
                                 (           1         �   � $                     �      )                  #GET *   %         @   @                           *                    
       #THIS +   #I ,   #J -             
                                +                   #SPARSE              
                                 ,                     
                                 -           1         �   � $                     �      .                  #GETNNZ /   %         @   @                           /                           #THIS 0             
                                0                   #SPARSE    1         �   � $                     �      1                  #GETN 2   %         @   @                           2                           #THIS 3             
                                3                   #SPARSE    1         �   � $                      �      4                  #PRINTVALUE 5   #         @     @                            5                    #THIS 6   #I 7   #J 8   #FILENAME 9             
                                6                   #SPARSE              
                                 7                     
                                 8                     
                               9                    1 1         �   � $                      �      :              	    #PRINTNONZEROS ;   #         @     @                            ;                    #THIS <   #FILENAME =             
                                <                   #SPARSE              
                               =                    1 1         �   � $                      �      >              
    #PRINTALL ?   #         @     @                            ?                    #THIS @   #FILENAME A             
                                @                   #SPARSE              
                               A                    1 1         �   � $                      �      B                  #DELETEROWANDCOL C   #         @     @                            C                    #THIS D   #ROW E   #COL F             
                                D                   #SPARSE              
                                 E                     
                                 F           1         �   � $                      �      G                  #FREE H   #         @     @                            H                    #THIS I             
                                I                   #SPARSE    1         �   � D                      �      J                  #HANDLEDUPLICATES K   #         @     @                            K                    #THIS L             
                                L                   #SPARSE                      @                                '                    #POINTID M   #COEF N   #TEMP O   #INIT P   #APPLY V                � $                             M                                � $                             N               
                � $                             O               
   1         �   � $                     �      P                  #INIT Q   #         @     @                            Q                    #THIS R   #POINTID S   #COEF T   #TEMP U             
D                                R                    #CONVECTIONPOINTTYPE              
                                 S                     
                                 T     
                
                                 U     
      1         �   � $                      �      V                  #APPLY W   #         @     @                             W                    #THIS X   #STIFFNESS Y   #RHS Z             
                                X                    #CONVECTIONPOINTTYPE              
D                                 Y                   #SPARSE              
D                                Z                   
               &                                              �   v      fn#fn (     8   b   uapp(CONVECTIONPOINTMOD    N  @   J  TOOLS    �  @   J  SPARSEKIT $   �  Q       gen@CONVECTIONPOINT      �      CONSTRUCTOR $   �  @   a   CONSTRUCTOR%POINTID !   �  @   a   CONSTRUCTOR%COEF !   )  @   a   CONSTRUCTOR%TEMP !   i  U      SPARSE+SPARSEKIT #   �  �   a   SPARSE%A+SPARSEKIT $   R  �   a   SPARSE%AI+SPARSEKIT $   �  �   a   SPARSE%AJ+SPARSEKIT ,   z  �   a   SPARSE%ROWCOUNTER+SPARSEKIT #     H   a   SPARSE%N+SPARSEKIT %   V  H   a   SPARSE%NNZ+SPARSEKIT )   �  H   a   SPARSE%COUNTER+SPARSEKIT )   �  ]   a   SPARSE%TRIPLET+SPARSEKIT "   C  i      TRIPLET+SPARSEKIT $   �  �   a   TRIPLET%A+SPARSEKIT &   @	  �   a   TRIPLET%ROW+SPARSEKIT &   �	  �   a   TRIPLET%COL+SPARSEKIT &   h
  R   a   SPARSE%INIT+SPARSEKIT    �
  e      INIT+SPARSEKIT $     T   a   INIT%THIS+SPARSEKIT #   s  @   a   INIT%NNZ+SPARSEKIT $   �  @   a   INIT%ROWS+SPARSEKIT (   �  T   a   SPARSE%UPDATE+SPARSEKIT !   G  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %      @   a   UPDATE%NNZ+SPARSEKIT &   @  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &   A  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %     @   a   APPEND%COL+SPARSEKIT )   U  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   
  T   a   MAKECRS%THIS+SPARSEKIT +   ^  @   a   MAKECRS%SORTROWS+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT    �  h      GET+SPARSEKIT #   W  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT     �  @   a   GET%J+SPARSEKIT (   +  T   a   SPARSE%GETNNZ+SPARSEKIT !     Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   -  R   a   SPARSE%GETN+SPARSEKIT      Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT ,   -  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *   �  T   a   PRINTVALUE%THIS+SPARSEKIT '   G  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   �  L   a   PRINTVALUE%FILENAME+SPARSEKIT /     [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   n  `      PRINTNONZEROS+SPARSEKIT -   �  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   "  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   n  V   a   SPARSE%PRINTALL+SPARSEKIT #   �  `      PRINTALL+SPARSEKIT (   $  T   a   PRINTALL%THIS+SPARSEKIT ,   x  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   !  d      DELETEROWANDCOL+SPARSEKIT /   �  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   �  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .     @   a   DELETEROWANDCOL%COL+SPARSEKIT &   Y  R   a   SPARSE%FREE+SPARSEKIT    �  R      FREE+SPARSEKIT $   �  T   a   FREE%THIS+SPARSEKIT C   Q  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �  R      HANDLEDUPLICATES+SPARSEKIT 0     T   a   HANDLEDUPLICATES%THIS+SPARSEKIT $   U  �       CONVECTIONPOINTTYPE ,   �  H   a   CONVECTIONPOINTTYPE%POINTID )   #  H   a   CONVECTIONPOINTTYPE%COEF )   k  H   a   CONVECTIONPOINTTYPE%TEMP )   �  R   a   CONVECTIONPOINTTYPE%INIT      s      INIT    x  a   a   INIT%THIS    �  @   a   INIT%POINTID      @   a   INIT%COEF    Y  @   a   INIT%TEMP *   �  S   a   CONVECTIONPOINTTYPE%APPLY    �  j      APPLY    V   a   a   APPLY%THIS     �   T   a   APPLY%STIFFNESS    !  �   a   APPLY%RHS 