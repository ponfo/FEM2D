  �!  ]   k820309    9          19.0        }$�]                                                                                                          
       /home/facundo/Documents/FEM2D-LastVersion/FEM2D/FEM2DThermal.gid/FEM2D/Source/BoundaryCondition/ConvectionPoint.f90 CONVECTIONPOINTMOD              CONVECTIONPOINTTYPE gen@CONVECTIONPOINT                                                     
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
               &                                              �   �      fn#fn (   /  8   b   uapp(CONVECTIONPOINTMOD    g  @   J  TOOLS    �  @   J  SPARSEKIT $   �  Q       gen@CONVECTIONPOINT    8  �      CONSTRUCTOR $   �  @   a   CONSTRUCTOR%POINTID !     @   a   CONSTRUCTOR%COEF !   B  @   a   CONSTRUCTOR%TEMP !   �  U      SPARSE+SPARSEKIT #   �  �   a   SPARSE%A+SPARSEKIT $   k  �   a   SPARSE%AI+SPARSEKIT $   �  �   a   SPARSE%AJ+SPARSEKIT ,   �  �   a   SPARSE%ROWCOUNTER+SPARSEKIT #   '  H   a   SPARSE%N+SPARSEKIT %   o  H   a   SPARSE%NNZ+SPARSEKIT )   �  H   a   SPARSE%COUNTER+SPARSEKIT )   �  ]   a   SPARSE%TRIPLET+SPARSEKIT "   \  i      TRIPLET+SPARSEKIT $   �  �   a   TRIPLET%A+SPARSEKIT &   Y	  �   a   TRIPLET%ROW+SPARSEKIT &   �	  �   a   TRIPLET%COL+SPARSEKIT &   �
  R   a   SPARSE%INIT+SPARSEKIT    �
  e      INIT+SPARSEKIT $   8  T   a   INIT%THIS+SPARSEKIT #   �  @   a   INIT%NNZ+SPARSEKIT $   �  @   a   INIT%ROWS+SPARSEKIT (     T   a   SPARSE%UPDATE+SPARSEKIT !   `  e      UPDATE+SPARSEKIT &   �  T   a   UPDATE%THIS+SPARSEKIT %     @   a   UPDATE%NNZ+SPARSEKIT &   Y  @   a   UPDATE%ROWS+SPARSEKIT (   �  T   a   SPARSE%APPEND+SPARSEKIT !   �  m      APPEND+SPARSEKIT &   Z  T   a   APPEND%THIS+SPARSEKIT %   �  @   a   APPEND%VAL+SPARSEKIT %   �  @   a   APPEND%ROW+SPARSEKIT %   .  @   a   APPEND%COL+SPARSEKIT )   n  U   a   SPARSE%MAKECRS+SPARSEKIT "   �  `      MAKECRS+SPARSEKIT '   #  T   a   MAKECRS%THIS+SPARSEKIT +   w  @   a   MAKECRS%SORTROWS+SPARSEKIT %   �  Q   a   SPARSE%GET+SPARSEKIT      h      GET+SPARSEKIT #   p  T   a   GET%THIS+SPARSEKIT     �  @   a   GET%I+SPARSEKIT       @   a   GET%J+SPARSEKIT (   D  T   a   SPARSE%GETNNZ+SPARSEKIT !   �  Z      GETNNZ+SPARSEKIT &   �  T   a   GETNNZ%THIS+SPARSEKIT &   F  R   a   SPARSE%GETN+SPARSEKIT    �  Z      GETN+SPARSEKIT $   �  T   a   GETN%THIS+SPARSEKIT ,   F  X   a   SPARSE%PRINTVALUE+SPARSEKIT %   �  n      PRINTVALUE+SPARSEKIT *     T   a   PRINTVALUE%THIS+SPARSEKIT '   `  @   a   PRINTVALUE%I+SPARSEKIT '   �  @   a   PRINTVALUE%J+SPARSEKIT .   �  L   a   PRINTVALUE%FILENAME+SPARSEKIT /   ,  [   a   SPARSE%PRINTNONZEROS+SPARSEKIT (   �  `      PRINTNONZEROS+SPARSEKIT -   �  T   a   PRINTNONZEROS%THIS+SPARSEKIT 1   ;  L   a   PRINTNONZEROS%FILENAME+SPARSEKIT *   �  V   a   SPARSE%PRINTALL+SPARSEKIT #   �  `      PRINTALL+SPARSEKIT (   =  T   a   PRINTALL%THIS+SPARSEKIT ,   �  L   a   PRINTALL%FILENAME+SPARSEKIT 1   �  ]   a   SPARSE%DELETEROWANDCOL+SPARSEKIT *   :  d      DELETEROWANDCOL+SPARSEKIT /   �  T   a   DELETEROWANDCOL%THIS+SPARSEKIT .   �  @   a   DELETEROWANDCOL%ROW+SPARSEKIT .   2  @   a   DELETEROWANDCOL%COL+SPARSEKIT &   r  R   a   SPARSE%FREE+SPARSEKIT    �  R      FREE+SPARSEKIT $     T   a   FREE%THIS+SPARSEKIT C   j  ^   %   SPARSE%HANDLEDUPLICATES+SPARSEKIT=HANDLEDUPLICATES +   �  R      HANDLEDUPLICATES+SPARSEKIT 0     T   a   HANDLEDUPLICATES%THIS+SPARSEKIT $   n  �       CONVECTIONPOINTTYPE ,   �  H   a   CONVECTIONPOINTTYPE%POINTID )   <  H   a   CONVECTIONPOINTTYPE%COEF )   �  H   a   CONVECTIONPOINTTYPE%TEMP )   �  R   a   CONVECTIONPOINTTYPE%INIT      s      INIT    �  a   a   INIT%THIS    �  @   a   INIT%POINTID    2  @   a   INIT%COEF    r  @   a   INIT%TEMP *   �  S   a   CONVECTIONPOINTTYPE%APPLY       j      APPLY    o   a   a   APPLY%THIS     �   T   a   APPLY%STIFFNESS    $!  �   a   APPLY%RHS 