  h$  j   k820309    9          19.0        ��]                                                                                                          
       /home/marco/Gid/problemtypes/FEM2D/FEM2DThermal.gid/FEM2D/Source/Domain/Domain.f90 DOMAINMOD              DOMAINTYPE                                                     
                       @�                                  
                            @                              
                            @                              
                      �  @                               '                     #ID    #X    #Y    #Z 	   #INIT 
   #SETID    #GETID    #SETX    #GETX    #SETY    #GETY #   #SETZ &   #GETZ *                � D                                                             � D                                            
                � D                                            
                � D                             	               
   1         �   � $                      �      
                  #INIT    #         @     @                                                #THIS    #ID    #X    #Y    #Z              
                                                     #POINTTYPE              
                                                      
                                      
                
                                      
                
                                      
      1         �   � $                      �                        #SETID    #         @     @                                                #THIS    #ID              
                                                     #POINTTYPE              
                                            1         �   � $                     �                        #GETID    %         @   @                                                      #THIS              
                                                     #POINTTYPE    1         �   � $                      �                        #SETX    #         @     @                                                #THIS    #X              
                                                     #POINTTYPE              
                                      
      1         �   � $                     �                   	     #GETX    %         @   @                                               
       #THIS              
                                                     #POINTTYPE    1         �   � $                      �                   
     #SETY     #         @     @                                                 #THIS !   #Y "             
                                !                     #POINTTYPE              
                                 "     
      1         �   � $                     �      #                  #GETY $   %         @   @                           $                    
       #THIS %             
                                %                     #POINTTYPE    1         �   � $                      �      &                  #SETZ '   #         @     @                            '                    #THIS (   #Z )             
                                (                     #POINTTYPE              
                                 )     
      1         �   � $                     �      *              	    #GETZ +   %         @   @                           +                    
       #THIS ,             
                                ,                     #POINTTYPE                   @  @                           -     '              
      #PTR .   #ALLOCATE /   #SETID 3   #GETID 7   #SETX :   #GETX >   #SETY A   #GETY E   #SETZ H   #GETZ L                �$                              .                            #POINTTYPE    1         �   � $                     �      /                  #ALLOCATE 0   #         @     @                           0                    #THIS 1   #POINT 2             
                                1                    #POINTPTRTYPE -             
                                  2                    #POINTTYPE    1         �   � $                      �      3                  #SETID 4   #         @     @                            4                    #THIS 5   #ID 6             
                                5                    #POINTPTRTYPE -             
                                 6           1         �   � $                     �      7                  #GETID 8   %         @   @                           8                           #THIS 9             
                                9                    #POINTPTRTYPE -   1         �   � $                      �      :                  #SETX ;   #         @     @                            ;                    #THIS <   #X =             
                                <                    #POINTPTRTYPE -             
                                 =     
      1         �   � $                     �      >                  #GETX ?   %         @   @                           ?                    
       #THIS @             
                                @                    #POINTPTRTYPE -   1         �   � $                      �      A                  #SETY B   #         @     @                            B                    #THIS C   #Y D             
                                C                    #POINTPTRTYPE -             
                                 D     
      1         �   � $                     �      E                  #GETY F   %         @   @                           F                    
       #THIS G             
                                G                    #POINTPTRTYPE -   1         �   � $                      �      H             	     #SETZ I   #         @     @                            I                    #THIS J   #Z K             
                                J                    #POINTPTRTYPE -             
                                 K     
      1         �   � $                     �      L             
 	    #GETZ M   %         @   @                           M                    
       #THIS N             
                                N                    #POINTPTRTYPE -                     @               @           O     '`                    #NPOINT P   #NLINE Q   #NTRIANG R   #NQUAD S   #NELEM T   #POINT U   #GETNPOINT V   #GETNLINE Y   #GETNTRIANG \   #GETNQUAD _   #GETNELEM b   #GETPOINT e                � $                             P                                � $                             Q                               � $                             R                               � $                             S                               � $                             T                             � $                              U                                 #POINTTYPE              &                                           1         �   � $                     �      V                  #GETNPOINT W   %         @   @                            W                           #THIS X             
                                X     `               #DOMAINTYPE O   1         �   � $                     �      Y                  #GETNLINE Z   %         @   @                            Z                           #THIS [             
                                [     `               #DOMAINTYPE O   1         �   � $                     �      \             	     #GETNTRIANG ]   %         @   @                            ]                           #THIS ^             
                                ^     `               #DOMAINTYPE O   1         �   � $                     �      _             
     #GETNQUAD `   %         @   @                            `                           #THIS a             
                                a     `               #DOMAINTYPE O   1         �   � $                     �      b                  #GETNELEM c   %         @   @                            c                           #THIS d             
                                d     `               #DOMAINTYPE O   1         �   � $                     �      e                  #GETPOINT f   &         @   @                            f                           #THIS g   #IDX h   #POINTPTRTYPE -             
                                g     `               #DOMAINTYPE O             
                                 h              �   e      fn#fn         b   uapp(DOMAINMOD       @   J  TOOLS    `  @   J  DEBUGGERMOD    �  @   J  POINTMOD    �  @   J  POINTPTRMOD #      �       POINTTYPE+POINTMOD )   �  H   %   POINTTYPE%ID+POINTMOD=ID '   1  H   %   POINTTYPE%X+POINTMOD=X '   y  H   %   POINTTYPE%Y+POINTMOD=Y '   �  H   %   POINTTYPE%Z+POINTMOD=Z (   	  R   a   POINTTYPE%INIT+POINTMOD    [  o      INIT+POINTMOD #   �  W   a   INIT%THIS+POINTMOD !   !  @   a   INIT%ID+POINTMOD     a  @   a   INIT%X+POINTMOD     �  @   a   INIT%Y+POINTMOD     �  @   a   INIT%Z+POINTMOD )   !  S   a   POINTTYPE%SETID+POINTMOD    t  Z      SETID+POINTMOD $   �  W   a   SETID%THIS+POINTMOD "   %  @   a   SETID%ID+POINTMOD )   e  S   a   POINTTYPE%GETID+POINTMOD    �  Z      GETID+POINTMOD $     W   a   GETID%THIS+POINTMOD (   i  R   a   POINTTYPE%SETX+POINTMOD    �  Y      SETX+POINTMOD #   	  W   a   SETX%THIS+POINTMOD     k	  @   a   SETX%X+POINTMOD (   �	  R   a   POINTTYPE%GETX+POINTMOD    �	  Z      GETX+POINTMOD #   W
  W   a   GETX%THIS+POINTMOD (   �
  R   a   POINTTYPE%SETY+POINTMOD       Y      SETY+POINTMOD #   Y  W   a   SETY%THIS+POINTMOD     �  @   a   SETY%Y+POINTMOD (   �  R   a   POINTTYPE%GETY+POINTMOD    B  Z      GETY+POINTMOD #   �  W   a   GETY%THIS+POINTMOD (   �  R   a   POINTTYPE%SETZ+POINTMOD    E  Y      SETZ+POINTMOD #   �  W   a   SETZ%THIS+POINTMOD     �  @   a   SETZ%Z+POINTMOD (   5  R   a   POINTTYPE%GETZ+POINTMOD    �  Z      GETZ+POINTMOD #   �  W   a   GETZ%THIS+POINTMOD )   8  �      POINTPTRTYPE+POINTPTRMOD -   �  _   a   POINTPTRTYPE%PTR+POINTPTRMOD 2   P  V   a   POINTPTRTYPE%ALLOCATE+POINTPTRMOD %   �  ]      ALLOCATE+POINTPTRMOD *     Z   a   ALLOCATE%THIS+POINTPTRMOD +   ]  W   a   ALLOCATE%POINT+POINTPTRMOD /   �  S   a   POINTPTRTYPE%SETID+POINTPTRMOD "     Z      SETID+POINTPTRMOD '   a  Z   a   SETID%THIS+POINTPTRMOD %   �  @   a   SETID%ID+POINTPTRMOD /   �  S   a   POINTPTRTYPE%GETID+POINTPTRMOD "   N  Z      GETID+POINTPTRMOD '   �  Z   a   GETID%THIS+POINTPTRMOD .     R   a   POINTPTRTYPE%SETX+POINTPTRMOD !   T  Y      SETX+POINTPTRMOD &   �  Z   a   SETX%THIS+POINTPTRMOD #     @   a   SETX%X+POINTPTRMOD .   G  R   a   POINTPTRTYPE%GETX+POINTPTRMOD !   �  Z      GETX+POINTPTRMOD &   �  Z   a   GETX%THIS+POINTPTRMOD .   M  R   a   POINTPTRTYPE%SETY+POINTPTRMOD !   �  Y      SETY+POINTPTRMOD &   �  Z   a   SETY%THIS+POINTPTRMOD #   R  @   a   SETY%Y+POINTPTRMOD .   �  R   a   POINTPTRTYPE%GETY+POINTPTRMOD !   �  Z      GETY+POINTPTRMOD &   >  Z   a   GETY%THIS+POINTPTRMOD .   �  R   a   POINTPTRTYPE%SETZ+POINTPTRMOD !   �  Y      SETZ+POINTPTRMOD &   C  Z   a   SETZ%THIS+POINTPTRMOD #   �  @   a   SETZ%Z+POINTPTRMOD .   �  R   a   POINTPTRTYPE%GETZ+POINTPTRMOD !   /  Z      GETZ+POINTPTRMOD &   �  Z   a   GETZ%THIS+POINTPTRMOD    �  �       DOMAINTYPE "   �  H   a   DOMAINTYPE%NPOINT !     H   a   DOMAINTYPE%NLINE #   _  H   a   DOMAINTYPE%NTRIANG !   �  H   a   DOMAINTYPE%NQUAD !   �  H   a   DOMAINTYPE%NELEM !   7  �   a   DOMAINTYPE%POINT %   �  W   a   DOMAINTYPE%GETNPOINT    1  Z      GETNPOINT    �  X   a   GETNPOINT%THIS $   �  V   a   DOMAINTYPE%GETNLINE    9  Z      GETNLINE    �  X   a   GETNLINE%THIS &   �  X   a   DOMAINTYPE%GETNTRIANG    C   Z      GETNTRIANG     �   X   a   GETNTRIANG%THIS $   �   V   a   DOMAINTYPE%GETNQUAD    K!  Z      GETNQUAD    �!  X   a   GETNQUAD%THIS $   �!  V   a   DOMAINTYPE%GETNELEM    S"  Z      GETNELEM    �"  X   a   GETNELEM%THIS $   #  V   a   DOMAINTYPE%GETPOINT    [#  u      GETPOINT    �#  X   a   GETPOINT%THIS    ($  @   a   GETPOINT%IDX 