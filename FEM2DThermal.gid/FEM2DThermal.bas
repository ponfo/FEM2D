###########################################################
            PROGRAM POISSON'S 1D AND 2D ANALYSIS
###########################################################

##################### PROBLEM DATA ########################

Problem Name: *gendata(Problem_Name)

###################### Mesh Data ##########################
Elements_Number........................: *nelem
Nodes_Number...........................: *npoin
Are_Elements_Quadratic.................: *isQuadratic
*#---------------------------------------------------------
*set var i=0
*set var j=0
*set var k=0
*loop elems
*#ElemsTypeName
*if(strcasecmp(ElemsTypeName(),"Linear")==0)
*set var i=operation(i+1)
*endif
*if(strcasecmp(ElemsTypeName(),"Triangle")==0)
*set var j=operation(j+1)
*endif
*if(strcasecmp(ElemsTypeName(),"Quadrilateral")==0)
*set var k=operation(k+1)
*endif
*end elems
*#---------------------------------------------------------
Linear_Elements_Number.................: *i
Triangular_Elements_Number.............: *j
Rectangular_Elements_Number............: *k
*#---------------------------------------------------------
Materials_Number.......................: *nmats
Gauss_Order............................: *GenData(Gauss_Order)
*#---------------------------------------------------------
*Set Cond Fix_Temperature_On_Lines *nodes
*Set var a = condnumentities
*Set Cond Fix_Temperature_On_Points *nodes
*Set var b = condnumentities
*Set var c = operation(a+b)
*#---------------------------------------------------------
Dirichlet_Conditions_Number............: *c
*#---------------------------------------------------------
*Set Cond Normal_Flux_On_Points *nodes
Normal_Flux_Points_Condition_Nodes.....: *condnumentities
*#---------------------------------------------------------
*Set Cond Convection_On_Points *nodes
Convection_Points_Condition_Nodes......: *condnumentities
*#---------------------------------------------------------
*Set Cond Normal_Flux_On_Lines *elems
Normal_Flux_Line_Condition_Elements....: *condnumentities
*#--------------------------------------------------------- REVISAR ESTAS DOS CONDICIONES DE LINEAR, *ELEMS ESTÁ MAL
*Set Cond Convection_On_Lines *elems
Convection_Lines_Condition_Elements....: *condnumentities
*#---------------------------------------------------------
*Set Cond Source_On_Points *nodes
Points_With_Point_Source...............: *condnumentities
*#---------------------------------------------------------
Source_Number_On_Points................: *Gendata(Source_Number_On_Points,int)
*#---------------------------------------------------------
*Set Cond Source_On_Lines *nodes
*#---------------------------------------------------------
Points_With_Line_Source.................: *condnumentities
*#---------------------------------------------------------
Source_Number_On_Lines.................: *Gendata(Source_Number_On_Lines,int)
*#---------------------------------------------------------
*Set Cond Source_On_Surfaces *elems
Surfaces_With_Surface_Source...........: *condnumentities
*#---------------------------------------------------------
Source_Number_On_Surfaces..............: *Gendata(Source_Number_On_Surfaces,int)
*#---------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coordinates:

  Node   |     X       |       Y        |        Z         |
------------------------------------------------------------  
*set elems(all)
*loop nodes
*format "%5i%10.4e%10.4e"
*NodesNum       *NodesCoord(1,real)     *NodesCoord(2,real)   *NodesCoord(3,real)
*end nodes

######################## Materials ########################

Materials List:

Material | Thermal conductivity X | Thermal conductivity Y
-----------------------------------------------------------
*loop materials
*format "%5i%10.4e%10.4e"
*matnum            *matprop(Thermal_Conductivity_X)               *matprop(Thermal_Conductivity_Y)
*end

################### Sources On Points ######################

Conditions List:

Source   |  Function
----------------------
*Set cond Source_On_Points *nodes
*for(i=1;i<=Gendata(Source_Number_On_Points,int);i=i+1))
*loop nodes *OnlyInCond
*if(i==cond(Source_Number_On_Points,real))
*format "%5i%10s"
*cond(Source_Number_On_Points) *cond(SourceOP)
*break
*endif
*end
*end for

################### Sources On Lines #######################

Conditions List:

Source   |  Function
----------------------
*Set cond Source_On_Lines *nodes
*for(i=1;i<=Gendata(Source_Number_On_Lines,int);i=i+1))
*loop nodes *OnlyInCond
*if(i==cond(Source_Number_On_Lines,real))
*format "%5i%10s"
*cond(Source_Number_On_Lines) *cond(SourceOL)
*break
*endif
*end
*end for

################## Sources On Surfaces #####################

Conditions List:

Source   |  Function
----------------------
*Set cond Source_On_Surfaces
*for(i=1;i<=Gendata(Source_Number_On_Surfaces,int);i=i+1))
*loop elems
*if(i==cond(Source_Number_On_Surfaces,real))
*format "%5i%10s"
*cond(Source_Number_On_Surfaces) *cond(SourceOS) 
*break
*endif
*end
*end for

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Element List:

      Element  |      Type      |  Material  |   Nodes  |      Conectivities
-----------------------------------------------------------------------------------
*Set Cond Source_On_Surfaces *elems
*loop elems
*format "%10i%10i%9i%9i%9i"
*elemsnum         *ElemstypeName  *elemsmat  *ElemsNnode  *elemsconec
*end elems

###################### point Sources ######################

Conditions List:

  Node   |  Source
--------------------------
*Set Cond Source_On_Points *nodes
*loop nodes *OnlyInCond
*format "%5i%5i"
*NodesNum      *cond(Source_Number_On_Points) 
*end

###################### line Sources ######################

Conditions List:

  Node   |  Source
--------------------------
*Set Cond Source_On_Lines *nodes
*loop nodes *OnlyInCond
*format "%5i%5i"
*NodesNum     *cond(Source_Number_On_Lines) 
*end

###################### surface Sources ######################

Conditions List:

  Element   |  Source
--------------------------
*Set Cond Source_On_Surfaces *elems
*loop elems *OnlyInCond
*format "%8i%8i"
*elemsnum  *cond(Source_Number_On_Surfaces)
*end

####################### Temperatures ######################

Conditions List:

  Node    |    Temperature
------------------------------
*Set Cond Fix_Temperature_On_Lines *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Temperature) 
*end
*Set Cond Fix_Temperature_On_Points *nodes
*loop nodes *OnlyInCond
*format "%5i%10.4e"
*NodesNum           *cond(Temperature) 
*end

########################## Normal Fluxs On Points ##########################

Conditions List:

  Nodes   |   Normal Flux
-----------------------------------
*Set Cond Normal_Flux_On_Points *nodes *canrepeat
*loop nodes *OnlyInCond
  *nodesnum   *cond(Flux)
*end

########################## Normal Fluxs On Lines ##########################

Conditions List:

 Element |       Nodes      |   Normal Flux
--------------------------------------------
*Set Cond Normal_Flux_On_Lines *elems *canrepeat
*loop elems *OnlyInCond
*format "%5i%7i%7i"
*elemsnum  *localnodes  *cond(Flux,real)
*end

########################### Convection On Points ###########################

Conditions List:

  Nodes   |   Coeficient  |      Temperature
--------------------------------------------------
*Set Cond Convection_On_Points *nodes *canrepeat
*loop nodes *OnlyInCond
*nodesnum         *cond(Coeficient)          *cond(Temperature)
*end 

########################### Convection On Lines ###########################

Conditions List:

 Element |       Nodes      |   Coeficient  |  Temperature
-------------------------------------------------------------
*Set Cond Convection_On_Lines *elems *canrepeat
*loop elems *OnlyInCond
*elemsnum   *localnodes  *cond(Coeficient)  *cond(Temperature)
*end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
