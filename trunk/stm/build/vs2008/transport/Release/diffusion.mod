	  �3  U   k820309    �
          11.1        3�L                                                                                                       
       D:\delta\trunk\stm\src\transport\diffusion.f90 DIFFUSION %     @                                                       #USE_DIFFUSION%ASSOCIATED                                                ASSOCIATED #     @                                                     #CONC    #CONC_PREV    #AREA    #AREA_PREV 	   #AREA_LO 
   #AREA_HI    #AREA_LO_PREV    #AREA_HI_PREV    #DISP_COEF_LO    #DISP_COEF_HI    #DISP_COEF_LO_PREV    #DISP_COEF_HI_PREV    #NCELL    #NVAR    #TIME    #THETA_STM    #DT    #DX                                                                       D @                                          
       p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
      p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                              	            
    p      5 � p    r        5 � p    r                   
  @                              
            
    p      5 � p    r        5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                                          
 	     p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
 
     p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
      p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
      p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                    
  @                                            
  @                                            
  @                                   
        
  @                                   
        
  @                                   
        
  @                                   
  #     @                                                    #EXPLICIT_DIFFUSE_OP    #CONC    #AREA_LO    #AREA_HI    #DISP_COEF_LO    #DISP_COEF_HI    #NCELL    #NVAR    #TIME    #DX     #DT !                                                                                                                                  D                                            
       p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
      p    5 � p    r    p      5 � p    r      5 � p    r        5 � p    r      5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                   
  @                                          
    p      5 � p    r        5 � p    r                    
  @                                            
  @                                            
  @                                   
        
  @                                    
        
  @                              !     
  #     @        
                         "                   #RIGHT_HAND_SIDE #   #EXPLICIT_DIFFUSE_OP &   #AREA_PREV '   #AREA_LO_PREV (   #AREA_HI_PREV )   #DISP_COEF_LO_PREV *   #DISP_COEF_HI_PREV +   #CONC_PREV ,   #THETA -   #NCELL $   #TIME .   #NVAR %   #DX /   #DT 0                                                                                                                            D                                #            
 !      p    5 � p 
   r $   p      5 � p 
   r $     5 � p    r %       5 � p 
   r $     5 � p    r %                  
                                 &            
 "     p    5 � p 
   r $   p      5 � p 
   r $     5 � p    r %       5 � p 
   r $     5 � p    r %                  
                                 '            
 #   p      5 � p 
   r $       5 � p 
   r $                  
                                 (            
 %   p      5 � p 
   r $       5 � p 
   r $                  
                                 )            
 &   p      5 � p 
   r $       5 � p 
   r $                  
                                 *            
 '     p    5 � p 
   r $   p      5 � p 
   r $     5 � p    r %       5 � p 
   r $     5 � p    r %                  
                                 +            
 (     p    5 � p 
   r $   p      5 � p 
   r $     5 � p    r %       5 � p 
   r $     5 � p    r %                  
                                 ,            
 $     p    5 � p 
   r $   p      5 � p 
   r $     5 � p    r %       5 � p 
   r $     5 � p    r %                   
                                 -     
        
                                  $             
                                 .     
        
                                  %             
                                 /     
        
                                 0     
  #     @                                 1                   #CENTER_DIAG 2   #UP_DIAG 5   #DOWN_DIAG 6   #AREA 7   #AREA_LO 8   #AREA_HI 9   #DISP_COEF_LO :   #DISP_COEF_HI ;   #THETA_STM <   #NCELL 3   #TIME =   #NVAR 4   #DX >   #DT ?                                                                                                                               D                                2            
 *      p    5 � p 
   r 3   p      5 � p 
   r 3     5 � p    r 4       5 � p 
   r 3     5 � p    r 4                  D                                5            
 +      p    5 � p 
   r 3   p      5 � p 
   r 3     5 � p    r 4       5 � p 
   r 3     5 � p    r 4                  D                                6            
 )      p    5 � p 
   r 3   p      5 � p 
   r 3     5 � p    r 4       5 � p 
   r 3     5 � p    r 4                  
                                 7            
 ,   p      5 � p 
   r 3       5 � p 
   r 3                  
                                 8            
 -   p      5 � p 
   r 3       5 � p 
   r 3                  
                                 9            
 .   p      5 � p 
   r 3       5 � p 
   r 3                  
                                 :            
 /   p      5 � p 
   r 3       5 � p 
   r 3                  
                                 ;            
 0   p      5 � p 
   r 3       5 � p 
   r 3                   
                                 <     
        
                                  3             
                                 =     
        
                                  4             
                                 >     
        
                                 ?     
  #     @                                 @                   #CENTER_DIAG A   #UP_DIAG D   #DOWN_DIAG E   #RIGHT_HAND_SIDE F   #CONC G   #NCELL B   #NVAR C                                                                
                                 A            
 3     p    5 � p    r B   p      5 � p    r B     5 � p    r C       5 � p    r B     5 � p    r C                  
                                 D            
 4     p    5 � p    r B   p      5 � p    r B     5 � p    r C       5 � p    r B     5 � p    r C                  
                                 E            
 2     p    5 � p    r B   p      5 � p    r B     5 � p    r C       5 � p    r B     5 � p    r C                  
                                 F            
 5     p    5 � p    r B   p      5 � p    r B     5 � p    r C       5 � p    r B     5 � p    r C                  D @                              G            
 6      p    5 � p    r B   p      5 � p    r B     5 � p    r C       5 � p    r B     5 � p    r C                   
  @                               B             
                                  C       #     @                                 H                   #DIFFUSIVE_FLUX_LO I   #DIFFUSIVE_FLUX_HI L   #CONC M   #AREA_LO N   #AREA_HI O   #DISP_COEF_LO P   #DISP_COEF_HI Q   #NCELL J   #NVAR K   #TIME R   #DX S   #DT T                                                                                           D @                              I            
       p    5 � p    r J   p      5 � p    r J     5 � p 	   r K       5 � p    r J     5 � p 	   r K                  D @                              L            
       p    5 � p    r J   p      5 � p    r J     5 � p 	   r K       5 � p    r J     5 � p 	   r K                  
  @                              M            
      p    5 � p    r J   p      5 � p    r J     5 � p 	   r K       5 � p    r J     5 � p 	   r K                  
  @                              N            
    p      5 � p    r J       5 � p    r J                  
  @                              O            
    p      5 � p    r J       5 � p    r J                  
  @                              P            
    p      5 � p    r J       5 � p    r J                  
  @                              Q            
     p      5 � p    r J       5 � p    r J                   
  @                               J             
  @                               K             
  @                              R     
        
  @                              S     
        
  @                              T     
     �   A      fn#fn    �   f       USE_DIFFUSION )   C  ?      USE_DIFFUSION%ASSOCIATED    �  �      DIFFUSE      �   a   DIFFUSE%CONC "   �  �   a   DIFFUSE%CONC_PREV    �  �   a   DIFFUSE%AREA "   �  �   a   DIFFUSE%AREA_PREV       �   a   DIFFUSE%AREA_LO     �  �   a   DIFFUSE%AREA_HI %   <  �   a   DIFFUSE%AREA_LO_PREV %   �  �   a   DIFFUSE%AREA_HI_PREV %   d  �   a   DIFFUSE%DISP_COEF_LO %   X	  �   a   DIFFUSE%DISP_COEF_HI *   L
  �   a   DIFFUSE%DISP_COEF_LO_PREV *   @  �   a   DIFFUSE%DISP_COEF_HI_PREV    4  8   a   DIFFUSE%NCELL    l  8   a   DIFFUSE%NVAR    �  8   a   DIFFUSE%TIME "   �  8   a   DIFFUSE%THETA_STM      8   a   DIFFUSE%DT    L  8   a   DIFFUSE%DX ,   �  N      EXPLICIT_DIFFUSION_OPERATOR @   �  �   a   EXPLICIT_DIFFUSION_OPERATOR%EXPLICIT_DIFFUSE_OP 1   �  �   a   EXPLICIT_DIFFUSION_OPERATOR%CONC 4   �  �   a   EXPLICIT_DIFFUSION_OPERATOR%AREA_LO 4   N  �   a   EXPLICIT_DIFFUSION_OPERATOR%AREA_HI 9   �  �   a   EXPLICIT_DIFFUSION_OPERATOR%DISP_COEF_LO 9   v  �   a   EXPLICIT_DIFFUSION_OPERATOR%DISP_COEF_HI 2   
  8   a   EXPLICIT_DIFFUSION_OPERATOR%NCELL 1   B  8   a   EXPLICIT_DIFFUSION_OPERATOR%NVAR 1   z  8   a   EXPLICIT_DIFFUSION_OPERATOR%TIME /   �  8   a   EXPLICIT_DIFFUSION_OPERATOR%DX /   �  8   a   EXPLICIT_DIFFUSION_OPERATOR%DT *   "  �      CONSTRUCT_RIGHT_HAND_SIDE :   �  �   a   CONSTRUCT_RIGHT_HAND_SIDE%RIGHT_HAND_SIDE >   �  �   a   CONSTRUCT_RIGHT_HAND_SIDE%EXPLICIT_DIFFUSE_OP 4   �  �   a   CONSTRUCT_RIGHT_HAND_SIDE%AREA_PREV 7   .  �   a   CONSTRUCT_RIGHT_HAND_SIDE%AREA_LO_PREV 7   �  �   a   CONSTRUCT_RIGHT_HAND_SIDE%AREA_HI_PREV <   V  �   a   CONSTRUCT_RIGHT_HAND_SIDE%DISP_COEF_LO_PREV <   J  �   a   CONSTRUCT_RIGHT_HAND_SIDE%DISP_COEF_HI_PREV 4   >  �   a   CONSTRUCT_RIGHT_HAND_SIDE%CONC_PREV 0   2  8   a   CONSTRUCT_RIGHT_HAND_SIDE%THETA 0   j  8   a   CONSTRUCT_RIGHT_HAND_SIDE%NCELL /   �  8   a   CONSTRUCT_RIGHT_HAND_SIDE%TIME /   �  8   a   CONSTRUCT_RIGHT_HAND_SIDE%NVAR -     8   a   CONSTRUCT_RIGHT_HAND_SIDE%DX -   J  8   a   CONSTRUCT_RIGHT_HAND_SIDE%DT +   �  n      CONSTRUCT_DIFFUSION_MATRIX 7   �  �   a   CONSTRUCT_DIFFUSION_MATRIX%CENTER_DIAG 3   �  �   a   CONSTRUCT_DIFFUSION_MATRIX%UP_DIAG 5   �   �   a   CONSTRUCT_DIFFUSION_MATRIX%DOWN_DIAG 0   �!  �   a   CONSTRUCT_DIFFUSION_MATRIX%AREA 3   `"  �   a   CONSTRUCT_DIFFUSION_MATRIX%AREA_LO 3   �"  �   a   CONSTRUCT_DIFFUSION_MATRIX%AREA_HI 8   �#  �   a   CONSTRUCT_DIFFUSION_MATRIX%DISP_COEF_LO 8   $  �   a   CONSTRUCT_DIFFUSION_MATRIX%DISP_COEF_HI 5   �$  8   a   CONSTRUCT_DIFFUSION_MATRIX%THETA_STM 1   �$  8   a   CONSTRUCT_DIFFUSION_MATRIX%NCELL 0    %  8   a   CONSTRUCT_DIFFUSION_MATRIX%TIME 0   X%  8   a   CONSTRUCT_DIFFUSION_MATRIX%NVAR .   �%  8   a   CONSTRUCT_DIFFUSION_MATRIX%DX .   �%  8   a   CONSTRUCT_DIFFUSION_MATRIX%DT     &  �       SOLVE "   �&  �   a   SOLVE%CENTER_DIAG    �'  �   a   SOLVE%UP_DIAG     �(  �   a   SOLVE%DOWN_DIAG &   �)  �   a   SOLVE%RIGHT_HAND_SIDE    �*  �   a   SOLVE%CONC    �+  8   a   SOLVE%NCELL    �+  8   a   SOLVE%NVAR    ,  <      DIFFUSIVE_FLUX 1   M-  �   a   DIFFUSIVE_FLUX%DIFFUSIVE_FLUX_LO 1   A.  �   a   DIFFUSIVE_FLUX%DIFFUSIVE_FLUX_HI $   5/  �   a   DIFFUSIVE_FLUX%CONC '   )0  �   a   DIFFUSIVE_FLUX%AREA_LO '   �0  �   a   DIFFUSIVE_FLUX%AREA_HI ,   Q1  �   a   DIFFUSIVE_FLUX%DISP_COEF_LO ,   �1  �   a   DIFFUSIVE_FLUX%DISP_COEF_HI %   y2  8   a   DIFFUSIVE_FLUX%NCELL $   �2  8   a   DIFFUSIVE_FLUX%NVAR $   �2  8   a   DIFFUSIVE_FLUX%TIME "   !3  8   a   DIFFUSIVE_FLUX%DX "   Y3  8   a   DIFFUSIVE_FLUX%DT 