	  �-  z   k820309    �
          11.1        ���L                                                                                                       
       D:\delta\trunk\stm\test\transport\test_explicit_diffusion_operator.f90 TEST_EXPLICIT_DIFFUSION_OPERATOR                                                
                        @                             
                                                      
                                                                                      8                                             
         
            T4o�A        1.23456789D8#     @                                                    #EXPLICIT_DIFFUSE_OP    #CONC    #AREA_LO 	   #AREA_HI 
   #DISP_COEF_LO    #DISP_COEF_HI    #NCELL    #NVAR    #TIME    #DX    #DT                                                     
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                             
      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                 	            
    p      5 O p        5 O p                   
                                 
            
    p      5 O p        5 O p                   
                                             
    p      5 O p        5 O p                   
                                             
    p      5 O p        5 O p                    
                                               
                                               
                                      
        
                                      
        
                                      
                                               
         
         ��&�.>        1.D-9#     @     @                                              #VAR1    #VAR2    #MESSAGE          
                                               
                                               
                                            1 #     @     @                                              #VAR1    #VAR2    #MESSAGE          
                                       
        
                                       
        
                                            1 #     @     @                                              #VAR1    #VAR2    #MESSAGE          
                                       	        
                                       	        
                                            1 #     @     @                                              #VAR1     #VAR2 !   #MESSAGE "         
                                                
                                  !             
                                "            1 #     @     @                           #                  #ASSERT_EQUALS_STRING_%TRIM $   #VAR1 %   #VAR2 &   #MESSAGE '                                           $     TRIM       
                                 %            1       
                                 &            1       
                                '            1 #     @     @                           (                   #VAR1 )   #VAR2 *   #MESSAGE +         
                                 )             
                                 *             
                                +            1 #     @     @                           ,                  #ASSERT_EQUALS_REAL_WITHIN_RANGE_%ABS -   #VAR1 .   #VAR2 /   #VAR3 0   #MESSAGE 1                                           -     ABS       
                                  .     	        
                                  /     	        
                                  0     	        
                                1            1 #     @     @                          2                  #ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%ABS 3   #VAR1 4   #VAR2 5   #VAR3 6   #MESSAGE 7                                           3     ABS       
                                  4     
        
                                  5     
        
                                  6     
        
                                7            1 #     @     @                           8                   #VAR1 9   #VAR2 :   #N ;   #MESSAGE <        
                                  9                p      5 O p        5 O p                   
                                  :                p      5 O p        5 O p                    
                                  ;             
                                <            1 #     @     @                           =                   #VAR1 >   #VAR2 ?   #N @   #MESSAGE A        
                                  >            
    p      5 O p        5 O p                   
                                  ?            
    p      5 O p        5 O p                    
                                  @             
                                A            1 #     @     @                           B                   #VAR1 C   #VAR2 D   #N E   #MESSAGE F        
                                  C            	    p      5 O p        5 O p                   
                                  D            	    p      5 O p        5 O p                    
                                  E             
                                F            1 #     @     @                           G                   #VAR1 H   #VAR2 I   #N J   #MESSAGE K   ,     
                                 H                 p      5 O p        5 O p              1 ,     
                                 I                 p      5 O p        5 O p              1       
                                  J             
                                K            1 #     @     @                           L                   #VAR1 M   #VAR2 N   #N O   #MESSAGE P        
                                 M                p      5 O p        5 O p                   
                                 N                p      5 O p        5 O p                    
                                  O             
                                P            1 #     @     @                           Q                  #ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%ABS R   #ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%MAXVAL S   #VAR1 T   #VAR2 U   #N V   #VAR3 W   #MESSAGE X                                           R     ABS                                         S     MAXVAL      
                                  T            	    p      5 O p        5 O p                   
                                  U            	    p      5 O p        5 O p                    
                                  V             
                                  W     	        
                                X            1 #     @     @                           Y                  #ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%ABS Z   #ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%MAXVAL [   #VAR1 \   #VAR2 ]   #N ^   #VAR3 _   #MESSAGE `                                           Z     ABS                                         [     MAXVAL      
                                  \            
    p      5 O p        5 O p                   
                                  ]            
 	   p      5 O p        5 O p                    
                                  ^             
                                  _     
        
                                `            1 #     @     @                           a                   #VAR1 b   #VAR2 c   #N d   #M e   #MESSAGE f        
                                  b                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                  c                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  d             
                                  e             
                                f            1 #     @     @                           g                   #VAR1 h   #VAR2 i   #N j   #M k   #MESSAGE l        
                                  h            
      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                  i            
      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  j             
                                  k             
                                l            1 #     @     @                           m                   #VAR1 n   #VAR2 o   #N p   #M q   #MESSAGE r        
                                  n            	      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                  o            	      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  p             
                                  q             
                                r            1 #     @     @                           s                   #VAR1 t   #VAR2 u   #N v   #M w   #MESSAGE x        
                                 t                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                 u                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  v             
                                  w             
                                x            1 #     @                                  y                       �   p      fn#fn      <   J   DIFFUSION    H  <   J   FRUIT    �  <   J   STM_PRECISION '   �  ]       STM_REAL+STM_PRECISION (     h       LARGEREAL+STM_PRECISION 6   �  �       EXPLICIT_DIFFUSION_OPERATOR+DIFFUSION J   Y  �   a   EXPLICIT_DIFFUSION_OPERATOR%EXPLICIT_DIFFUSE_OP+DIFFUSION ;   %  �   a   EXPLICIT_DIFFUSION_OPERATOR%CONC+DIFFUSION >   �  �   a   EXPLICIT_DIFFUSION_OPERATOR%AREA_LO+DIFFUSION >   u  �   a   EXPLICIT_DIFFUSION_OPERATOR%AREA_HI+DIFFUSION C   �  �   a   EXPLICIT_DIFFUSION_OPERATOR%DISP_COEF_LO+DIFFUSION C   }  �   a   EXPLICIT_DIFFUSION_OPERATOR%DISP_COEF_HI+DIFFUSION <     8   a   EXPLICIT_DIFFUSION_OPERATOR%NCELL+DIFFUSION ;   9  8   a   EXPLICIT_DIFFUSION_OPERATOR%NVAR+DIFFUSION ;   q  8   a   EXPLICIT_DIFFUSION_OPERATOR%TIME+DIFFUSION 9   �  8   a   EXPLICIT_DIFFUSION_OPERATOR%DX+DIFFUSION 9   �  8   a   EXPLICIT_DIFFUSION_OPERATOR%DT+DIFFUSION '     a       WEAK_EPS+STM_PRECISION )   z  e      ASSERT_EQUALS_INT_+FRUIT .   �  8   a   ASSERT_EQUALS_INT_%VAR1+FRUIT .   	  8   a   ASSERT_EQUALS_INT_%VAR2+FRUIT 1   O	  @   a   ASSERT_EQUALS_INT_%MESSAGE+FRUIT ,   �	  e      ASSERT_EQUALS_DOUBLE_+FRUIT 1   �	  8   a   ASSERT_EQUALS_DOUBLE_%VAR1+FRUIT 1   ,
  8   a   ASSERT_EQUALS_DOUBLE_%VAR2+FRUIT 4   d
  @   a   ASSERT_EQUALS_DOUBLE_%MESSAGE+FRUIT *   �
  e      ASSERT_EQUALS_REAL_+FRUIT /   	  8   a   ASSERT_EQUALS_REAL_%VAR1+FRUIT /   A  8   a   ASSERT_EQUALS_REAL_%VAR2+FRUIT 2   y  @   a   ASSERT_EQUALS_REAL_%MESSAGE+FRUIT -   �  e      ASSERT_EQUALS_LOGICAL_+FRUIT 2     8   a   ASSERT_EQUALS_LOGICAL_%VAR1+FRUIT 2   V  8   a   ASSERT_EQUALS_LOGICAL_%VAR2+FRUIT 5   �  @   a   ASSERT_EQUALS_LOGICAL_%MESSAGE+FRUIT ,   �  �      ASSERT_EQUALS_STRING_+FRUIT 1   S  9      ASSERT_EQUALS_STRING_%TRIM+FRUIT 1   �  @   a   ASSERT_EQUALS_STRING_%VAR1+FRUIT 1   �  @   a   ASSERT_EQUALS_STRING_%VAR2+FRUIT 4     @   a   ASSERT_EQUALS_STRING_%MESSAGE+FRUIT -   L  e      ASSERT_EQUALS_COMPLEX_+FRUIT 2   �  8   a   ASSERT_EQUALS_COMPLEX_%VAR1+FRUIT 2   �  8   a   ASSERT_EQUALS_COMPLEX_%VAR2+FRUIT 5   !  @   a   ASSERT_EQUALS_COMPLEX_%MESSAGE+FRUIT 7   a  �      ASSERT_EQUALS_REAL_WITHIN_RANGE_+FRUIT ;   �  8      ASSERT_EQUALS_REAL_WITHIN_RANGE_%ABS+FRUIT <   2  8   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%VAR1+FRUIT <   j  8   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%VAR2+FRUIT <   �  8   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%VAR3+FRUIT ?   �  @   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%MESSAGE+FRUIT 9     �      ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_+FRUIT =   �  8      ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%ABS+FRUIT >   �  8   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%VAR1+FRUIT >   %  8   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%VAR2+FRUIT >   ]  8   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%VAR3+FRUIT A   �  @   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%MESSAGE+FRUIT ,   �  l      ASSERT_EQUALS_1D_INT_+FRUIT 1   A  �   a   ASSERT_EQUALS_1D_INT_%VAR1+FRUIT 1   �  �   a   ASSERT_EQUALS_1D_INT_%VAR2+FRUIT .   I  8   a   ASSERT_EQUALS_1D_INT_%N+FRUIT 4   �  @   a   ASSERT_EQUALS_1D_INT_%MESSAGE+FRUIT /   �  l      ASSERT_EQUALS_1D_DOUBLE_+FRUIT 4   -  �   a   ASSERT_EQUALS_1D_DOUBLE_%VAR1+FRUIT 4   �  �   a   ASSERT_EQUALS_1D_DOUBLE_%VAR2+FRUIT 1   5  8   a   ASSERT_EQUALS_1D_DOUBLE_%N+FRUIT 7   m  @   a   ASSERT_EQUALS_1D_DOUBLE_%MESSAGE+FRUIT -   �  l      ASSERT_EQUALS_1D_REAL_+FRUIT 2     �   a   ASSERT_EQUALS_1D_REAL_%VAR1+FRUIT 2   �  �   a   ASSERT_EQUALS_1D_REAL_%VAR2+FRUIT /   !  8   a   ASSERT_EQUALS_1D_REAL_%N+FRUIT 5   Y  @   a   ASSERT_EQUALS_1D_REAL_%MESSAGE+FRUIT /   �  l      ASSERT_EQUALS_1D_STRING_+FRUIT 4     �   a   ASSERT_EQUALS_1D_STRING_%VAR1+FRUIT 4   �  �   a   ASSERT_EQUALS_1D_STRING_%VAR2+FRUIT 1     8   a   ASSERT_EQUALS_1D_STRING_%N+FRUIT 7   M  @   a   ASSERT_EQUALS_1D_STRING_%MESSAGE+FRUIT 0   �  l      ASSERT_EQUALS_1D_COMPLEX_+FRUIT 5   �  �   a   ASSERT_EQUALS_1D_COMPLEX_%VAR1+FRUIT 5   }  �   a   ASSERT_EQUALS_1D_COMPLEX_%VAR2+FRUIT 2     8   a   ASSERT_EQUALS_1D_COMPLEX_%N+FRUIT 8   9  @   a   ASSERT_EQUALS_1D_COMPLEX_%MESSAGE+FRUIT :   y  �      ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_+FRUIT >   L  8      ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%ABS+FRUIT A   �  ;      ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%MAXVAL+FRUIT ?   �  �   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%VAR1+FRUIT ?   C  �   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%VAR2+FRUIT <   �  8   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%N+FRUIT ?   �  8   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%VAR3+FRUIT B   7  @   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%MESSAGE+FRUIT <   w  �      ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_+FRUIT @   N   8      ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%ABS+FRUIT C   �   ;      ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%MAXVAL+FRUIT A   �   �   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%VAR1+FRUIT A   E!  �   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%VAR2+FRUIT >   �!  8   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%N+FRUIT A   "  8   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%VAR3+FRUIT D   9"  @   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%MESSAGE+FRUIT ,   y"  s      ASSERT_EQUALS_2D_INT_+FRUIT 1   �"  �   a   ASSERT_EQUALS_2D_INT_%VAR1+FRUIT 1   �#  �   a   ASSERT_EQUALS_2D_INT_%VAR2+FRUIT .   �$  8   a   ASSERT_EQUALS_2D_INT_%N+FRUIT .   �$  8   a   ASSERT_EQUALS_2D_INT_%M+FRUIT 4   �$  @   a   ASSERT_EQUALS_2D_INT_%MESSAGE+FRUIT /   4%  s      ASSERT_EQUALS_2D_DOUBLE_+FRUIT 4   �%  �   a   ASSERT_EQUALS_2D_DOUBLE_%VAR1+FRUIT 4   s&  �   a   ASSERT_EQUALS_2D_DOUBLE_%VAR2+FRUIT 1   ?'  8   a   ASSERT_EQUALS_2D_DOUBLE_%N+FRUIT 1   w'  8   a   ASSERT_EQUALS_2D_DOUBLE_%M+FRUIT 7   �'  @   a   ASSERT_EQUALS_2D_DOUBLE_%MESSAGE+FRUIT -   �'  s      ASSERT_EQUALS_2D_REAL_+FRUIT 2   b(  �   a   ASSERT_EQUALS_2D_REAL_%VAR1+FRUIT 2   .)  �   a   ASSERT_EQUALS_2D_REAL_%VAR2+FRUIT /   �)  8   a   ASSERT_EQUALS_2D_REAL_%N+FRUIT /   2*  8   a   ASSERT_EQUALS_2D_REAL_%M+FRUIT 5   j*  @   a   ASSERT_EQUALS_2D_REAL_%MESSAGE+FRUIT 0   �*  s      ASSERT_EQUALS_2D_COMPLEX_+FRUIT 5   +  �   a   ASSERT_EQUALS_2D_COMPLEX_%VAR1+FRUIT 5   �+  �   a   ASSERT_EQUALS_2D_COMPLEX_%VAR2+FRUIT 2   �,  8   a   ASSERT_EQUALS_2D_COMPLEX_%N+FRUIT 2   �,  8   a   ASSERT_EQUALS_2D_COMPLEX_%M+FRUIT 8   %-  @   a   ASSERT_EQUALS_2D_COMPLEX_%MESSAGE+FRUIT 4   e-  D       TEST_EXPLICIT_INTERIOR_DIFFUSION_OP 