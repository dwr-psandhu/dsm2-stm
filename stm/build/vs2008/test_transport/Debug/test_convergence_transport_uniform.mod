	  M  L   k820309    �
          11.1        ӗL                                                                                                       
       D:\delta\trunk\stm\test\transport\test_convergence_transport_uniform.f90 TEST_CONVERGENCE_TRANSPORT_UNIFORM                                                
                                                                                      8                                             
         
               �?        1.D0                                             
         
                         0.D0                                             
         
            T4o�A        1.23456789D8                                             
         
               0@        1.6D1                                             
         
                @        2.D0                                             
         
               @        3.D0      @                               	     
                                           
     
          
                         0.D0                                             
         
               �@        25600.D0      @@                                   
         @@                                   
                                                
          
               �?        1.D0       @ @                                   
          @                                     
                                               
   #     @                                                     #VERBOSE          
  @                                      #     @                                                   #CONVERGE_TRANSPORT_UNIFORM%PRESENT    #CONVERGE_TRANSPORT_UNIFORM%NULL    #VERBOSE    #LABEL    #TEST_FLOW    #TEST_DIFFUSE    #TEST_DECAY    #BOUNDARY_REMOTE    #DETAIL_RESULT                                                PRESENT                                            NULL       
  @                                            
  @                                          1       
  @                                   
        
                                      
        
                                      
        
 @                                            
 @                                      #     @    @                                           	   #GAUSSIAN_DATA%SQRT    #BC_DATA     #XLOC "   #CONC #   #NCELL $   #NVAR !   #ORIGIN %   #TIME &   #DX '   #DT (                                               SQRT      D                                             
 (    p      5 � p    r !       5 � p    r !                   
  @                              "     
       
                                 #            
 )     p    5 � p    r $   p      5 � p    r $     5 � p    r !       5 � p    r $     5 � p    r !                   
                                  $             
                                  !             
                                 %     
        
                                 &     
        
                                 '     
        
                                 (     
  #     @    @                            )                	   #BC_DATA *   #XLOC ,   #CONC -   #NCELL .   #NVAR +   #ORIGIN /   #TIME 0   #DX 1   #DT 2        D                                *            
 ,    p      5 � p    r +       5 � p    r +                   
                                 ,     
       
                                 -            
 -     p    5 � p    r .   p      5 � p    r .     5 � p    r +       5 � p    r .     5 � p    r +                   
                                  .             
                                  +             
                                 /     
        
                                 0     
        
                                 1     
        
                                 2     
  #     @                                 3                  #INITIAL_FINAL_SOLUTION_UNIFORM%EXP 4   #INITIAL_FINAL_SOLUTION_UNIFORM%SQRT 5   #FINE_INITIAL_CONC 6   #FINE_SOLUTION_CONC 9   #IC_CENTER :   #IC_PEAK ;   #CONST_VELOCITY <   #DECAY_RATE =   #TOTAL_TIME >   #ORIGIN ?   #DOMAIN_LENGTH @   #NX_BASE 7   #NCONC 8                                          4     EXP                                        5     SQRT      D @                              6            
 &      p    5 � p 
   r 7   p      5 � p 
   r 7     5 � p    r 8       5 � p 
   r 7     5 � p    r 8                  D @                              9            
 '      p    5 � p 
   r 7   p      5 � p 
   r 7     5 � p    r 8       5 � p 
   r 7     5 � p    r 8                   
  @                              :     
        
  @                              ;     
        
                                 <     
        
                                 =     
        
                                 >     
        
  @                              ?     
        
                                 @     
        
  @                               7             
                                  8       #     @                                  A               	   #GAUSSIAN_GRADIENT_DATA%SQRT B   #BC_DATA C   #XLOC E   #CONC F   #NCELL G   #NVAR D   #ORIGIN H   #TIME I   #DX J   #DT K                                          B     SQRT      D                                C            
 *    p      5 � p    r D       5 � p    r D                   
  @                              E     
       
                                 F            
 +     p    5 � p    r G   p      5 � p    r G     5 � p    r D       5 � p    r G     5 � p    r D                   
                                  G             
                                  D             
                                 H     
        
                                 I     
        
                                 J     
        
                                 K     
     �   t      fn#fn      <   J   STM_PRECISION '   L  ]       STM_REAL+STM_PRECISION "   �  `       ONE+STM_PRECISION #   	  `       ZERO+STM_PRECISION (   i  h       LARGEREAL+STM_PRECISION &   �  a       SIXTEEN+STM_PRECISION "   2  `       TWO+STM_PRECISION $   �  `       THREE+STM_PRECISION     �  8       CONST_DISP_COEF    *  `       ORIGIN #   �  d       BASE_DOMAIN_LENGTH    �  8       DOMAIN_LENGTH    &  8       IC_CENTER    ^  `       IC_PEAK    �  8       CONST_VELOCITY #   �  8       DIFFUSE_START_TIME !   .  8       DIFFUSE_END_TIME 0   f  Q       TEST_CONVERGE_TRANSPORT_UNIFORM 8   �  8   a   TEST_CONVERGE_TRANSPORT_UNIFORM%VERBOSE +   �        CONVERGE_TRANSPORT_UNIFORM 3   �  <      CONVERGE_TRANSPORT_UNIFORM%PRESENT <   -  9      CONVERGE_TRANSPORT_UNIFORM%NULL+SOURCE_SINK 3   f  8   a   CONVERGE_TRANSPORT_UNIFORM%VERBOSE 1   �  @   a   CONVERGE_TRANSPORT_UNIFORM%LABEL 5   �  8   a   CONVERGE_TRANSPORT_UNIFORM%TEST_FLOW 8   	  8   a   CONVERGE_TRANSPORT_UNIFORM%TEST_DIFFUSE 6   N	  8   a   CONVERGE_TRANSPORT_UNIFORM%TEST_DECAY ;   �	  8   a   CONVERGE_TRANSPORT_UNIFORM%BOUNDARY_REMOTE 9   �	  8   a   CONVERGE_TRANSPORT_UNIFORM%DETAIL_RESULT    �	  �       GAUSSIAN_DATA #   �
  9      GAUSSIAN_DATA%SQRT &   �
  �   a   GAUSSIAN_DATA%BC_DATA #   {  8   a   GAUSSIAN_DATA%XLOC #   �  �   a   GAUSSIAN_DATA%CONC $   �  8   a   GAUSSIAN_DATA%NCELL #   �  8   a   GAUSSIAN_DATA%NVAR %     8   a   GAUSSIAN_DATA%ORIGIN #   O  8   a   GAUSSIAN_DATA%TIME !   �  8   a   GAUSSIAN_DATA%DX !   �  8   a   GAUSSIAN_DATA%DT -   �  �       EXTRAPOLATE_HI_BOUNDARY_DATA 5   �  �   a   EXTRAPOLATE_HI_BOUNDARY_DATA%BC_DATA 2   +  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%XLOC 2   c  �   a   EXTRAPOLATE_HI_BOUNDARY_DATA%CONC 3   W  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%NCELL 2   �  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%NVAR 4   �  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%ORIGIN 2   �  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%TIME 0   7  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%DX 0   o  8   a   EXTRAPOLATE_HI_BOUNDARY_DATA%DT /   �  K      INITIAL_FINAL_SOLUTION_UNIFORM 3   �  8      INITIAL_FINAL_SOLUTION_UNIFORM%EXP 4   *  9      INITIAL_FINAL_SOLUTION_UNIFORM%SQRT A   c  �   a   INITIAL_FINAL_SOLUTION_UNIFORM%FINE_INITIAL_CONC B   W  �   a   INITIAL_FINAL_SOLUTION_UNIFORM%FINE_SOLUTION_CONC 9   K  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%IC_CENTER 7   �  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%IC_PEAK >   �  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%CONST_VELOCITY :   �  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%DECAY_RATE :   +  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%TOTAL_TIME 6   c  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%ORIGIN =   �  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%DOMAIN_LENGTH 7   �  8   a   INITIAL_FINAL_SOLUTION_UNIFORM%NX_BASE 5     8   a   INITIAL_FINAL_SOLUTION_UNIFORM%NCONC '   C  �       GAUSSIAN_GRADIENT_DATA ,     9      GAUSSIAN_GRADIENT_DATA%SQRT /   =  �   a   GAUSSIAN_GRADIENT_DATA%BC_DATA ,   �  8   a   GAUSSIAN_GRADIENT_DATA%XLOC ,   	  �   a   GAUSSIAN_GRADIENT_DATA%CONC -   �  8   a   GAUSSIAN_GRADIENT_DATA%NCELL ,   5  8   a   GAUSSIAN_GRADIENT_DATA%NVAR .   m  8   a   GAUSSIAN_GRADIENT_DATA%ORIGIN ,   �  8   a   GAUSSIAN_GRADIENT_DATA%TIME *   �  8   a   GAUSSIAN_GRADIENT_DATA%DX *     8   a   GAUSSIAN_GRADIENT_DATA%DT 