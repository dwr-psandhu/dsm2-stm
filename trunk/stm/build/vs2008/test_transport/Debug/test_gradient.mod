	  4+  x   k820309    �
          11.1        ���L                                                                                                       
       D:\delta\trunk\stm\test\transport\test_gradient.f90 TEST_GRADIENT                  @                             
                                                      
                                                      
                                                                                      8#     @                                                    #A_NCELL    #A_NVAR          
                                               
                                              @@                                         
        &           &                                                               	     
         
            T4o�A        1.23456789D8                                        
     
         
                @        2.D0                                             
           
                �                                                     
         
               �?        5.D-1                                             
         
               �?        1.D0                                             
         
                         0.D0#     @                                                     #     @     @                                              #VAR1    #VAR2    #MESSAGE          
                                               
                                               
                                            1 #     @     @                                             #VAR1    #VAR2    #MESSAGE          
                                       
        
                                       
        
                                            1 #     @     @                                              #VAR1    #VAR2    #MESSAGE          
                                       	        
                                       	        
                                            1 #     @     @                                              #VAR1    #VAR2    #MESSAGE          
                                               
                                               
                                            1 #     @     @                                              #ASSERT_EQUALS_STRING_%TRIM !   #VAR1 "   #VAR2 #   #MESSAGE $                                           !     TRIM       
                                 "            1       
                                 #            1       
                                $            1 #     @     @                           %                   #VAR1 &   #VAR2 '   #MESSAGE (         
                                 &             
                                 '             
                                (            1 #     @     @                           )                  #ASSERT_EQUALS_REAL_WITHIN_RANGE_%ABS *   #VAR1 +   #VAR2 ,   #VAR3 -   #MESSAGE .                                           *     ABS       
                                  +     	        
                                  ,     	        
                                  -     	        
                                .            1 #     @     @                           /                  #ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%ABS 0   #VAR1 1   #VAR2 2   #VAR3 3   #MESSAGE 4                                           0     ABS       
                                  1     
        
                                  2     
        
                                  3     
        
                                4            1 #     @     @                           5                   #VAR1 6   #VAR2 7   #N 8   #MESSAGE 9        
                                  6                p      5 O p        5 O p                   
                                  7                p      5 O p        5 O p                    
                                  8             
                                9            1 #     @     @                           :                   #VAR1 ;   #VAR2 <   #N =   #MESSAGE >        
                                  ;            
    p      5 O p        5 O p                   
                                  <            
    p      5 O p        5 O p                    
                                  =             
                                >            1 #     @     @                           ?                   #VAR1 @   #VAR2 A   #N B   #MESSAGE C        
                                  @            	    p      5 O p        5 O p                   
                                  A            	    p      5 O p        5 O p                    
                                  B             
                                C            1 #     @     @                           D                   #VAR1 E   #VAR2 F   #N G   #MESSAGE H   ,     
                                 E                 p      5 O p        5 O p              1 ,     
                                 F                 p      5 O p        5 O p              1       
                                  G             
                                H            1 #     @     @                           I                   #VAR1 J   #VAR2 K   #N L   #MESSAGE M        
                                 J                p      5 O p        5 O p                   
                                 K                p      5 O p        5 O p                    
                                  L             
                                M            1 #     @     @                           N                  #ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%ABS O   #ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%MAXVAL P   #VAR1 Q   #VAR2 R   #N S   #VAR3 T   #MESSAGE U                                           O     ABS                                         P     MAXVAL      
                                  Q            	    p      5 O p        5 O p                   
                                  R            	    p      5 O p        5 O p                    
                                  S             
                                  T     	        
                                U            1 #     @     @                           V                  #ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%ABS W   #ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%MAXVAL X   #VAR1 Y   #VAR2 Z   #N [   #VAR3 \   #MESSAGE ]                                           W     ABS                                         X     MAXVAL      
                                  Y            
    p      5 O p        5 O p                   
                                  Z            
 	   p      5 O p        5 O p                    
                                  [             
                                  \     
        
                                ]            1 #     @     @                           ^                   #VAR1 _   #VAR2 `   #N a   #M b   #MESSAGE c        
                                  _                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                  `                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  a             
                                  b             
                                c            1 #     @     @                           d                   #VAR1 e   #VAR2 f   #N g   #M h   #MESSAGE i        
                                  e            
      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                  f            
      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  g             
                                  h             
                                i            1 #     @     @                           j                   #VAR1 k   #VAR2 l   #N m   #M n   #MESSAGE o        
                                  k            	      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                  l            	      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  m             
                                  n             
                                o            1 #     @     @                           p                   #VAR1 q   #VAR2 r   #N s   #M t   #MESSAGE u        
                                 q                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                 r                  p    5 O p    p      5 O p      5 O p        5 O p      5 O p                    
                                  s             
                                  t             
                                u            1 #     @                                  v                    #     @                                  w                       �   J      fn#fn    �   <   J   FRUIT    "  <   J   STM_PRECISION     ^  <   J   STATE_VARIABLES '   �  ]       STM_REAL+STM_PRECISION /   �  ]       ALLOCATE_STATE+STATE_VARIABLES 7   T  8   a   ALLOCATE_STATE%A_NCELL+STATE_VARIABLES 6   �  8   a   ALLOCATE_STATE%A_NVAR+STATE_VARIABLES %   �  t       CONC+STATE_VARIABLES (   8  h       LARGEREAL+STM_PRECISION "   �  `       TWO+STM_PRECISION $      \       MINUS+STM_PRECISION #   \  a       HALF+STM_PRECISION "   �  `       ONE+STM_PRECISION #     `       ZERO+STM_PRECISION 1   }  D       DEALLOCATE_STATE+STATE_VARIABLES )   �  e      ASSERT_EQUALS_INT_+FRUIT .   &  8   a   ASSERT_EQUALS_INT_%VAR1+FRUIT .   ^  8   a   ASSERT_EQUALS_INT_%VAR2+FRUIT 1   �  @   a   ASSERT_EQUALS_INT_%MESSAGE+FRUIT ,   �  e      ASSERT_EQUALS_DOUBLE_+FRUIT 1   ;  8   a   ASSERT_EQUALS_DOUBLE_%VAR1+FRUIT 1   s  8   a   ASSERT_EQUALS_DOUBLE_%VAR2+FRUIT 4   �  @   a   ASSERT_EQUALS_DOUBLE_%MESSAGE+FRUIT *   �  e      ASSERT_EQUALS_REAL_+FRUIT /   P  8   a   ASSERT_EQUALS_REAL_%VAR1+FRUIT /   �  8   a   ASSERT_EQUALS_REAL_%VAR2+FRUIT 2   �  @   a   ASSERT_EQUALS_REAL_%MESSAGE+FRUIT -    	  e      ASSERT_EQUALS_LOGICAL_+FRUIT 2   e	  8   a   ASSERT_EQUALS_LOGICAL_%VAR1+FRUIT 2   �	  8   a   ASSERT_EQUALS_LOGICAL_%VAR2+FRUIT 5   �	  @   a   ASSERT_EQUALS_LOGICAL_%MESSAGE+FRUIT ,   
  �      ASSERT_EQUALS_STRING_+FRUIT 1   �
  9      ASSERT_EQUALS_STRING_%TRIM+FRUIT 1   �
  @   a   ASSERT_EQUALS_STRING_%VAR1+FRUIT 1     @   a   ASSERT_EQUALS_STRING_%VAR2+FRUIT 4   S  @   a   ASSERT_EQUALS_STRING_%MESSAGE+FRUIT -   �  e      ASSERT_EQUALS_COMPLEX_+FRUIT 2   �  8   a   ASSERT_EQUALS_COMPLEX_%VAR1+FRUIT 2   0  8   a   ASSERT_EQUALS_COMPLEX_%VAR2+FRUIT 5   h  @   a   ASSERT_EQUALS_COMPLEX_%MESSAGE+FRUIT 7   �  �      ASSERT_EQUALS_REAL_WITHIN_RANGE_+FRUIT ;   A  8      ASSERT_EQUALS_REAL_WITHIN_RANGE_%ABS+FRUIT <   y  8   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%VAR1+FRUIT <   �  8   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%VAR2+FRUIT <   �  8   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%VAR3+FRUIT ?   !  @   a   ASSERT_EQUALS_REAL_WITHIN_RANGE_%MESSAGE+FRUIT 9   a  �      ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_+FRUIT =   �  8      ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%ABS+FRUIT >   4  8   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%VAR1+FRUIT >   l  8   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%VAR2+FRUIT >   �  8   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%VAR3+FRUIT A   �  @   a   ASSERT_EQUALS_DOUBLE_WITHIN_RANGE_%MESSAGE+FRUIT ,     l      ASSERT_EQUALS_1D_INT_+FRUIT 1   �  �   a   ASSERT_EQUALS_1D_INT_%VAR1+FRUIT 1     �   a   ASSERT_EQUALS_1D_INT_%VAR2+FRUIT .   �  8   a   ASSERT_EQUALS_1D_INT_%N+FRUIT 4   �  @   a   ASSERT_EQUALS_1D_INT_%MESSAGE+FRUIT /     l      ASSERT_EQUALS_1D_DOUBLE_+FRUIT 4   t  �   a   ASSERT_EQUALS_1D_DOUBLE_%VAR1+FRUIT 4   �  �   a   ASSERT_EQUALS_1D_DOUBLE_%VAR2+FRUIT 1   |  8   a   ASSERT_EQUALS_1D_DOUBLE_%N+FRUIT 7   �  @   a   ASSERT_EQUALS_1D_DOUBLE_%MESSAGE+FRUIT -   �  l      ASSERT_EQUALS_1D_REAL_+FRUIT 2   `  �   a   ASSERT_EQUALS_1D_REAL_%VAR1+FRUIT 2   �  �   a   ASSERT_EQUALS_1D_REAL_%VAR2+FRUIT /   h  8   a   ASSERT_EQUALS_1D_REAL_%N+FRUIT 5   �  @   a   ASSERT_EQUALS_1D_REAL_%MESSAGE+FRUIT /   �  l      ASSERT_EQUALS_1D_STRING_+FRUIT 4   L  �   a   ASSERT_EQUALS_1D_STRING_%VAR1+FRUIT 4   �  �   a   ASSERT_EQUALS_1D_STRING_%VAR2+FRUIT 1   \  8   a   ASSERT_EQUALS_1D_STRING_%N+FRUIT 7   �  @   a   ASSERT_EQUALS_1D_STRING_%MESSAGE+FRUIT 0   �  l      ASSERT_EQUALS_1D_COMPLEX_+FRUIT 5   @  �   a   ASSERT_EQUALS_1D_COMPLEX_%VAR1+FRUIT 5   �  �   a   ASSERT_EQUALS_1D_COMPLEX_%VAR2+FRUIT 2   H  8   a   ASSERT_EQUALS_1D_COMPLEX_%N+FRUIT 8   �  @   a   ASSERT_EQUALS_1D_COMPLEX_%MESSAGE+FRUIT :   �  �      ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_+FRUIT >   �  8      ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%ABS+FRUIT A   �  ;      ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%MAXVAL+FRUIT ?     �   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%VAR1+FRUIT ?   �  �   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%VAR2+FRUIT <     8   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%N+FRUIT ?   F  8   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%VAR3+FRUIT B   ~  @   a   ASSERT_EQUALS_1D_REAL_WITHIN_RANGE_%MESSAGE+FRUIT <   �  �      ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_+FRUIT @   �  8      ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%ABS+FRUIT C   �  ;      ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%MAXVAL+FRUIT A     �   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%VAR1+FRUIT A   �  �   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%VAR2+FRUIT >     8   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%N+FRUIT A   H  8   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%VAR3+FRUIT D   �  @   a   ASSERT_EQUALS_1D_DOUBLE_WITHIN_RANGE_%MESSAGE+FRUIT ,   �  s      ASSERT_EQUALS_2D_INT_+FRUIT 1   3   �   a   ASSERT_EQUALS_2D_INT_%VAR1+FRUIT 1   �   �   a   ASSERT_EQUALS_2D_INT_%VAR2+FRUIT .   �!  8   a   ASSERT_EQUALS_2D_INT_%N+FRUIT .   "  8   a   ASSERT_EQUALS_2D_INT_%M+FRUIT 4   ;"  @   a   ASSERT_EQUALS_2D_INT_%MESSAGE+FRUIT /   {"  s      ASSERT_EQUALS_2D_DOUBLE_+FRUIT 4   �"  �   a   ASSERT_EQUALS_2D_DOUBLE_%VAR1+FRUIT 4   �#  �   a   ASSERT_EQUALS_2D_DOUBLE_%VAR2+FRUIT 1   �$  8   a   ASSERT_EQUALS_2D_DOUBLE_%N+FRUIT 1   �$  8   a   ASSERT_EQUALS_2D_DOUBLE_%M+FRUIT 7   �$  @   a   ASSERT_EQUALS_2D_DOUBLE_%MESSAGE+FRUIT -   6%  s      ASSERT_EQUALS_2D_REAL_+FRUIT 2   �%  �   a   ASSERT_EQUALS_2D_REAL_%VAR1+FRUIT 2   u&  �   a   ASSERT_EQUALS_2D_REAL_%VAR2+FRUIT /   A'  8   a   ASSERT_EQUALS_2D_REAL_%N+FRUIT /   y'  8   a   ASSERT_EQUALS_2D_REAL_%M+FRUIT 5   �'  @   a   ASSERT_EQUALS_2D_REAL_%MESSAGE+FRUIT 0   �'  s      ASSERT_EQUALS_2D_COMPLEX_+FRUIT 5   d(  �   a   ASSERT_EQUALS_2D_COMPLEX_%VAR1+FRUIT 5   0)  �   a   ASSERT_EQUALS_2D_COMPLEX_%VAR2+FRUIT 2   �)  8   a   ASSERT_EQUALS_2D_COMPLEX_%N+FRUIT 2   4*  8   a   ASSERT_EQUALS_2D_COMPLEX_%M+FRUIT 8   l*  @   a   ASSERT_EQUALS_2D_COMPLEX_%MESSAGE+FRUIT #   �*  D       TEST_GRADIENT_CALC    �*  D       TEST_LIMITER 