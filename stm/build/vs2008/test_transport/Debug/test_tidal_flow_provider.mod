	  ^  7   k820309    ³
          11.1        \uíL                                                                                                       
       D:\delta\trunk\stm\test\transport\test_tidal_hydro.f90 TEST_TIDAL_FLOW_PROVIDER                                                
                                                      
                                                      
                        @                             
                                                      
                                                                                      8#     @                                      	            
   #FLOW    #FLOW_LO 	   #FLOW_HI 
   #AREA    #AREA_LO    #AREA_HI    #NCELL    #TIME    #DX    #DT                                                    
     p      5 O p        5 O p                                                  	            
     p      5 O p        5 O p                                                  
            
     p      5 O p        5 O p                                                              
     p      5 O p        5 O p                                                              
     p      5 O p        5 O p                                                              
     p      5 O p        5 O p                    
                                              
                                     
        
                                     
        
                                     
  #     @                                                
   #TIDAL_FLOW_MODIFIED%DBLE    #TIDAL_FLOW_MODIFIED%SIN    #TIDAL_FLOW_MODIFIED%SQRT    #TIDAL_FLOW_MODIFIED%COS    #FLOW    #FLOW_LO    #FLOW_HI    #AREA    #AREA_LO    #AREA_HI    #NCELL    #TIME    #DX    #DT                                                  DBLE                                              SIN                                              SQRT                                              COS                                                  
     p      5 O p        5 O p                                                               
     p      5 O p        5 O p                                                               
     p      5 O p        5 O p                                                               
     p      5 O p        5 O p                                                               
     p      5 O p        5 O p                                                               
     p      5 O p        5 O p                    
                                               
                                      
        
                                      
        
                                       
                                          !     
          
                         0.D0                                        "     
           
               Ìå@                                                #     
           
                ù@                                                $     
         
                         0.D0#     @                                 %                  #ERROR_NORM%ABS &   #ERROR_NORM%DBLE '   #ERROR_NORM%SQRT (   #NORM_1 )   #NORM_2 *   #NORM_INF +   #WHICH_CELL ,   #VALS -   #REFERENCE .   #NCELL /   #DX 0                                           &     ABS                                         '     DBLE                                         (     SQRT                                       )     
                                         *     
                                         +     
                                          ,             
                                 -            
    p      5 O p        5 O p                   
                                 .            
    p      5 O p        5 O p                    
                                  /                                             0     
                                           1     
         
         Ö&è.>        1.D-9#     @     @                          2                   #VAR1 3   #MESSAGE 4         
                                  3             
                                4            1 #     @                                  5                   #TEST_TIDAL_HYDRO_PROVIDER%MAXVAL 6                                          6     MAXVAL       X      fn#fn    ô   <   J   STM_PRECISION 0   0  <   J   TEST_ADVECTION_TIDAL_EXPERIENCE    l  <   J   HYDRO_DATA    ¨  <   J   FRUIT    ä  <   J   ERROR_METRIC '      ]       STM_REAL+STM_PRECISION )   }  ±       HYDRO_DATA_IF+HYDRO_DATA .   .     a   HYDRO_DATA_IF%FLOW+HYDRO_DATA 1   ²     a   HYDRO_DATA_IF%FLOW_LO+HYDRO_DATA 1   6     a   HYDRO_DATA_IF%FLOW_HI+HYDRO_DATA .   º     a   HYDRO_DATA_IF%AREA+HYDRO_DATA 1   >     a   HYDRO_DATA_IF%AREA_LO+HYDRO_DATA 1   Â     a   HYDRO_DATA_IF%AREA_HI+HYDRO_DATA /   F  8   a   HYDRO_DATA_IF%NCELL+HYDRO_DATA .   ~  8   a   HYDRO_DATA_IF%TIME+HYDRO_DATA ,   ¶  8   a   HYDRO_DATA_IF%DX+HYDRO_DATA ,   î  8   a   HYDRO_DATA_IF%DT+HYDRO_DATA D   &  '      TIDAL_FLOW_MODIFIED+TEST_ADVECTION_TIDAL_EXPERIENCE I   M  9      TIDAL_FLOW_MODIFIED%DBLE+TEST_ADVECTION_TIDAL_EXPERIENCE H     8      TIDAL_FLOW_MODIFIED%SIN+TEST_ADVECTION_TIDAL_EXPERIENCE I   ¾  9      TIDAL_FLOW_MODIFIED%SQRT+TEST_ADVECTION_TIDAL_EXPERIENCE H   ÷  8      TIDAL_FLOW_MODIFIED%COS+TEST_ADVECTION_TIDAL_EXPERIENCE I   /	     a   TIDAL_FLOW_MODIFIED%FLOW+TEST_ADVECTION_TIDAL_EXPERIENCE L   ³	     a   TIDAL_FLOW_MODIFIED%FLOW_LO+TEST_ADVECTION_TIDAL_EXPERIENCE L   7
     a   TIDAL_FLOW_MODIFIED%FLOW_HI+TEST_ADVECTION_TIDAL_EXPERIENCE I   »
     a   TIDAL_FLOW_MODIFIED%AREA+TEST_ADVECTION_TIDAL_EXPERIENCE L   ?     a   TIDAL_FLOW_MODIFIED%AREA_LO+TEST_ADVECTION_TIDAL_EXPERIENCE L   Ã     a   TIDAL_FLOW_MODIFIED%AREA_HI+TEST_ADVECTION_TIDAL_EXPERIENCE J   G  8   a   TIDAL_FLOW_MODIFIED%NCELL+TEST_ADVECTION_TIDAL_EXPERIENCE I     8   a   TIDAL_FLOW_MODIFIED%TIME+TEST_ADVECTION_TIDAL_EXPERIENCE G   ·  8   a   TIDAL_FLOW_MODIFIED%DX+TEST_ADVECTION_TIDAL_EXPERIENCE G   ï  8   a   TIDAL_FLOW_MODIFIED%DT+TEST_ADVECTION_TIDAL_EXPERIENCE ;   '  `       START_TIME+TEST_ADVECTION_TIDAL_EXPERIENCE ;     \       TOTAL_TIME+TEST_ADVECTION_TIDAL_EXPERIENCE >   ã  \       DOMAIN_LENGTH+TEST_ADVECTION_TIDAL_EXPERIENCE #   ?  `       ZERO+STM_PRECISION (     ä       ERROR_NORM+ERROR_METRIC ,     8      ERROR_NORM%ABS+ERROR_METRIC -   »  9      ERROR_NORM%DBLE+ERROR_METRIC -   ô  9      ERROR_NORM%SQRT+ERROR_METRIC /   -  8   a   ERROR_NORM%NORM_1+ERROR_METRIC /   e  8   a   ERROR_NORM%NORM_2+ERROR_METRIC 1     8   a   ERROR_NORM%NORM_INF+ERROR_METRIC 3   Õ  8   a   ERROR_NORM%WHICH_CELL+ERROR_METRIC -        a   ERROR_NORM%VALS+ERROR_METRIC 2        a   ERROR_NORM%REFERENCE+ERROR_METRIC .     8   a   ERROR_NORM%NCELL+ERROR_METRIC +   M  8   a   ERROR_NORM%DX+ERROR_METRIC '     a       WEAK_EPS+STM_PRECISION +   æ  [      ASSERT_TRUE_LOGICAL_+FRUIT 0   A  8   a   ASSERT_TRUE_LOGICAL_%VAR1+FRUIT 3   y  @   a   ASSERT_TRUE_LOGICAL_%MESSAGE+FRUIT *   ¹  j       TEST_TIDAL_HYDRO_PROVIDER 1   #  ;      TEST_TIDAL_HYDRO_PROVIDER%MAXVAL 