	  �+  [   k820309    �
          11.1        0I�L                                                                                                       
       D:\delta\trunk\stm\test\transport\test_convergence_transport.f90 TEST_CONVERGENCE_TRANSPORT                                                 
       STM_REAL ZERO                                                                                8                                             
         
                         0.D0       @                                   
   #     @                                                    #TEST_CONVERGENCE%MINVAL    #TEST_CONVERGENCE%ABS    #TEST_CONVERGENCE%MAXVAL    #TEST_CONVERGENCE%TRIM 	   #TEST_CONVERGENCE%DBLE 
   #TEST_CONVERGENCE%ASSOCIATED    #TEST_CONVERGENCE%PRESENT    #LABEL    #HYDRO    #BC_ADVECTIVE_FLUX    #BC_DIFFUSIVE_FLUX %   #BC_DIFFUSIVE_MATRIX 2   #SOURCE_TERM D   #DOMAIN_LENGTH L   #TOTAL_TIME M   #START_TIME N   #FINE_INITIAL_CONC O   #FINE_SOLUTION R   #NSTEP_BASE S   #NX_BASE P   #NCONC Q   #VERBOSE T   #DETAIL_PRINTOUT U                                                                                                                                                               MINVAL                                             ABS                                             MAXVAL                                        	     TRIM                                        
     DBLE                                             ASSOCIATED                                             PRESENT       
  @                                          1 #     @
                        `            	            
   #FLOW    #FLOW_LO    #FLOW_HI    #AREA    #AREA_LO    #AREA_HI    #NCELL    #TIME    #DX    #DT                                                   
     p      5 O p        5 O p                                                             
     p      5 O p        5 O p                                                             
     p      5 O p        5 O p                                                             
     p      5 O p        5 O p                                                             
     p      5 O p        5 O p                                                             
     p      5 O p        5 O p                    
                                              
                                    
         
                                    
         
                                    
   #     @
                        `            	               #FLUX_LO    #FLUX_HI    #CONC_LO    #CONC_HI    #FLOW_LO    #FLOW_HI    #NCELL     #NVAR !   #TIME "   #DT #   #DX $        
                                          
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                          
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                           
 	      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                           
 
      p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                                           
     p      5 O p        5 O p                   
                                           
     p      5 O p        5 O p                    
                                               
                                !              
                               "     
         
                               #     
         
                               $     
   #     @
                         `       %     	               #DIFFUSIVE_FLUX_LO &   #DIFFUSIVE_FLUX_HI '   #CONC (   #AREA_LO )   #AREA_HI *   #DISP_COEF_LO +   #DISP_COEF_HI ,   #NCELL -   #NVAR .   #TIME /   #DX 0   #DT 1        
                              &            
       p    5 O p    p      5 O p      5 O p 	       5 O p      5 O p 	                  
                              '            
       p    5 O p    p      5 O p      5 O p 	       5 O p      5 O p 	                  
                               (            
       p    5 O p    p      5 O p      5 O p 	       5 O p      5 O p 	                  
                               )            
     p      5 O p        5 O p                   
                               *            
     p      5 O p        5 O p                   
                               +            
     p      5 O p        5 O p                   
                               ,            
     p      5 O p        5 O p                    
                                -              
                                .              
                               /     
         
                               0     
         
                               1     
   #     @
                         `       2     	               #CENTER_DIAG 3   #UP_DIAG 4   #DOWN_DIAG 5   #RIGHT_HAND_SIDE 6   #CONC 7   #EXPLICIT_DIFFUSE_OP 8   #AREA 9   #AREA_LO :   #AREA_HI ;   #DISP_COEF_LO <   #DISP_COEF_HI =   #THETA_STM >   #NCELL ?   #TIME @   #NVAR A   #DX B   #DT C        
                              3            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                              4            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                              5            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                              6            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                               7            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                               8            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                               9            
     p      5 O p        5 O p                   
                               :            
     p      5 O p        5 O p                   
                               ;            
     p      5 O p        5 O p                   
                               <            
     p      5 O p        5 O p                   
                               =            
     p      5 O p        5 O p                    
                               >     
         
                                ?              
                               @     
         
                                A              
                               B     
         
                               C     
   #     @
   @                     `       D     	               #SOURCE E   #CONC F   #AREA G   #FLOW H   #NCELL I   #NVAR J   #TIME K        
                              E            
       p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                               F            
        p    5 O p    p      5 O p      5 O p        5 O p      5 O p                   
                               G            
 !    p      5 O p        5 O p                   
                               H            
 "    p      5 O p        5 O p                    
                                I              
                                J              
                               K     
         
                                 L     
        
                                 M     
        
                                 N     
       
  @                              O            
 #     p    5 � p    r P   p      5 � p    r P     5 � p    r Q       5 � p    r P     5 � p    r Q                  
  @                              R            
 $     p    5 � p    r P   p      5 � p    r P     5 � p    r Q       5 � p    r P     5 � p    r Q                   
  @                               S             
  @                               P             
  @                               Q             
                                  T             
 @                               U       #     @                                 V                   #DISP_COEF_LO W   #DISP_COEF_HI Y   #NCELL X   #TIME Z                                                                                                                                    D                                W            
 1    p      5 � p    r X       5 � p    r X                  D                                Y            
 2    p      5 � p    r X       5 � p    r X                                          @         X              
                                 Z     
     �   d      fn#fn       J   J  STM_PRECISION '   J  ]       STM_REAL+STM_PRECISION #   �  `       ZERO+STM_PRECISION !     8       CONST_DISPERSION !   ?  �      TEST_CONVERGENCE (   �  ;      TEST_CONVERGENCE%MINVAL %     8      TEST_CONVERGENCE%ABS (   D  ;      TEST_CONVERGENCE%MAXVAL &     9      TEST_CONVERGENCE%TRIM &   �  9      TEST_CONVERGENCE%DBLE ,   �  ?      TEST_CONVERGENCE%ASSOCIATED )   0  <      TEST_CONVERGENCE%PRESENT '   l  @   a   TEST_CONVERGENCE%LABEL '   �  �      TEST_CONVERGENCE%HYDRO 7   ]  �   a   TEST_CONVERGENCE%HYDRO%FLOW+HYDRO_DATA :   �  �   a   TEST_CONVERGENCE%HYDRO%FLOW_LO+HYDRO_DATA :   e  �   a   TEST_CONVERGENCE%HYDRO%FLOW_HI+HYDRO_DATA 7   �  �   a   TEST_CONVERGENCE%HYDRO%AREA+HYDRO_DATA :   m	  �   a   TEST_CONVERGENCE%HYDRO%AREA_LO+HYDRO_DATA :   �	  �   a   TEST_CONVERGENCE%HYDRO%AREA_HI+HYDRO_DATA 8   u
  8   a   TEST_CONVERGENCE%HYDRO%NCELL+HYDRO_DATA 7   �
  8   a   TEST_CONVERGENCE%HYDRO%TIME+HYDRO_DATA 5   �
  8   a   TEST_CONVERGENCE%HYDRO%DX+HYDRO_DATA 5     8   a   TEST_CONVERGENCE%HYDRO%DT+HYDRO_DATA 3   U  �      TEST_CONVERGENCE%BC_ADVECTIVE_FLUX N     �   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%FLUX_LO+BOUNDARY_ADVECTION N   �  �   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%FLUX_HI+BOUNDARY_ADVECTION N   �  �   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%CONC_LO+BOUNDARY_ADVECTION N   z  �   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%CONC_HI+BOUNDARY_ADVECTION N   F  �   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%FLOW_LO+BOUNDARY_ADVECTION N   �  �   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%FLOW_HI+BOUNDARY_ADVECTION L   N  8   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%NCELL+BOUNDARY_ADVECTION K   �  8   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%NVAR+BOUNDARY_ADVECTION K   �  8   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%TIME+BOUNDARY_ADVECTION I   �  8   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%DT+BOUNDARY_ADVECTION I   .  8   a   TEST_CONVERGENCE%BC_ADVECTIVE_FLUX%DX+BOUNDARY_ADVECTION 3   f  �      TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX X   O  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%DIFFUSIVE_FLUX_LO+BOUNDARY_DIFFUSION X     �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%DIFFUSIVE_FLUX_HI+BOUNDARY_DIFFUSION K   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%CONC+BOUNDARY_DIFFUSION N   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%AREA_LO+BOUNDARY_DIFFUSION N   7  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%AREA_HI+BOUNDARY_DIFFUSION S   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%DISP_COEF_LO+BOUNDARY_DIFFUSION S   ?  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%DISP_COEF_HI+BOUNDARY_DIFFUSION L   �  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%NCELL+BOUNDARY_DIFFUSION K   �  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%NVAR+BOUNDARY_DIFFUSION K   3  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%TIME+BOUNDARY_DIFFUSION I   k  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%DX+BOUNDARY_DIFFUSION I   �  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_FLUX%DT+BOUNDARY_DIFFUSION 5   �  /     TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX T   
  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%CENTER_DIAG+BOUNDARY_DIFFUSION P   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%UP_DIAG+BOUNDARY_DIFFUSION R   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%DOWN_DIAG+BOUNDARY_DIFFUSION X   n  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%RIGHT_HAND_SIDE+BOUNDARY_DIFFUSION M   :  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%CONC+BOUNDARY_DIFFUSION \     �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%EXPLICIT_DIFFUSE_OP+BOUNDARY_DIFFUSION M   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%AREA+BOUNDARY_DIFFUSION P   V  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%AREA_LO+BOUNDARY_DIFFUSION P   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%AREA_HI+BOUNDARY_DIFFUSION U   ^  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%DISP_COEF_LO+BOUNDARY_DIFFUSION U   �  �   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%DISP_COEF_HI+BOUNDARY_DIFFUSION R   f   8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%THETA_STM+BOUNDARY_DIFFUSION N   �   8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%NCELL+BOUNDARY_DIFFUSION M   �   8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%TIME+BOUNDARY_DIFFUSION M   !  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%NVAR+BOUNDARY_DIFFUSION K   F!  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%DX+BOUNDARY_DIFFUSION K   ~!  8   a   TEST_CONVERGENCE%BC_DIFFUSIVE_MATRIX%DT+BOUNDARY_DIFFUSION -   �!  �      TEST_CONVERGENCE%SOURCE_TERM @   C"  �   a   TEST_CONVERGENCE%SOURCE_TERM%SOURCE+SOURCE_SINK >   #  �   a   TEST_CONVERGENCE%SOURCE_TERM%CONC+SOURCE_SINK >   �#  �   a   TEST_CONVERGENCE%SOURCE_TERM%AREA+SOURCE_SINK >   _$  �   a   TEST_CONVERGENCE%SOURCE_TERM%FLOW+SOURCE_SINK ?   �$  8   a   TEST_CONVERGENCE%SOURCE_TERM%NCELL+SOURCE_SINK >   %  8   a   TEST_CONVERGENCE%SOURCE_TERM%NVAR+SOURCE_SINK >   S%  8   a   TEST_CONVERGENCE%SOURCE_TERM%TIME+SOURCE_SINK /   �%  8   a   TEST_CONVERGENCE%DOMAIN_LENGTH ,   �%  8   a   TEST_CONVERGENCE%TOTAL_TIME ,   �%  8   a   TEST_CONVERGENCE%START_TIME 3   3&  �   a   TEST_CONVERGENCE%FINE_INITIAL_CONC /   ''  �   a   TEST_CONVERGENCE%FINE_SOLUTION ,   (  8   a   TEST_CONVERGENCE%NSTEP_BASE )   S(  8   a   TEST_CONVERGENCE%NX_BASE '   �(  8   a   TEST_CONVERGENCE%NCONC )   �(  8   a   TEST_CONVERGENCE%VERBOSE 1   �(  8   a   TEST_CONVERGENCE%DETAIL_PRINTOUT $   3)  �       GET_DISPERSION_COEF 1   ,*  �   a   GET_DISPERSION_COEF%DISP_COEF_LO 1   �*  �   a   GET_DISPERSION_COEF%DISP_COEF_HI *   T+  8   a   GET_DISPERSION_COEF%NCELL )   �+  8   a   GET_DISPERSION_COEF%TIME 