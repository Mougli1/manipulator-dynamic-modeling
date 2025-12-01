function T = T_Matrix(th,alph,d,r)

T = [ cos(th)             , - sin(th)           , 0          , d                  ; 
      cos(alph) * sin(th) , cos(alph) * cos(th) , -sin(alph) , - r * sin(alph) ;
      sin(alph) * sin(th) , sin(alph) * cos(th) , cos(alph)  , r * cos(alph)   ;             
      0                   , 0                   , 0          , 1                 ];

end