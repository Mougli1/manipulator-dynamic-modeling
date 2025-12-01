function dT = dT_Matrix ( th , alph , type )

    if type == 1
        dT = [ - sin(th)             , - cos(th)             , 0 , 0 ; 
                 cos(alph) * cos(th) , - cos(alph) * sin(th) , 0 , 0 ;
                 sin(alph) * cos(th) , - sin(alph) * sin(th) , 0 , 0 ;             
                 0                   , 0                     , 0 , 1 ];
    elseif type == 0
        dT = [ 0 , 0 , 0 , 0 ; 
               0 , 0 , 0 , - sin(alph) ;
               0 , 0 , 0 , cos(alph) ;             
               0 , 0 , 0 , 1 ];
    end

end