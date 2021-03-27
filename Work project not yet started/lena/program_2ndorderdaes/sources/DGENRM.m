function  NRM=DGENRM(NEQ, V, RWT)

%     Compute the weighted root-mean-square vector norm for DGELDA [2].

      NRM = 0.0;
      VMAX = 0.0;
      for i = 1:NEQ
        if(abs(V(i)*RWT(i)) > VMAX) 
            VMAX = abs(V(i)*RWT(i));
        end
      end
      if(VMAX <= 0.0) 
          return;
      end
      SUM = 0.0;
      for i = 1:NEQ
         SUM = SUM + ((V(i)*RWT(i))/VMAX)^2;
      end
      NRM = VMAX*sqrt(SUM/NEQ);


      
      
      
      