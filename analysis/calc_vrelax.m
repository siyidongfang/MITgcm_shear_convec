
omega = 2*pi/43200;
Hmax = -zz(end)+dz(end)/2;

if(useTanhShear)
  vrelax = h_shear*Shear *(1+ tanh( (zz    +Hmax/2) / (h_shear/2) )) /2;
end

 if(useLinearShear)

      Nshear_smooth_half = 100;
      for i=1:Nr
           if((zz(i)+Hmax)<h_shear) 
               shearProfile(i)=(zz(i)+Hmax)/h_shear;
           else
               shearProfile(i)=1.;
           end 
       end 
    
      vrelax2 = Shear*h_shear*shearProfile; %%% un-smoothed shear
    
      %--- smooth the velocity shear
      for kLev = 1:Nr
          if(kLev>Nshear_smooth_half) 
              if((Nr-kLev)>=Nshear_smooth_half) 
                  NsmoothStart = kLev-Nshear_smooth_half;
                  NsmoothEnd = kLev+Nshear_smooth_half;
              end 
          end 
    
          if((Nr-kLev)<Nshear_smooth_half) 
              NsmoothStart = kLev-(Nr-kLev);
              NsmoothEnd = Nr;
          end 
    
          if(kLev<=Nshear_smooth_half) 
              NsmoothStart = 1;
              NsmoothEnd = kLev+(kLev-1);
          end 
    
          shearRatio(kLev) = 0.;
          Ndivide(kLev) = 0.;
    
          if(NsmoothEnd>NsmoothStart) 
              for i= NsmoothStart:NsmoothEnd
                   shearRatio(kLev) = shearRatio(kLev) + shearProfile(i);
                   Ndivide(kLev) = Ndivide(kLev) + 1;
              end 
              shearRatio(kLev) = shearRatio(kLev)/Ndivide(kLev);
          else
              shearRatio(kLev) = shearProfile(kLev);
          end 
      end
    
      vrelax = Shear*h_shear*shearRatio;

  end

  