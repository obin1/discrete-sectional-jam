###################################################################################
# This function outputs the collision rate frequency function
# Based on a MATLAB tutorial, Li and Cai 2020: https://doi.org/10.1016/j.jaerosci.2020.105615
# "fuchs": Fuchs transition regime limiting sphere collision theory
# i & j are the number of monomers in the two colliding particles

# Note: all other collision kernels besides "fuchs" are deprecated from the MATLAB version
####################################################################################
function makecollisionkernel(kerneltype,T,P,rho,mMono)

    boltz = 1.38e-23
    v_mono = mMono/rho

    function c_kernel_fuchs(i,j)

        # reference temperature; mean free path; pressure & viscosity
        T0 = 300
        mfp0 = 67E-9
        P0 = 1.01E5
        mu0 = 1.82E-5
        
        mu = mu0*(T0+110.4)/(T+110.4)*(T/T0)^(1.5);  # corrected viscosity
        mfp = mfp0*(T/T0)*(P0/P)*(1+110.4/T0)/(1+110.4/T); # corrected mean free path()
        
        
        rpi = ((i.*mMono)/rho*6/pi).^(1/3)/2
        rpj = ((j.*mMono)/rho*6/pi).^(1/3)/2; # size in m
        
        kni = mfp./rpi
        knj = mfp./rpj
        
        a = 1.257; # Cunningham correction coefficients
        b = 0.4
        c = 1.1
        
        Cc_i = 1 + kni.*(a + b*exp(-c./kni))
        Cc_j = 1 + knj.*(a + b*exp(-c./knj))
        
        Bi = Cc_i./(6*pi*mu.*rpi); #s/kg
        Bj = Cc_j./(6*pi*mu.*rpj)
        Di = Bi.*boltz*T; # m^2/s; close to Table 2.1 in Friedlander's book
        Dj = Bj.*boltz*T
        
        mi = mMono*i;  # mass
        mj = mMono*j;  # mass
        
        ci = (8*boltz*T./(pi*mi)).^0.5; # m/s, mean thermal velocity of the particles
        cj = (8*boltz*T./(pi*mj)).^0.5
        
        # mean traveling time
        tau_i = (i.*mMono.*Cc_i)/6/pi./mu./rpi
        tau_j = (j.*mMono.*Cc_j)/6/pi./mu./rpj
        
        # "mean free path" of the particles
        lbi = tau_i.*ci
        lbj = tau_j.*cj
        
        # delta values
        deltai = 1/6/rpi./lbi.*((2*rpi + lbi).^3 -(4*rpi.^2+lbi.^2).^1.5)-2*rpi
        deltaj = 1/6/rpj./lbj.*((2*rpj + lbj).^3 -(4*rpj.^2+lbj.^2).^1.5)-2*rpj
        
        r_mean = (rpi+rpj)./2
        beta_corr = r_mean./(r_mean+ sqrt(deltai.^2 + deltaj.^2)/2) + 4*(Di+Dj)./2/sqrt(ci.^2+cj.^2)./r_mean
        beta = 2*8*pi.*(rpi+rpj)./2*(Di+Dj)./2/beta_corr*1e6
    
        return beta
    end
    return c_kernel_fuchs
end
