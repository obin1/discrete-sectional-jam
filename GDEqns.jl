#=************************************************************************** 
This function represents the general dynamic equations (Gelbard and Seinfeld, 1979)
in the canonical form for ordinary differential equations for the Julia suite 
DifferentialEquations.jl: https://diffeq.sciml.ai/stable/

Gelbard, F., & Seinfeld, J. H. (1979). The general dynamic equation for aerosols. 
Theory and application to aerosol formation and growth. Journal of Colloid and
Interface Science, 68(2), 363–382. https://doi.org/10.1016/0021-9797(79)90289-3.

Based on a function of the same name from the MATLAB discrete-sectional aerosol model:
 Li, C., & Cai, R. (2020). Tutorial: The discrete-sectional method to simulate an evolving aerosol. 
 Journal of Aerosol Science, 150, 105615. https://doi.org/10.1016/j.jaerosci.2020.105615
 Original MATLAB model available at: https://github.com/chenxi20JT/discrete-sectional-code
**************************************************************************=#

function GDEqns!(du,u,p,t)
    

    R,κ,ND,NT,coagSnkCoefs,coagSrcIds,coagSrcCoefs,wall_coeffs,l0_coeffs,e0_coeffs,w,l0,e0,m0,CONST_MONOMER,COAG_OFF = p
    
    for i = 1:NT
        du[i] = 0
    end
    # du = zeros(NT)
    
    if COAG_OFF
        coagSnkCoefs[2:end,2:end] .= 0
        coagSrcCoefs[2:end,2:end] .= 0
    end

    # coagulation and condensation sink
    for i = 1:NT
        du[i] = du[i] - sum( coagSnkCoefs[i,:].*u ) * u[i]
    end

    # coagulation and condensation source
    for i = 1:NT
        for j = 1:i
            if j == i 
                preFactor = 0.5 # to avoid double counting 
            else
                preFactor = 1.0
            end
            if coagSrcIds[i,j] !=0
                du[ coagSrcIds[i,j] ] = du[ coagSrcIds[i,j] ] + preFactor * coagSrcCoefs[i,j]*u[i]*u[j]
                if coagSrcIds[i,j] < NT
                    du[coagSrcIds[i,j]+1] = du[coagSrcIds[i,j]+1] + preFactor * (coagSnkCoefs[i,j]+coagSnkCoefs[j,i]-coagSrcCoefs[i,j]) * u[i]*u[j]
                end
            end
        end
    end

    # evaporation, monomer production, wall loss, and scavenging by pre-existing particles

    # monomer rate
    if CONST_MONOMER
        du[1] = 0
    else
        du[1] = du[1] + R - w*wall_coeffs[1]*u[1] - l0*l0_coeffs[1]*u[1] - m0*u[1] + sum(u.*e0_coeffs*e0)
    end

    # other bin rates
    du[2] = du[2] - w*wall_coeffs[2]*u[2] - l0*l0_coeffs[2]*u[2] - m0*u[2] - u[2]*e0_coeffs[2]*e0 + u[3]*e0_coeffs[3]*e0*2

    for i = 3:ND-1
        du[i] = du[i] -w*wall_coeffs[i]*u[i]-l0*l0_coeffs[i]*u[i] - m0*u[i] - u[i]*e0*e0_coeffs[i]*i + u[i+1]*e0*e0_coeffs[i+1]*i
    end
    
    if ND < NT # if there are sectional bins
        du[ND] = du[ND] - w*wall_coeffs[ND]*u[ND] - l0*l0_coeffs[ND]*u[ND] - m0*u[ND] - u[ND]*e0*e0_coeffs[ND]*ND + (e0*e0_coeffs[ND+1]*u[ND+1])/(κ-1)
        for i = ND+1:NT-1
            du[i] = du[i] - w*wall_coeffs[i]*u[i] - l0*l0_coeffs[i]*u[i] - m0*u[i] - u[i]*e0*e0_coeffs[i]*(1+1/(κ-1)) + u[i+1]*e0*e0_coeffs[i+1]/(κ-1)
        end
        du[NT] = du[NT] - w*wall_coeffs[NT]*u[NT] - l0*l0_coeffs[NT]*u[NT] - m0*u[NT] - u[NT]*e0*e0_coeffs[NT]*(1+1/(κ-1))
    else # if there are only discrete bins
        du[ND] = du[ND] - w*wall_coeffs[ND]*u[ND] - l0*l0_coeffs[ND]*u[ND] - m0*u[ND] - u[ND]*e0*e0_coeffs[ND]*ND
    end

    return du

end
