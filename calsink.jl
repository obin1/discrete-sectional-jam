#=************************************************************************** 
Given model parameters, this function calculates scaling factors for 
rate constants of the following processes for both discrete & sectional bins: 
a) loss to pre-existing particles: l0_coeffs
b) wall loss: wall_coeffs
c) evaporation: e0_coeffs

Based on a function of the same name from the MATLAB discrete-sectional aerosol model:
 Li, C., & Cai, R. (2020). Tutorial: The discrete-sectional method to simulate an evolving aerosol. 
 Journal of Aerosol Science, 150, 105615. https://doi.org/10.1016/j.jaerosci.2020.105615
 Original MATLAB model available at: https://github.com/chenxi20JT/discrete-sectional-code
**************************************************************************=#



function calsink(NT,ND,slim,slim0,slen,A,kernel,c_kernel,binsPer2,Prefix)

    wall_coeffs, l0_coeffs, e0_coeffs = [zeros(NT) for _ = 1:3]

    # calculate scaling of wall loss rates 
    if ND > 0
        wall_coeffs[1:ND] = 1.0 ./ ((1:ND).^(1/3))
    end
    if ND+1 <= NT
        wall_coeffs[ND+1] = 1.5/slen[ND+1] * (slim[ND+1]^(2/3) - slim0^(2/3))
    end
    if ND+2 <= NT
        wall_coeffs[ND+2:NT] = 1.5./slen[ND+2:NT].*(slim[ND+2:NT].^(2/3) - slim[ND+1:NT-1].^(2/3))
    end
    println("wall loss parameters completed")

    # calculate scaling of loss rates to preexisting particle
    if ND > 0
        l0_coeffs[1:ND] = 1.0 ./((1:ND).^(1/2))
    end
    if ND+1 <= NT
        l0_coeffs[ND+1] = 2/slen[ND+1] * (slim[ND+1]^(1/2) - slim0^(1/2))
    end
    if ND+2 <= NT
        l0_coeffs[ND+2:NT] = 2.0 ./slen[ND+2:NT].* (slim[ND+2:NT].^(1/2) - slim[ND+1:NT-1].^(1/2))
    end
    println("loss to preexisting particles parameters completed")


    # calculate scaling of evaporation rates
    e0_coeffs[1] = 0 # for monomer, evaporation does not happen (monomer is gas phase)
    # e0_coeffs[2] is treated differently
    # *2 is to account for the factor of 2 in Eqn. 1 in the tutorial
    # 0.5 is to account for the factor of 0.5 in the collision rate constant
    # when two monomers collide
    # *2 & *0.5 cancels out,but e0_coeffs[2] is written out here explicitly for clarity
    if ND > 1
        e0_coeffs[2] = 1/2*calevap(c_kernel,A,2)*2*0.5
    end
    if ND > 2
        for i = 3:ND
            e0_coeffs[i] = 1/i*calevap(c_kernel,A,i)
        end
    end
    # sectional evaporation
    if ND+1 <= NT
        e0_coeffs[ND+1] = quadgk(x->(calevap(c_kernel,A,x)/x),slim0,slim[ND+1],atol=0,rtol=1e-6)[1]/slen[ND+1]
    end
    if ND+2 <= NT
        for i = ND+2:NT
            e0_coeffs[i] = quadgk(x->(calevap(c_kernel,A,x)/x),slim[i-1],slim[i],atol=0,rtol=1e-6)[1]/slen[i]
        end
    end
    println("evaporation parameters completed")

    experiment_name = Prefix*"coag_"*kernel*"_ND"*string(ND)*"_NT"*string(NT)*"_"*string(binsPer2)

    if any(experiment_name.==readdir())
        if ~any("coefficients".==readdir(experiment_name))
            mkdir(experiment_name*"/coefficients")
        end
        writedlm( pwd()*"/"*experiment_name*"/coefficients/wall_coeffs.txt",  wall_coeffs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/l0_coeffs.txt",  l0_coeffs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/e0_coeffs.txt",  e0_coeffs, ',')
    else
        mkdir(experiment_name)
        mkdir(experiment_name*"/coefficients")
        writedlm( pwd()*"/"*experiment_name*"/coefficients/wall_coeffs.txt",  wall_coeffs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/l0_coeffs.txt",  l0_coeffs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/e0_coeffs.txt",  e0_coeffs, ',')
    end

    return wall_coeffs, l0_coeffs, e0_coeffs


end
