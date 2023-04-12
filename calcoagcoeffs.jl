#=************************************************************************** 
Based on a function of the same name from the MATLAB discrete-sectional aerosol model:
 Li, C., & Cai, R. (2020). Tutorial: The discrete-sectional method to simulate an evolving aerosol. 
 Journal of Aerosol Science, 150, 105615. https://doi.org/10.1016/j.jaerosci.2020.105615
 Original MATLAB model available at: https://github.com/chenxi20JT/discrete-sectional-code

Original intro reproduced below:
This function calculate the mass based coagulation coefficient for every
two size bin
Discrete & sectional bins are calculated separately but not
distinguished in the output matrix [because ODE do not distinguish them]
Output: coagSnkCoefs; coagSrcIds; coagSrcCoefs
coagSnkCoefs:stores aij described in section 2.2.2 of the tutorial
coagScrIds: stores u described in section 2.2.2 of the tutorial 
coagSrcCoefs: stores bij->u described in section 2.2.2 of the tutorial
Note that coagSnkCoefs also incoporates condensation described in 
section 2.2.3 of the tutorial to make the code compact

**************************************************************************=#


function calcoagcoeffs(ND, DIS_SECT, sLim, sLim0, sLen, κ, binsPer2, KERNEL, coag_kernel, GROW_BEY_BOUND, Prefix)
    NT = length(sLim)

    # initalize output variables
    coagSnkCoefs, coagSrcCoefs, = [-9999*ones(NT,NT) for _ = 1:2]
    coagSrcIds = -999*ones(Integer,NT,NT)

    # D-D coagulation
    if (DIS_SECT == 1) || (DIS_SECT == 3) # there are discrete bins beyond monomers
        print("Calculating D-D coagulation...")
        for ii = 1:ND
            for jj = 1:ii
                β = coag_kernel(ii,jj)
                coagSnkCoefs[ii,jj] = β./jj # loss for ii bin
                coagSnkCoefs[jj,ii] = β./ii # loss for jj bin

                nk = ii + jj
                # determine the corresponding bin of the newly formed particle
                if nk <= NT
                    idk = findfirst(nk .<= sLim)   # the right limit belongs to a given sectional bin, but not the left limit
                    coagSrcIds[ii,jj] = idk
                    coagSrcCoefs[ii,jj] = coagSnkCoefs[ii,jj] + coagSnkCoefs[jj,ii]
                elseif GROW_BEY_BOUND == false # remain in the largest bin?
                    coagSrcIds[ii,jj] = NT
                    coagSrcCoefs[ii,jj] = coagSnkCoefs[ii,jj] + coagSnkCoefs[jj,ii]
                    # else coagSrcIds[ii,jj] remains 0
                end
                coagSrcIds[jj,ii] = coagSrcIds[ii,jj]
                coagSrcCoefs[jj,ii] = coagSrcCoefs[ii,jj]
            end # of jj loop
        end     # of ii loop
    println("\tcompleted")
    end

    # D-S coagulation
    if (DIS_SECT == 1) || (DIS_SECT == 2)   # there are sectional bins
        print("Calculating D-S coagulation...")
        # Condensation: when ii is 1
        infoDisSec = "1/"*string(ND)
        print(infoDisSec)
        for jj = ND+1:NT # all sectional bins
            if jj-1==ND
                coagSnkCoefs[1,jj] = quadgk(x->(coag_kernel(1,x)/x),sLim0,sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
                # @btime(quadgk(x->(coag_kernel(1,x)/x),100.5, 107.71323298489746,atol=0,rtol=1e-6)[1]/7.213232984897459) 5 μs compared to 200 μs for integral in THATLAB®
                coagSnkCoefs[jj,1] = quadgk(x->(coag_kernel(1,x)),sLim0,sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
            else
                coagSnkCoefs[1,jj] = quadgk(x->(coag_kernel(1,x)/x),sLim[jj-1],sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
                coagSnkCoefs[jj,1] = quadgk(x->(coag_kernel(1,x)),sLim[jj-1],sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
            end
            coagSrcIds[1,jj] = jj
            if jj == NT && GROW_BEY_BOUND == false
                coagSrcCoefs[1,jj] = coagSnkCoefs[1,jj] + coagSnkCoefs[jj,1]
            else
                coagSrcCoefs[1,jj] = coagSnkCoefs[1,jj] + coagSnkCoefs[jj,1] - coagSnkCoefs[1,jj]/(1-1/κ)
            end
            coagSrcIds[jj,1] = coagSrcIds[1,jj]
            coagSrcCoefs[jj,1] = coagSrcCoefs[1,jj]
        end

        if DIS_SECT == 1    # there are both discrete and sectional bins
            # coagulation
            for ii = 2:ND # all discrete-sectional coagulation beyond ii=1 (that was condensation)
                for nn = 1 : length(infoDisSec)
                    print("\b")
                end
                infoDisSec = string(ii)*"/"*string(ND)
                print(infoDisSec)
                for jj = ND+1:NT
                        if jj-1 == ND
                            coagSnkCoefs[ii,jj] =  quadgk(x->(coag_kernel(ii,x)/x),sLim0,sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
                            coagSnkCoefs[jj,ii] =  quadgk(x->(coag_kernel(ii,x)/ii),sLim0,sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
                        else
                            coagSnkCoefs[ii,jj] =  quadgk(x->(coag_kernel(ii,x)/x),sLim[jj-1],sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
                            coagSnkCoefs[jj,ii] =  quadgk(x->(coag_kernel(ii,x)/ii),sLim[jj-1],sLim[jj],atol=0,rtol=1e-6)[1]/sLen[jj]
                        end            
                    # determine the size range of new particles, (nMin, nMax] and find the new bin(s)
                    if jj-1 == ND
                        nMin = ii + sLim0
                    else
                        nMin = ii + sLim[jj-1]
                    end
                    nMax = ii + sLim[jj]
                    if nMin <= sLim[NT] # within the whole size range?
                        idk = findfirst(nMin .<= sLim) 
                        coagSrcIds[ii,jj] = idk
                        if (nMax <= sLim[idk]) ||                     # all new particles in one size bin
                               (idk == NT && GROW_BEY_BOUND == false)   # or they are all put in the largest size bin
                            coagSrcCoefs[ii,jj] = coagSnkCoefs[ii,jj] + coagSnkCoefs[jj,ii]
                        else # new particles get sorted into two adjacent bins
                            if jj-1 == ND
                                coagSrcCoefs[ii,jj] = quadgk(x->(coag_kernel(ii,x)*(ii+x)/ii/x),sLim0,sLim[idk]-ii,atol=0,rtol=1e-6)[1]/sLen[jj]
                            else
                                coagSrcCoefs[ii,jj] = quadgk(x->(coag_kernel(ii,x)*(ii+x)/ii/x),sLim[jj-1],sLim[idk]-ii,atol=0,rtol=1e-6)[1]/sLen[jj]
                            end
                        end
                    else 
                        idk = NT
                        if GROW_BEY_BOUND == false
                            coagSrcIds[ii,jj] = NT
                            coagSrcCoefs[ii,jj] = coagSnkCoefs[ii,jj] + coagSnkCoefs[jj,ii]
                        end
                    end 
                    coagSrcIds[jj,ii] = coagSrcIds[ii,jj]
                    coagSrcCoefs[jj,ii] = coagSrcCoefs[ii,jj]
                end    # the jj loop
            end        # the ii loop
            for nn = 1 : length(infoDisSec)
                print("\b")
            end
        end    # discrete and sectional bin interactions
        print("\b \b")
        println("\b \t completed")
    end

    # S-S coagulation
    if (DIS_SECT == 1) || (DIS_SECT == 2)   # there are sectional bins
        print("Calculating S-S coagulation...")
        infoSecSec = "1/"*string(NT-ND)
        for ii = ND+1:NT # upper triangle
            infoSecSec = string(ii-ND)*"/"*string(NT-ND)
            if ii > ND + 1
                for nn = 1:length(infoSecSec)
                    print("\b")
                end
            end
            print(infoSecSec)
            for jj = ND+1:ii
                if ii-1 == ND
                    coagSnkCoefs[ii,jj] = hcubature(xy->coag_kernel(xy[1],xy[2])/xy[2], (sLim0,sLim0), (sLim[ii],sLim[jj]), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]
                elseif (ii-1 > ND) && (jj-1 == ND)
                    coagSnkCoefs[ii,jj] = hcubature(xy->coag_kernel(xy[1],xy[2])/xy[2], (sLim[ii-1],sLim0), (sLim[ii],sLim[jj]), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]      
                    coagSnkCoefs[jj,ii] = hcubature(xy->coag_kernel(xy[1],xy[2])/xy[1], (sLim[ii-1],sLim0), (sLim[ii],sLim[jj]), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]   
                else
                    coagSnkCoefs[ii,jj] = hcubature(xy->coag_kernel(xy[1],xy[2])/xy[2], (sLim[ii-1],sLim[jj-1]), (sLim[ii],sLim[jj]), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]      
                    coagSnkCoefs[jj,ii] = hcubature(xy->coag_kernel(xy[1],xy[2])/xy[1], (sLim[ii-1],sLim[jj-1]), (sLim[ii],sLim[jj]), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]                      
                end
                # determine the molecule number range of new particles, (nMin, nMax], and find the new bins
                if ii-1 == ND
                    nMin = 2*sLim0
                elseif (ii-1 > ND) && (jj-1 == ND)
                    nMin = sLim[ii-1] + sLim0
                else
                    nMin = sLim[ii-1] + sLim[jj-1]
                end
                nMax = sLim[ii] + sLim[jj]
                if nMin < sLim[NT] # within the whole size range?
                    idk = findfirst(nMin .< sLim)
                    coagSrcIds[ii,jj] = idk
                    if (nMax <= sLim[idk]) ||                     # all new particles in one size bin
                        (idk == NT && GROW_BEY_BOUND == false)   # or they are all put in the largest size bin
                     coagSrcCoefs[ii,jj] = coagSnkCoefs[ii,jj] + coagSnkCoefs[jj,ii]
                    else # new particles get sorted into two adjacent bins
                        if ii-1 == ND
                            coagSrcCoefs[ii,jj] = hcubature(xy->coag_kernel(xy[1],xy[2])*(xy[1]+xy[2])/xy[1]/xy[2], (sLim0,sLim0), (min(sLim[idk]-sLim0,sLim[ii]),min(sLim[idk]-sLim0,sLim[jj])), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]
                        elseif (ii-1 > ND) && (jj-1 == ND)
                            coagSrcCoefs[ii,jj] = hcubature(xy->coag_kernel(xy[1],xy[2])*(xy[1]+xy[2])/xy[1]/xy[2], (sLim[ii-1],sLim0), (min(sLim[idk]-sLim0,sLim[ii]),min(sLim[idk]-sLim[ii-1],sLim[jj])), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]
                        else 
                            coagSrcCoefs[ii,jj] = hcubature(xy->coag_kernel(xy[1],xy[2])*(xy[1]+xy[2])/xy[1]/xy[2], (sLim[ii-1],sLim[jj-1]), (min(sLim[idk]-sLim[jj-1],sLim[ii]),min(sLim[idk]-sLim[ii-1],sLim[jj])), rtol=1e-6, atol=0)[1]/sLen[ii]/sLen[jj]
                        end
                    end
                else
                    if GROW_BEY_BOUND == false
                        coagSrcIds[ii,jj] = NT
                        coagSrcCoefs[ii,jj] = coagSnkCoefs[ii,jj] + coagSnkCoefs[jj,ii]
                    end
                end
                coagSrcIds[jj,ii] = coagSrcIds[ii,jj]
                coagSrcCoefs[jj,ii] = coagSrcCoefs[ii,jj]
            end # jj loop
        end     # ii loop
        for nn = 1:length(infoSecSec)
            print("\b")
        end
        print("\t completed")
    end

    experiment_name = Prefix*"coag_"*KERNEL*"_ND"*string(ND)*"_NT"*string(NT)*"_"*string(binsPer2)
    println("\n current directory: "*pwd())

    if any(experiment_name.==readdir())
        if ~any("coefficients".==readdir(experiment_name))
            mkdir(experiment_name*"/coefficients")
        end
        writedlm( pwd()*"/"*experiment_name*"/coefficients/coagSnkCoefs.txt",  coagSnkCoefs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/coagSrcCoefs.txt",  coagSrcCoefs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/coagSrcIds.txt",  coagSrcIds, ',')
    else
        mkdir(experiment_name)
        mkdir(experiment_name*"/coefficients")
        writedlm( pwd()*"/"*experiment_name*"/coefficients/coagSnkCoefs.txt",  coagSnkCoefs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/coagSrcCoefs.txt",  coagSrcCoefs, ',')
        writedlm( pwd()*"/"*experiment_name*"/coefficients/coagSrcIds.txt",  coagSrcIds, ',')
    end

    return coagSnkCoefs, coagSrcIds, coagSrcCoefs

end

