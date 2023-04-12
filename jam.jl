
#=**************************************************************************

The Julia Aerosol Module (JAM) acts as the main function for a 
discrete-sectional aerosol model.  The model simulates the evolution
of a single condensable species into a particle side distribution,
including representation of nucleation, coagulation, condensation,
external gas production, evaporation, wall loss, and other external losses.

JAM in its current form was created for the final project ENE 527 Climate Change and Aerosols
and was written to replicate a MATLAB discrete-sectional aerosol model
Li, C., & Cai, R. (2020). Tutorial: The discrete-sectional method to simulate an evolving aerosol. 
Journal of Aerosol Science, 150, 105615. https://doi.org/10.1016/j.jaerosci.2020.105615
Original MATLAB model available at: https://github.com/chenxi20JT/discrete-sectional-code
Translated to Julia by Obin Sturm (psturm@usc.edu) November 2022

Jam may change in the future to no longer replicate the MATLAB discrete-sectional model.
**************************************************************************=#

module Jam
include("makecollisionkernel.jl")
include("calcoagcoeffs.jl")
include("calsink.jl")
include("calevap.jl")
include("GDEqns.jl")

using QuadGK
using HCubature
using DelimitedFiles
using DifferentialEquations

function main(MP)
    # construct the collision kernel function
    collisionKernel = makecollisionkernel(MP.KERNEL,MP.T,MP.P,MP.ρ,MP.mMono)
    experiment_name = MP.Prefix*"coag_"*MP.KERNEL*"_ND"*string(MP.ND)*"_NT"*string(MP.NT)*"_"*string(MP.binsPer2)

    # if desired, calculate the coagulation coefficients in the current simulation
    if MP.CAL_COEFFS 
        coagSnkCoefs, coagSrcIds, coagSrcCoefs = calcoagcoeffs(MP.ND,MP.DIS_SECT,MP.sLim,MP.sLim0,
                                                               MP.sLen,MP.κ,MP.binsPer2,MP.KERNEL,collisionKernel,
                                                               MP.GROW_BEY_BOUND,MP.Prefix)
    else
        if isfile(experiment_name*"/coefficients/coagSnkCoefs.txt") && isfile(experiment_name*"/coefficients/coagSrcCoefs.txt") && isfile(experiment_name*"/coefficients/coagSrcIds.txt")
            coagSnkCoefs = readdlm(experiment_name*"/coefficients/coagSnkCoefs.txt",',')
            coagSrcCoefs = readdlm(experiment_name*"/coefficients/coagSrcCoefs.txt",',')
            coagSrcIds   = readdlm(experiment_name*"/coefficients/coagSrcIds.txt",',',Int)
            println("Coagulation coefficients loaded")
        else
            println("Coagulation coefficient not found. Recalculating...")
            coagSnkCoefs, coagSrcIds, coagSrcCoefs = calcoagcoeffs(MP.ND,MP.DIS_SECT,MP.sLim,MP.sLim0,
                                                                   MP.sLen,MP.κ,MP.binsPer2,MP.KERNEL,collisionKernel,
                                                                   MP.GROW_BEY_BOUND,MP.Prefix)
        end
    end

    if MP.CAL_SINK_PARAMS
        wall_coeffs,l0_coeffs,e0_coeffs = calsink(MP.NT,MP.ND,MP.sLim,MP.sLim0,MP.sLen,MP.A,MP.KERNEL,collisionKernel,MP.binsPer2,MP.Prefix)

    else
        if isfile(experiment_name*"/coefficients/e0_coeffs.txt") && isfile(experiment_name*"/coefficients/l0_coeffs.txt") && isfile(experiment_name*"/coefficients/wall_coeffs.txt")
            e0_coeffs = readdlm(experiment_name*"/coefficients/e0_coeffs.txt",',')
            l0_coeffs = readdlm(experiment_name*"/coefficients/l0_coeffs.txt",',')
            wall_coeffs   = readdlm(experiment_name*"/coefficients/wall_coeffs.txt",',')
            println("Loss coefficients loaded")
        else
            println("Loss coefficients not found. Recalculating...")
            wall_coeffs,l0_coeffs,e0_coeffs = calsink(MP.NT,MP.ND,MP.sLim,MP.sLim0,MP.sLen,MP.A,MP.KERNEL,collisionKernel,MP.binsPer2,MP.Prefix)
        end
    end

    #-------------------------compute distribution----------------------------#
    tspan = (0.0,MP.totalTime)
    u0 = MP.initConc
    p = ( MP.RR,MP.κ,MP.ND,MP.NT,coagSnkCoefs,coagSrcIds,coagSrcCoefs,wall_coeffs,l0_coeffs,e0_coeffs,MP.wLoss,MP.svgLoss,MP.satPressure,MP.dilution,MP.CONST_MONOMER,MP.COAG_OFF )
    
    odeFunc = ODEProblem(GDEqns!,u0,tspan,p)
    print("Solving general dynamic equations...")
    sol = solve(odeFunc,QNDF(autodiff=false),abstol=1e-6,reltol=1e-6,saveat=10)
    println("\t completed")

    resultsMass = hcat(sol.u...)'
    resultsNum = zeros(MP.NT,MP.NT)
    for i = 1:MP.nStep
        resultsNum[i,:] = resultsMass[i,1:MP.NT]./MP.sAvg;
    end

    if ~any("results".==readdir(experiment_name)) ;  mkdir(experiment_name*"/results") end
    writedlm( pwd()*"/"*experiment_name*"/results/resultsMass.txt",  resultsMass, ',')
    writedlm( pwd()*"/"*experiment_name*"/results/resultsNum.txt",  resultsNum, ',')
    # cp( "simulationEntry.jl", pwd()*"/"*experiment_name*"/results/"*MP.Prefix*"simulationEntry.jl",force=true)
    
    return resultsMass, resultsNum
end

export main

end


