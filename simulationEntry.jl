#=************************************************************************** 
 This Julia module is based on the input of a MATLAB discrete-sectional aerosol model:
 Li, C., & Cai, R. (2020). Tutorial: The discrete-sectional method to simulate an evolving aerosol. 
 Journal of Aerosol Science, 150, 105615. https://doi.org/10.1016/j.jaerosci.2020.105615
 Original MATLAB model available at: https://github.com/chenxi20JT/discrete-sectional-code

 INPUT:
   All the input parameters are contained in the Julia module MP, short for model parameters.
   The fields of MP that require user input are listed below:
      resultsFileName : filename in which the output is saved
      Prefix          : a string prefix that is added to filename for further distinction between simulations
      ND              : the number of discrete bins
      NS              : number of sectional bins
      DIS_SECT        : 1 -> discrete-sectional model; 2 -> pure sectional model; 3 -> pure discrete model
      sLim0           : the left limit of sectional bin. The input value is used only for the pure sectional
                        model; otherwise; sLim0 is overwritten by ND+0.5
      binsPer2        : (2)**(1/binsPer2) is the geometric factor between neighbouring sectional bins
      totalTime       : total simulation time; in sec()
      nStep           : total number of output time step
      wLoss           : the constant part of wall loss rate;in s^-1
      scgLoss         : loss constant of scavenging by pre-existing particles
      dilution        : dilution constant; in s^-1
      satPressure     : saturation vapor pressure in units of #/cm^3
      wLoss; scgLoss; dilution; satPressure correspond to the second column of Table S1 in the tutorial. 
      KERNEL          : "fuchs": Fuchs transition regime limiting sphere collision theory
      CAL_COEFFS      : true-> recalculate coefficients between bins; false-> load previous coefficients
      CAL_SINK_PARAMS : true-> recalculate coefficients for sink processes; false-> load previous coefficients
      CONST_MONOMER   : true-> the monomer concentration will remain constant throught the simulation
      COAG_OFF        : true-> coagualtion between particles will be turned off    
      GROW_BEY_BOUND  : true -> particles are allowed to grow out of the upper boundary; 
                        false -> no mass loss of the whole system due to particle growth 
      odeOptions      : solver options for ode15s()
      initConc        ; initial mass concentration
      RR              : monomer production rate, in #/(cm^3*s)
      T               : Temperature in K
      P               : Pressure in Pa
      ρ               : bulk density in kg/m^3
      mMono           : monomer mass in kg
      surfaceTension  : suface tension of nulceation species; in N/m

 OUTPUT:
      resultsMass    : particle mass in each bin in molecule/cm^3
      resultsNum     : particle number in each bin in particle/cm^3
**************************************************************************=#

module MP

# specify the output filename
resultsFileName = "test0"
Prefix      = "test0_"

# Set up the simulation bins and section equations
ND       = 100
binsPer2 = 10 # how many bins until volume is doubled? note κ = 2^(1/binsPer2)
# note: MATLAB model says diameter instead of volume, this does not seem correct
NS       = 40*binsPer2
NT       = ND + NS # total number of bins
sLim0    = 0   # will be overwitten by ND+0.5 when ND > 1

# specify simulation time
totalTime = 3000 #60*60*4
nStep     = 301

# specify particle sink parameters
wLoss  = 0
svgLoss = 0
dilution = 0
satPressure  = 0

# true or false parameters
CONST_MONOMER   = false
COAG_OFF        = false
CAL_COEFFS      = true
CAL_SINK_PARAMS = true
GROW_BEY_BOUND  = false
                             
# choose a collision kernel for particles
KERNEL  = "fuchs"  # currently the only option

odeOptions = ["list of options"]
#set the initial number distribution of particles
initConc    = zeros(NT)
initConc[1] = 0 #1e8

#set the monomer generation term
RR = 5e6
                          
#Physical parameter of the system
T       = 293.15   
P       = 1.01e5   
ρ     = 1.47e3   
mMono   = 2.4e-25  
surfaceTension = 67.5e-3 

# calculate other model parameters (previously refinemodelparams.m)
if ND == 0
    error("There has to be at least one discrete bin for monomers, ND≥1")
elseif (ND > 1) && (NS > 0)
    DIS_SECT = 1   # discrete-sectional model
elseif (ND == 1) && (NS > 0)
    DIS_SECT = 2   # pure sectional model
elseif (ND > 0) && (NS == 0)
    DIS_SECT = 3   # pure discrete model
else
    error("Specify number of bins ND≥1 and NS≥0")
end


vMono  = mMono/ρ # monomer volume
dMono  = (vMono*6/pi)^(1/3)    # monomer diameter
A      = pi*dMono^2*surfaceTension/1.5/1.38e-23/T #calculate the dimensionless surface tension

outputTime = range(0, totalTime,nStep)

κ = 2^(1/binsPer2)
sLim  = zeros(NT)  # right limit of each bin, in number of molecules
sLen  = zeros(NT)  # length of each bin, in number of molecules
sAvg  = zeros(NT)  # average particle mass in each bin, in number of molecules
convertNumToLog = zeros(NT)  # conversion factor from number per bin to dN/dlog10(Dp)

if (DIS_SECT == 1) || (DIS_SECT == 3)
    for i = 1:ND
        sLim[i] = i
        sAvg[i] = i
        convertNumToLog[i] = 3/log10((i+0.5)/(i-0.5))
    end
    sLen[1:ND] .= 1
    sLim0 = sLim[ND] + 0.5  # overwrite sLim0
end
if (DIS_SECT == 2)  # Note: the original MATLAB model cannot handle a purely sectional model 
    sAvg[ND] = ND
    convertNumToLog[ND] = 3/log10((ND+0.5)/(ND-0.5))
    sLim[ND] = ND
    sLim0 = sLim[ND]
end


if NS > 0
    sLim[ND+1]  = sLim0*κ
    sLen[ND+1]  = sLim[ND+1]-sLim0
    sAvg[ND+1]  = sLen[ND+1]/log(κ)
    convertNumToLog[ND+1]= 3/log10(κ)
end
if NS > 1
    for i = ND+2:NT
        sLim[i]  = sLim[i-1]*κ
        sLen[i]  = sLim[i]-sLim[i-1]
        sAvg[i]  = sLen[i]/log(κ)
        convertNumToLog[i] = 3/log10(κ)
    end
end
dpBins = (sAvg*vMono).^(1/3).*(6/pi)^(1/3)
end


jam = include("jam.jl")

resultsMass, resultsNum = jam.main(MP)

# using BenchmarkTools
# @btime jam.main(MP) evals=5