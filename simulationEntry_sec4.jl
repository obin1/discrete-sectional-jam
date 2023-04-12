module MP

# specify the output filename
resultsFileName = "sec4"
Prefix      = "sec4_"

# Set up the simulation bins and section equations
ND       = 1
binsPer2 = log10(2)/log10(64) #10 # how many bins until volume is doubled? note κ = 2^(1/binsPer2)
# note: MATLAB model says diameter instead of volume, this does not seem correct
NS       = 4 #40*binsPer2
NT       = ND + NS # total number of bins
sLim0    = 0   # will be overwitten by ND+0.5 when ND > 1

MOSAIC_4bin = 1 # if 1, forego the usual model settings to use a 4 bin scheme:  
# (0.039–0.156, 0.156–0.625, 0.625–2.500, and 2.5–10.0 µm diameters)
# MOSAIC uses 4 bins in WRF-Chem, idea from Chang et al, ACP 2021, https://doi.org/10.5194/acp-21-4403-2021

# specify simulation time
totalTime = 60*60*4
nStep     = 301

# specify particle sink parameters
wLoss  = 0
svgLoss = 0
dilution = 0
satPressure  = 2e4

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
initConc[1] = 1e8

#set the monomer generation term
RR = 0
                          
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
sLim  = zeros(1, NT)  # right limit of each bin, in number of molecules
sLen  = zeros(1, NT)  # length of each bin, in number of molecules
sAvg  = zeros(1, NT)  # average particle mass in each bin, in number of molecules
convertNumToLog = zeros(1, NT)  # conversion factor from number per bin to dN/dlog10(Dp)

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
if MOSAIC_4bin == 1
    dpLim = [0.039, .156, 0.625, 2.5, 10]*1e-6
    sLim = vcat(1, ((pi/6).*dpLim[2:end].^3)./vMono)
    sLim0 = ((pi/6).*dpLim[1].^3)/vMono
    κ = (dpLim[end]/dpLim[end-1])^3
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
# resultsMass, resultsNum = jam.main(MP)
resultsMass, resultsNum = jam.main(MP)

