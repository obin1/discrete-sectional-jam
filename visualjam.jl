# visualjam.jl

using Plots
using DelimitedFiles


timing = [110,9]
p = bar([1,2],timing,xticks = ((1,2),["MATLAB","JAM"]),fill=["orange", "purple"],legend = false)
ylabel!("Wall time [seconds]")
savefig(p,"barplot.pdf")




# Mass balance experiment
q_massbal0 = readdlm("massbal0_coag_fuchs_ND100_NT500_10/results/resultsMass.txt",',')
q_massbal5e5 = readdlm("massbal5e5_coag_fuchs_ND100_NT500_10/results/resultsMass.txt",',')
q_wloss = readdlm("wloss_coag_fuchs_ND100_NT500_10/results/resultsMass.txt",',')
# q_evap0 = readdlm("evap0_coag_fuchs_ND100_NT500_10/results/resultsMass.txt",',')
plot_massbal = plot(0:10:1000,sum(q_massbal5e5[1:101,:],dims=2),label=s"$R_1 = 5 \times 10^5 \; cm^{-3}s^{-1}$",ylims=(0,6e8))
plot_massbal = plot!(0:10:1000,sum(q_massbal0[1:101,:],dims=2),label=s"$R_1 = 0 \; cm^{-3}s^{-1}$",legend=:topleft,xlabel="time [seconds]")
plot_massbal = plot!(0:10:1000,sum(q_wloss[1:101,:],dims=2),label=s"wall loss, $w_0 = 0.01$",ylabel=s"concentration [$molec \;\; cm^{3}$]")
savefig("massbal.pdf")
# plot_massbal = plot!(0:10:1000,sum(q_evap0[1:101,:],dims=2),label=s"evap, $w_0 = 2.98 $^{-1}$")


# Trimer growth from condensation
# dNdlog10dp = 

# 
q_evap0 = readdlm("evap0_coag_fuchs_ND100_NT500_10/results/resultsMass.txt",',')

psd = plot(MP.dpBins[2:end]*1e6,q_evap0[end,2:end],xaxis=:log,xlims=(1e-4,40),label="particle phase",xlabel=s"diameter [$μm$]",legend=:topleft)
psd = scatter!([MP.dpBins[1]*1e6],[q_evap0[end,1]],label="gas phase",ylabel=s"concentration [$molec \;\; cm^{3}$]",title="PSD after 4 hours")
savefig("PSD4hours.pdf")

anim = @animate for i ∈ 1:size(q_evap0)[1]
    plot(MP.dpBins[2:end]*1e6,q_evap0[i,2:end],xaxis=:log,xlims=(1e-4,40),ylims=(0,1e8),label="particle phase",xlabel=s"diameter [$μm$]",legend=:topleft)
    scatter!([MP.dpBins[1]*1e6],[q_evap0[i,1]],label="gas phase",ylabel=s"concentration [$molec \;\; cm^{3}$]",title="PSD")
end
gif(anim, "anim_fps150.gif",fps=300)