
#=************************************************************************** 
This function calculates beta_(i,j)*exp(1.5*A*(k**(2/3)-(k-1)**(2/3)))
which is the scaling factor of the evaporation rate

Based on a function of the same name from the MATLAB discrete-sectional aerosol model:
 Li, C., & Cai, R. (2020). Tutorial: The discrete-sectional method to simulate an evolving aerosol. 
 Journal of Aerosol Science, 150, 105615. https://doi.org/10.1016/j.jaerosci.2020.105615
 Original MATLAB model available at: https://github.com/chenxi20JT/discrete-sectional-code
**************************************************************************=#

function calevap(c_kernel,A,k)
    E_k = c_kernel(1,k-1)*exp(1.5*A*(k^(2/3)-(k-1)^(2/3)))
    return E_k
end