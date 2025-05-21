### Hard X-ray source heights above photosphere

Hard X-ray footpoint simulations are following [Kontar et al, 2010](http://dx.doi.org/10.1088/0004-637X/717/1/250). In the thick-target approximation, the electron flux spectrum versus distance can be written

$$\dot{N}(E, s)=\dot{N}_0\left(E_0\right) \frac{E}{E_0} $$

where $N_0(E)$ is the injected spectrum of energetic electrons, taken to be (Equation 16 in https://doi.org/10.3847/1538-4357/aafad3 and [Luo et al, 2024](https://doi.org/10.3847/1538-4357/ad6a59))

$$\dot{N}_0(E)=\left(\frac{\dot{N}_0}{k_bT}\right)\frac{(k-1)(k-2)}{k^2} \frac{E /k_bT}{\left(1+E/(k_bTk)\right)^k}$$

a power law of $\propto E^{-\delta}$ at $E>> kbT$ and

$$E_0(E, s)^2=E^2+2 K \int_0^s n\left(s'\right) ds'$$

where $K=2 \pi e^2 \ln \Lambda, \ln \Lambda$ is the Coulomb logarithm, and $e$ is the electron charge. 
For neutral chromosphere $\ln \Lambda=\ln \Lambda_{e H}=7$ and for fully ionised
$\ln \Lambda=\ln \Lambda_{e H}=20$ (e.g., [Brown 1973](http://dx.doi.org/10.1007/BF00152919); [Emslie 1978](http://dx.doi.org/10.1086/156371)).

The X-ray flux spectrum emitted by the energetic electrons in a magnetic flux tube of cross-sectional area $A$ and observed at 1 AU is given as

$$I(\epsilon, s)=\frac{1}{4 \pi R^2} n(s) \int_\epsilon^{\infty} \dot{N}(E, s) \sigma(E, \epsilon) dE $$

where $\sigma(E, \epsilon)$ is the isotropic bremsstrahlung cross section, $R$ is the flare-observer distance, so the total thick target flux spectrum at 1au is
$$I(\epsilon)=\int_0^{h_{max}} \dot{N}(E, s)ds $$

The density profile was taken from the previous observations see [Kontar et al 2008](http://dx.doi.org/10.1051/0004-6361:200810719) 
showing chromospheric neutral number hydrogen density as a function of height
![image](https://github.com/user-attachments/assets/eb4db3c7-5bd7-4725-be1b-273d06a88e31)
