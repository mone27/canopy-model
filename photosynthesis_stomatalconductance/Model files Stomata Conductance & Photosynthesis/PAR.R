PAR = function(leaf,fsds){
  
  
params = list()
ground = list()
  
# --- Waveband indices for visible and near-infrared
  
params$vis = 1; params$nir = 2
 
# Solar radiation (W/m2)

# switch leaftype
# case 'sun'
#fsds = 800;   # Sun leaf

# case 'shade'
# fsds = 300;   # Shade leaf
# end
fsds = Hainich5Days$SW_IN_F

# par to W m^-2 ?
# radiation replaced

atmos$swskyvis = 0.5 * fsds;   # short wave sky

# --- Ground variables

ground$albsoi[params$vis] = 0.1;      # Soil albedo (visible waveband)

# --- Radiation absorbed by leaf (from gourp 3)

# vis Solar radiation incident on leaf

flux$swincvis = atmos$swskyvis * (1 + ground$albsoi[params$vis]);

# vis Solar radiation absorbed by leaf

flux$swflxvis = flux$swincvis * (1 - leaf$rho[params$vis] - leaf$tau[params$vis]);
flux$apar = flux$swflxvis * 4.6;

return(flux$apar)

}