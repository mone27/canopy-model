install.packages("dplyr")
library(dplyr)

setwd("D:/Kathina/FES/2020-21 WS/Eco-Atmo Processes/Modelling/data")

#obtain data from the Hainich forest
data <-read.csv("FLX_DE-Hainich.csv",header=TRUE, sep=";", na.strings="NA", dec=",")
# select: air temperature TA_F (°C), vapor pressure deficit VPD_F (hPa) and air pressure PA_F (kPa)
# incoming solar SW_IN_F and longwave radiation LW_IN_F (W m-2), soil Temp.in 2cm depth TS_F_MDS_1 (°C)
# wind speed WS_F (M s-1)
mydata<- data%>%select(TIMESTAMP_START,TA_F,VPD_F,PA_F,SW_IN_F,LW_IN_F,TS_F_MDS_1,WS_F)

#change units from kPa to Pa and °C to K
mydata<- mydata%>% 
  mutate(tair = TA_F+273.15)%>% 
  mutate(VPD_F = VPD_F*100)%>% 
  mutate(patm = PA_F*1000)%>% 
  mutate(u = WS_F)%>% 
  mutate(Ts = TS_F_MDS_1+273.15)%>%
  select(TIMESTAMP_START,tair,VPD_F,patm,u,Ts,SW_IN_F,LW_IN_F)

#atmospheric parameters

#select only necessary columns
atmo <- mydata%>%
  select(TIMESTAMP_START,patm,tair,VPD_F)

#create list
atmos <- as.list(atmo) 
# Specific heat of air at constant pressure (J/mol/K)
atmos$cpair <- 29.2

#physical constants
tfrz <- 273.15                          # Freezing point of water (K)
mmh2o <- 0.01801528                     # Molecular mass of water (kg/mol)
sigma <- 5.670374419*10^-8              # Stefan-Boltzmann constant (W/m2/K4)
physcon <- data.frame(tfrz,mmh2o,sigma)
physcon <- as.list(physcon)

#leaf propperties
emiss <- 0.96                        #leaf emissity
leaf <- data.frame(emiss)
leaf <- as.list(leaf)

#fluxes

# Leaf radiative forcing qa (W/m2)
# radiative forcing Qa = incoming shortwave - 15% + incoming longwave radiation + longwave soil (Boltzman law)
# longwave radiation soil = 0.97*Boltzman constant (sigma)* soil temperature (K)^4a
qa<- mydata%>%
  select(TIMESTAMP_START,Ts,SW_IN_F,LW_IN_F)%>%
  mutate(LW_S = 0.97*physcon$sigma*Ts^4)%>% # Longwave radiation from soil
  mutate(qa = SW_IN_F-0.15*SW_IN_F + LW_IN_F + LW_S)%>% # radiative forcing
  select(TIMESTAMP_START,qa)



### boundary layer conductancee for heat gbh (mol/m2 leaf/s) for FORCED convection 
# gbh = (Nu*??m*Dh)/dl 

# required parameters:
# pm         #molar density (mol m-3) 
#ideal gas law: p = atmo$patm/R*atmo$tair   gas constant R: 0.167226 (J/Kg K)
# dl         #leaf dimension characteristic length (cm)
# Dh         #diffusivity of heat Dh (m^2 s-1)
# Nu         #Nusselt number for forced conv. 
#             Nu = 0,66* Pr^0,33* Re^0,5
# Pr         #Prandtl number (molecular diffusion properties of the fluid)
#             Pr = ??/Dh
# Re         #Reynolds number 
#             Re = u*dl/??
# u          #wind speed (M s-1)

# constants 
dl <- 0.1    #leaf dimension characteristic length (cm)

# constants needed for diffusivity calculation (from 29.0 appendices in Eco-Atmo Bonan)
T0<- 273.15  # Temperature 0°C in Kelvin
P0 <- 101325 # air Pressure (Pa) at 0°C
Dh0 <- 18.9*10^(-6) # reference diffusivity of heat Dh (m^2 s-1) at 0°C
v0 <- 13.3*10^(-6) # reference diffusivity of momentum v (m^2 s-1) at 0°C
R <- 0.0821        # gas constant R: 0.167226 (J/Kg K)
d <- 1.225         #density of dry air at 15°C and sea level (g/L)
conduct1 <- mydata%>%
  select(TIMESTAMP_START,tair,patm,u)%>%
  mutate(pm = d*R*tair)%>%                # molar density (mol m-3) according to ideal gas law ???
  mutate(Dh = Dh0*(P0/patm)*(tair/T0))%>% # diffusivity of heat (m^2 s-1) 
  mutate(v = v0*(P0/patm)*(tair/T0))      # diffusivity of momentum (m^2 s-1)

conduct2 <- conduct1%>%
  select(TIMESTAMP_START,tair,patm,u,pm,Dh,v)%>%
  mutate(Pr = v/Dh)%>%                    # Prandtl number
  mutate(Re = u*dl/v)%>%                  # Reynolds number 
  mutate(Nu = 0.66*Pr^0.33*Re^0.5)%>%     # Nusselt number for forced conv.
  mutate(gbh = (Nu*pm*Dh)/dl)             # Boundary layer conductance for heat (mol/m2 leaf/s)

#boundary layer conductance for water vapor gbw (mol m-2 s-1) 
#              gbw = 0,036* Sc^0,33* Re^0,50
# Sc         #Sherwood number
#             Sc = ??/Dj

#constants
DW0 <- 21.8*10^(-6)       #Molecular diffusivity of mass H2o at 0°C (Mm² s-1)

conduct3 <- conduct2%>%
  mutate(DW = DW0*(P0/patm)*(tair/T0)) #Molecular diffusivity of mass H2o (Mm² s-1)
conduct4 <- conduct3%>%
  mutate(Sc = v/DW)  #Sherwood number
conduct5 <- conduct4%>%
  mutate(gbw = 0.036*Sc^0.33*Re^0.50) #boundary layer conductance for water vapor (mol m-2 s-1) 

conduct <- conduct5%>%
  select(TIMESTAMP_START,gbh,gbw)

#  predefine leaf temperature
leaftemp<- mydata%>%
  select(TIMESTAMP_START,tair)%>%
  mutate(tleaf = tair)               

flux <- merge(qa, conduct, by="TIMESTAMP_START")
flux <- merge(flux, leaftemp, by="TIMESTAMP_START")
flux <- as.list(flux)

# Leaf stomatal conductance (mol H2O/m2 leaf/s)  
flux$gs <- seq(from=-0.2,to=0.1,length.out = 52608)


#function of Saturation vapor pressure and temperature derivative -- satvap()
satvap <- function (tc) {   
  # --- For water vapor (temperature range is 0C to 100C)  
  a0 =  6.11213476;        b0 =  0.444017302; 
  a1 =  0.444007856;       b1 =  0.286064092e-01;
  a2 =  0.143064234e-01;   b2 =  0.794683137e-03;
  a3 =  0.264461437e-03;   b3 =  0.121211669e-04; 
  a4 =  0.305903558e-05;   b4 =  0.103354611e-06; 
  a5 =  0.196237241e-07;   b5 =  0.404125005e-09; 
  a6 =  0.892344772e-10;   b6 = -0.788037859e-12; 
  a7 = -0.373208410e-12;   b7 = -0.114596802e-13;
  a8 =  0.209339997e-15;   b8 =  0.381294516e-16;  
  # --- For ice (temperature range is -75C to 0C)  
  c0 =  6.11123516;        d0 =  0.503277922; 
  c1 =  0.503109514;       d1 =  0.377289173e-01; 
  c2 =  0.188369801e-01;   d2 =  0.126801703e-02; 
  c3 =  0.420547422e-03;   d3 =  0.249468427e-04; 
  c4 =  0.614396778e-05;   d4 =  0.313703411e-06; 
  c5 =  0.602780717e-07;   d5 =  0.257180651e-08; 
  c6 =  0.387940929e-09;   d6 =  0.133268878e-10; 
  c7 =  0.149436277e-11;   d7 =  0.394116744e-13; 
  c8 =  0.262655803e-14;   d8 =  0.498070196e-16;    
  
  # --- Limit temperature to -75C to 100C
  
  tc = min(tc, 100);
  tc = max(tc, -75);
  
  #--- Saturation vapor pressure (esat, mb) and derivative (desat, mb)  
  if (tc >= 0){ 
    esat  = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8))))))); 
    desat = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8))))))); 
  } else {
    esat  = c0 + tc*(c1 + tc*(c2 + tc*(c3 + tc*(c4 + tc*(c5 + tc*(c6 + tc*(c7 + tc*c8))))))); 
    desat = d0 + tc*(d1 + tc*(d2 + tc*(d3 + tc*(d4 + tc*(d5 + tc*(d6 + tc*(d7 + tc*d8)))))));
  }
  esat  = esat  * 100 
  desat = desat * 100  
  return(list(esat, desat))
}


#function of Latent heat of vaporization -- latvap()
latvap<-function (tc, mmh2o){   ### same question here...
  # Latent heat of vaporization (J/mol) at temperature tc (degC)  
  val<- 2501.6 - 2.3773 * tc # Latent heat of vaporization (J/g) 
  val<- val * 1000           # Convert from J/g to J/kg 
  val<- val * mmh2o          # Convert from J/kg to J/mol  
  return(val) 
}

str(physcon)
str(atmos)
str(leaf)
str(flux)

LeafTemperature <- function(physcon, atmos, leaf, flux){
  
  # Leaf temperature and energy fluxes
  
  # ------------------------------------------------------
  # Input
  #   physcon.tfrz     ! Freezing point of water (K)
  #   physcon.mmh2o    ! Molecular mass of water (kg/mol)
  #   physcon.sigma    ! Stefan-Boltzmann constant (W/m2/K4)
  #   atmos.patm       ! Atmospheric pressure (Pa)
  #   atmos.cpair      ! Specific heat of air at constant pressure (representable value: 29.2 J/mol/K)
  #   atmos.tair       ! Air temperature (K)  
  #   atmos.VPD_F      ! vapor pressure deficit VPD_F (Pa)  
  #   leaf.emiss       ! Leaf emissivity
  #   flux.gbh         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
  #   flux.gbw         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
  #   flux.gs          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
  #   flux.qa          ! Leaf radiative forcing (W/m2 leaf)
  #
  # Input/ouput
  #   flux.tleaf       ! Leaf temperature (K)
  #
  # Output
  #   flux.rnet        ! Leaf net radiation (W/m2 leaf)
  #   flux.lwrad       ! Longwave radiation emitted from leaf (W/m2 leaf)
  #   flux.shflx       ! Leaf sensible heat flux (W/m2 leaf)
  #   flux.lhflx       ! Leaf latent heat flux (W/m2 leaf)
  #   flux.etflx       ! Leaf transpiration flux (mol H2O/m2 leaf/s)
  # ------------------------------------------------------
  
  # --- Latent heat of vaporization (J/mol)
  
  lambda_val <- latvap ((atmos$tair-physcon$tfrz), physcon$mmh2o)
  
  # --- Newton-Raphson iteration until leaf energy balance is less than
  # f0_max or to niter_max iterations
  
  niter <- 0 # Number of iterations
  f0 <- 1e36 # Leaf energy balance (W/m2)
  
  niter_max <- 100 # Maximum number of iterations
  f0_max <- 1e-06  # Maximum energy imbalance (W/m2)
  
  while ((niter < niter_max) & (abs(f0) > f0_max)){
    
    # Increment iteration counter
    
    niter <- niter + 1
    
    # Saturation vapor pressure (Pa) and temperature derivative (Pa/K)
    
    satvap_out <- satvap (flux$tleaf-physcon$tfrz)
    esat <- satvap_out[[1]]
    desat <- satvap_out[[2]];
    
    # Leaf conductance for water vapor (mol H2O/m2/s)
    
    gleaf <- flux$gs * flux$gbw / (flux$gs + flux$gbw)
    
    # Emitted longwave radiation (W/m2) and temperature derivative (W/m2/K)
    
    flux$lwrad <- 2 * leaf$emiss * physcon$sigma * flux$tleaf^4
    dlwrad <- 8 * leaf$emiss * physcon$sigma * flux$tleaf^3
    
    # Sensible heat flux (W/m2) and temperature derivative (W/m2/K)
    
    flux$shflx <- 2 * atmos$cpair * (flux$tleaf - atmos$tair) * flux$gbh
    dshflx <- 2 * atmos$cpair * flux$gbh
    
    # Latent heat flux (W/m2) and temperature derivative (W/m2/K)
    
    flux$lhflx <- lambda_val / atmos$patm * atmos$VPD_F * gleaf
    dlhflx <- lambda_val / atmos$patm * desat * gleaf
    
    # Energy balance (W/m2) and temperature derivative (W/m2/K)
    
    f0 <- flux$qa - flux$lwrad - flux$shflx - flux$lhflx
    df0 <- -dlwrad - dshflx - dlhflx
    
    # Change in leaf temperature
    
    dtleaf <- -f0 / df0
    
    # Update leaf temperature
    
    flux$tleaf <- flux$tleaf + dtleaf
    
  }
  
  # --- Net radiation
  
  flux$rnet <- flux$qa - flux$lwrad
  
  # --- Error check
  
  err <- flux$rnet - flux$shflx - flux$lhflx
  if (abs(err) > f0_max){
    cat('err  = ', err, '\n',
        'qa  = ',flux$qa, '\n',
        'lwrad  = ',flux$lwrad, '\n',
        'sh  = ',flux$shflx, '\n',
        'lh  = ',flux$lhflx)
    stop ('LeafTemperature error')
  }
  
  # Water vapor flux: W/m2 -> mol H2O/m2/s
  
  flux$etflx <- flux$lhflx / lambda_val
  
  return(list(flux$tleaf, flux$rnet, flux$lwrad, flux$shflx, flux$lhflx, flux$etflx))
}

plot(flux$tleaf, ylab = "leaf temperature (K)")
plot(flux$tair,ylab = "air temperature (K)")
plot(flux$rnet, ylab = "Leaf net radiation (W/m^2 leaf)")
plot(flux$lwrad, ylab = "Longwave radiation emitted from leaf (W/m2 leaf)")
plot(flux$shflx, ylab = "Leaf sensible heat flux (W/m2 leaf)")
plot(flux$lhflx, ylab = "Leaf latent heat flux (W/m2 leaf)")
plot(flux$etflx, ylab = " Leaf transpiration flux (mol H2O/m2 leaf/s)")
