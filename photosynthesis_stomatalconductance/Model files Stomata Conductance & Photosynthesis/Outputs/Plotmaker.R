

Hainich5Days = read.csv("5_days.csv",header = T, dec = ",", sep = ";")

####normal output####
d1 <- read.table("testoutputs2_with_PAR", header = T, sep=",", dec = ".")
pdf (file="Figures/Photosynthesis Model Outputs.pdf", title= "Photosynthesis Moddel Output (unmodified)")

plot(d1$an ~ c(1:240), xlab = "Timesteps", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main= "Net assimilation over time")
plot(d1$gs ~ c(1:240), xlab = "Timesteps", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main= "stomata conductance over time")
plot(d1$an ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over temperature")
plot(d1$gs ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main ="stomata conductance over temperature")
plot(d1$an ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over ambient CO2 concentration")
plot(d1$gs ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)" , main ="stomata conductance over ambient CO2 concentration")
plot(d1$an ~ Hainich5Days$SW_IN_F, xlab = "Incoming shortwave radiation (W/m2)", ylab ="Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over shortwave radiation")

dev.off()
####modified output (vcmaxha = 58520)#####

d2 <- read.table("testoutputs5_vxmaxhd_mod", header = T, sep=",", dec = ".")
pdf (file="Figures/Photosynthesis Model Outputs modified (vcmaxha=58520).pdf", title= "Photosynthesis Moddel Output (vcmaxha = 58520)")

plot(d2$an ~ c(1:240), xlab = "Timesteps", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main= "Net assimilation over time", sub= "Vcmax activation energyy lowered (58520 J/mol)")
plot(d2$gs ~ c(1:240), xlab = "Timesteps", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main= "stomata conductance over time", sub= "Vcmax activation energyy lowered (58520 J/mol)")
plot(d2$an ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over temperature", sub= "Vcmax activation energyy lowered (58520 J/mol)")
plot(d2$gs ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main ="stomata conductance over temperature", sub= "Vcmax activation energyy lowered (58520 J/mol)")
plot(d2$an ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over ambient CO2 concentration", sub= "Vcmax activation energyy lowered (58520 J/mol)")
plot(d2$gs ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)" , main ="stomata conductance over ambient CO2 concentration", sub= "Vcmax activation energyy lowered (58520 J/mol)")
plot(d2$an ~ Hainich5Days$SW_IN_F, xlab = "Incoming shortwave radiation (W/m2)", ylab ="Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over shortwave radiation", sub= "Vcmax activation energyy lowered (58520 J/mol)")

dev.off()

#### modified output (Vcmax25 = 120)####
d3 <- read.table("testoutputs6_vxmax25_mod_120", header = T, sep=",", dec = ".")
pdf (file="Figures/Photosynthesis Model Outputs modified (vcmax25 = 120).pdf", title= "Photosynthesis Moddel Output (vcmax25 = 120)")

plot(d3$an ~ c(1:240), xlab = "Timesteps", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main= "Net assimilation over time", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")
plot(d3$gs ~ c(1:240), xlab = "Timesteps", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main= "stomata conductance over time", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")
plot(d3$an ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over temperature", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")
plot(d3$gs ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main ="stomata conductance over temperature", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")
plot(d3$an ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over ambient CO2 concentration", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")
plot(d3$gs ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)" , main ="stomata conductance over ambient CO2 concentration", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")
plot(d3$an ~ Hainich5Days$SW_IN_F, xlab = "Incoming shortwave radiation (W/m2)", ylab ="Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over shortwave radiation", sub= "Vcmax25 doubled (120 μmol m–2 s–1)")

dev.off()

#### modified output (g1 = 100) ####
d4 <- read.table("testoutputs7_g1_100", header = T, sep=",", dec = ".")
pdf (file="Figures/Photosynthesis Model Outputs modified (g1=100).pdf", title= "Photosynthesis Moddel Output (g1=100)")

plot(d4$an ~ c(1:240), xlab = "Timesteps", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main= "Net assimilation over time", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")
plot(d4$gs ~ c(1:240), xlab = "Timesteps", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main= "stomata conductance over time", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")
plot(d4$an ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over temperature", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")
plot(d4$gs ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)", main ="stomata conductance over temperature", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")
plot(d4$an ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over ambient CO2 concentration", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")
plot(d4$gs ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)" , main ="stomata conductance over ambient CO2 concentration", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")
plot(d4$an ~ Hainich5Days$SW_IN_F, xlab = "Incoming shortwave radiation (W/m2)", ylab ="Leaf net photosynthesis (umol CO2/m2 leaf/s)", main ="Net assimilation over shortwave radiation", sub= "Ball-Berry slope of conductance-photosynthesis relationship(g1) = 100")

dev.off()
