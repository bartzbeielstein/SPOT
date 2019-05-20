###################################################################################
#' Cyclone Simulation: Mothes
#'
#' Calculate cyclone collection efficiency and pressure drop according to Mothes.
#'
#' @param x a single particle diameter (scalar)
#' @param cyclone list of cyclone's geometrical parameters
#' @param fluid list of fluid parameters
#'
#' @return returns a vector with first element being the collection efficiency, second element being the pressure drop
#'
#' @keywords internal
#' @export
###################################################################################
calculationMothes <- function(x,cyclone, fluid)
{
 	#---------------------------------------------
	#Umrechnung Geometriedaten f?r Modell Barth/ geometric data
	#-------------------------------------------
	ra = cyclone$Da/2 				#Zyklonaussenradius /cyclone radius [m]
	ri = cyclone$Dt/2  			#Tauchrohrradius /vortex finder radius [m]
	rx = ri	       			#Staubaustragsradius / discharge duct radius [m]
	h = cyclone$H        	      	#Gesamthoehe Zyklon / total cyclone height [m]
	ht = cyclone$Ht       	    		#Tauchrohreintauchtiefe / [m] Parameter 8.4
	#unused: hi = h - ht
	#epsilon = Cyclone$Epsilon*pi/180  		#Umrechnung in Bogenmass / calculation angle radians            
	#hc = Cyclone$Hc							#H?he konischer Teil / heigth conical section [m]
	he = cyclone$He       			#Schlitzeinlauhoehe in m / inlet section height [m]
	be = cyclone$Be      			#Schlitzeinlauhoehe in m / inlet section width  [m]
	hk = h-he								#Hoehe konischer Teil / heigth cylindrical section [m]
	  
	ra.alt <- ra  
	  
	e=atan((ra-rx) / hk)    #winkel epsilon, bogenmass


	vp = fluid$Vp     			#volumenstrom / volumetric gas flow rate [m/s]
	rhop = fluid$Rhop
	croh = fluid$Croh
	rhof = fluid$Rhof 
	mu = fluid$Mu
	dp = 0.0125   # Diffusionskoeff. n. Abraham - diffusionCoeffAbraham
  	rb = 0.0075  # Reibungskoeffizient n. Meissner - frictionCoeffMeissner
  
   
 
  #x = 1        # partikeldurchmesser - partileDiam
  
  
  # BN: Abbruch falls x = 1000 fehlt
  #x = x * 1e-6 
    
  #
  # fluid parameter
  #
  vis = mu
  #hz = h - (ra - rx) / tan(e)  
  hz= he
  #print(e)
  #unused: l  = ht + 0.9 * (h - ht)     # BN: Eins (1) ? NEIN: "l" charakteristische Laenge
  vk = pi * (h - hz) * (ra^2 + rx^2 + ra*rx) / 3
  vk = vk + pi * ra^2 * hz
  vr = vp / (2 * pi * ri * (h - ht))
  u  = 1 - be / ra
  # u  = acos(be / ra - l)
  ar = -1 * atan(u / sqrt(-u * u+1)) + pi/2
  vd = vp / (pi * ra^2)
  bt = -0.204 * be / ra + 0.889                    # Gl. 2.25
  ve = pi * ra^2 / (be * he * bt) * vd
  hc = (2 * pi - ar) / (2 * pi) - 1                # BN: oder l?!? Nein, eher nicht! Gl. 2.27
  hc = hz / ra + he / ra * hc                      # Gl. 2.27
  hc = vd / (rb * hc)
  ve = hc * (sqrt(0.25 + ve /hc) - 0.5)            # Gl. 2.26
  dm = ve / vd * (rb + rb / sin(e))                # Gl. 2.29
  vt = ve / ((ri / ra) * (1 + dm * (1 - ri / ra)))    # BN: dritte Klammer? Vllt l? Gl. 2.28
  rf = ra
  ra = sqrt(vk / (pi * h))                         
  ve = ve / (ra / rf * (1 + dm * (1 - ra / rf)))   # BN: vllt l? Gl. 2.28

  #
  # separation parameter
  #
  wi = rhop * x^2 * vt^2 / (18 * vis * ri) # Gl. 2.36
  wa = rhop * x^2 * ve^2 / (18 * vis * ra) # BN: "y" replaced by "x"! Gl. 2.34
  k0 = h - ht
  k1 = 2 * pi * ra * wa / vp                 # A nach Gl. 2.42
  k2 = 2 * pi * ri * dp / (vp * (ra - ri))   # B nach Gl. 2.42
  k3 = 2 * pi * ri * (wi - vr) / vp          # C nach Gl. 2.42
  
 
  #
  # cases
  # 
  if ((wi - vr) <= 0){
    a = -k0 * (-k1 + k3 - k2) - 1 #BN: "l" ?!?
    b = -k0 * k2
    c =  k0 * (-k3 + k2)
  } else {
    a = -1 + k0 * (k1 + k2) #BN: "l" ?!?
    b = k0 * (k2 - k3) # nach Beate , Gl. 
    c = k0 * k2    
  }
  d = b - 1 # BN: raus aus if-clause!

  m3 = (d + a) / 2
  # m4 = sqrt(abs(m3^2 - (a * d - b * c))) # BN: abs integrated to ensure positive numbers
  m4 = sqrt(m3^2 - (a * d - b * c)) 
  m1 = m3 + m4
  m2 = m3 - m4
  c0 = 1 #BN: "l" ?!?
  c2 = c0 * exp(-k1 * (ht - he / 2))

  
  # 
  # neglect re - up - turbulence (?!?) 
  #
  r1 = c2
  r2 = 0
  
  #
  # calculate separation degree
  #
  c4 = r1 * (m1 - a) / b + r2 * (m2 - a) / b    # BN: Gl. 2.41,2
  t = 1 - c4 / c0                               # BN: "l" ?!? Nein: Gl. 2.44
  
	#Pressure drop
 	BE = be/ra.alt
	re = ra.alt-be/2  

	Fe = be*he	            		#Flaeche Schlitzeinlauf  (Eintrittsquerschnitt) / inlet area                
	Fi = (pi*ri^2)   				#Tauchrohrquerschnittsflaeche (Austrittsquerschnitt) / 
	F = (Fe/Fi)            			#Verhaeltnis von Eintrittsflaeche Fe zu Tauchrohr     

	B = croh/rhof		  	# Gutbeladung 
	lambda = fluid$lambdag*(1+2*sqrt(B))  	# Wandreibungswert
	alpha = 1.0-(0.54-0.153/F)*BE^(1/3)    
	vi = (vp/(pi*ri^2)) 
	U = 1/(F*alpha*ri/re+lambda*h/ri)   #Umfangsgeschwindigkeit auf dem Aussenradius Gl. 2.18  


	#unused: xi1 = 0
	xi2 = (U^2*(ri/ra.alt))*(1-lambda*(h/ri)*U)^(-1)	# Widerstandsbeiwert Einlaufzone zwischen 1 und 5
	xi3 = 2+3*U^(4/3)+U^2   				# Widerstandsbeiwert Einlaufzone zwischen 10 und 50
	deltaP =(rhof/2)*vi^2*(xi2+xi3)			# N/m^2 = 16.99 hPa	
	##############################################
  return(c(t,deltaP))	
}
###################################################################################
#' Cyclone Simulation: Barth/Muschelknautz
#'
#' Calculate cyclone collection efficiency according to Barth/Muschelknautz.
#'
#' @param cyclone list of cyclone's geometrical parameters
#' @param fluid list of fluid parameters
#' @param xmean vector of middle points for the intervals of the particle size distribution
#' @param delta vector, amount of dust (percentage value) in each interval
#'
#' @return returns a vector, first element is the efficiency, second element is the pressure drop
#'
#' @keywords internal
#' @export
###################################################################################
calculationBarthMuschelknautz <- function(cyclone, fluid, xmean, delta){
	#---------------------------------------------
	#Umrechnung Geometriedaten fuer Modell Barth / geometric data
	#-------------------------------------------
	ra <- cyclone$Da/2 				#Zyklonaussenradius /cyclone radius [m]
	ri <- cyclone$Dt/2  			#Tauchrohrradius /vortex finder radius [m]
#unused	rx <- ri	       			#Staubaustragsradius / discharge duct radius [m]
#unused	hi <- ((cyclone$H-cyclone$Ht)/1000)   
	h <- cyclone$H        	      	#Gesamthoehe Zyklon / total cyclone height [m]
	ht <- cyclone$Ht       	    		#Tauchrohreintauchtiefe / [m] Parameter 8.4
#unused	epsilon <- cyclone$Epsilon*pi/180  		#Umrechnung in Bogenmass / calculation angle radians            
#unused	hc <- cyclone$Hc							#Hoehe konischer Teil / heigth conical section [m]
	he <- cyclone$He       			#Schlitzeinlauhoehe in m / inlet section height [m]
	be <- cyclone$Be       			#Schlitzeinlauhoehe in m / inlet section width  [m]
#unused	hk <- (h-he)								#Hoehe konischer Teil / heigth cylindrical section [m]
	  
  vp <- fluid$Vp     		
	croh <- fluid$Croh	
	rhof <- fluid$Rhof 
	rhop <- fluid$Rhop
	mu <- fluid$Mu

	BE <- be/ra
	re <- ra-be/2  

	Fe <- be*he	            		#Flaeche Schlitzeinlauf  (Eintrittsquerschnitt) / inlet area                
	Fi <- (pi*ri^2)   				#Tauchrohrquerschnittsflaeche (Austrittsquerschnitt) / 
	F <- (Fe/Fi)            			#Verhaeltnis von Eintrittsflaeche Fe zu Tauchrohr     

	B <- croh/rhof		  	# Gutbeladung 
	lambda <- fluid$lambdag*(1+2*sqrt(B))  	# Wandreibungswert
	alpha <- 1.0-(0.54-0.153/F)*BE^(1/3)    
	vi <- (vp/(pi*ri^2)) 
	vr <- (vp/(2*pi*ri*(h-ht)))    	#Radialgeschwindigkeit am Innenzylinder / radial gas velocity
	U <- 1/(F*alpha*ri/re+lambda*h/ri)   #Umfangsgeschwindigkeit auf dem Aussenradius Gl. 2.18 
	vphii <- U*vi  

  xGr <- (sqrt((18*mu*vr*ri)/((rhop-rhof)*vphii^2)))#*1e6;# Grenzkorn Gl. 2.4 Seite

	Tf <-  function(x){  
			  f <-(1+2/(x/xGr)^(3.564))^(-1.235);  
			  return(f)    
			  }
	#unused	xi1 <- 0
	xi2 <- (U^2*(ri/ra))*(1-lambda*(h/ri)*U)^(-1)	# Widerstandsbeiwert Einlaufzone zwischen 1 und 5
	xi3 <- 2+3*U^(4/3)+U^2   				# Widerstandsbeiwert Einlaufzone zwischen 10 und 50
	deltaP <-(rhof/2)*vi^2*(xi2+xi3)			# N/m^2 = 16.99 hPa

	
	Ew <- sum(Tf(xmean)*delta)	
	
	#x50 <- xGr*1.3 #wrong, see replacement below
	x503 <- xmean[which(cumsum(delta) >= 0.5) [1]] #median of the particle size distribution 
	#(based on intervals given by xmean and percentages given by delta)				
	ve <- (vp/Fe)	                	# Einlaufgeschwindigkeit 
	vphia <- ve*((re/ra)*(1/alpha)) 		
	BGr <-(lambda*mu*sqrt(ra*ri))/ ((1-ri/ra)*rhop*(x503)^2*sqrt(vphia*vphii)) #Grenzbeladung / critical load	
	### Gutbeladung ist groesser als die Grenzbeladung:
	E <- Ew
	if (B > BGr){
		E <- 1-BGr/B+(BGr/B)*Ew
	}
	Calculation <- c(E, Ew,deltaP)
	return(Calculation)
}


###################################################################################
#' Objective function - Cyclone Simulation: Barth/Muschelknautz
#'
#' Calculate cyclone collection efficiency. A simple, physics-based
#' optimization problem (potentially bi-objective). See the references [1,2].
#'
#' @param x vector of length at least one and up to six, specifying non-default geometrical parameters in [m]: Da, H, Dt, Ht, He, Be
#' @param deterministic binary vector. First element specifies whether volume flow is deterministic or not. Second element specifies whether particle density is deterministic or not. Third element specifies whether particle diameters are deterministic or not. Default: All are deterministic (TRUE).
#' @param cyclone list of a default cyclone's geometrical parameters: fluid$Da, fluid$H, fluid$Dt, fluid$Ht, fluid$He and fluid$Be
#' @param fluid list of default fluid parameters: fluid$Mu, fluid$Vp, fluid$Rhop, fluid$Rhof and fluid$Croh
#' @param noiseLevel list of noise levels for volume flow (noiseLevel$Vp) and particle density (noiseLevel$Rhop), only used if non-deterministic.
#' @param model type of the model (collection efficiency only): either "Barth-Muschelknautz" or "Mothes"
#' @param intervals vector specifying the particle size interval bounds.
#' @param delta vector of densities in each interval (specified by intervals). Should have one element less than the intervals parameter.
#'
#' @return returns a function that calculates the fractional efficiency for the specified diameter, see example.
#'
#' @references 
#' [1] Zaefferer, M.; Breiderhoff, B.; Naujoks, B.; Friese, M.; Stork, J.; Fischbach, A.; Flasch, O.; Bartz-Beielstein, T. Tuning Multi-objective Optimization Algorithms for Cyclone Dust Separators Proceedings of the 2014 Conference on Genetic and Evolutionary Computation, ACM, 2014, 1223-1230 \cr\cr
#' [2] Breiderhoff, B.; Bartz-Beielstein, T.; Naujoks, B.; Zaefferer, M.; Fischbach, A.; Flasch, O.; Friese, M.; Mersmann, O.; Stork, J.; Simulation and Optimization of Cyclone Dust Separators Proceedings 23. Workshop Computational Intelligence, 2013, 177-196
#'
#' @examples
#' ## Call directly
#' funCyclone(c(1.26,2.5))
#' ## create vectorized target funcion, vectorized, first objective only
#' ## Also: negated, since SPOT always does minimization.
#' tfunvecF1 <-function(x){-apply(x,1,funCyclone)[1,]}
#' tfunvecF1(matrix(c(1.26,2.5,1,2),2,2,byrow=TRUE))
#' ## optimize with spot
#' res <- spot(fun=tfunvecF1,lower=c(1,2),upper=c(2,3),
#'    control=list(modelControl=list(target="ei"),
#'    model=buildKriging,optimizer=optimLBFGSB,plots=TRUE)) 
#' ## best found solution ...
#' res$xbest
#' ## ... and its objective function value
#' res$ybest
#'
#' @export
###################################################################################
funCyclone <- function(x,deterministic=c(T,T,T),
					cyclone=list(Da=1.260,H=2.500,Dt=0.420,Ht=0.650,He=0.600,Be=0.200),
					fluid=list(Mu=1.85e-5,
							Ve=(50/36)/0.12,#Vp=5000,
							lambdag=1/200,
							Rhop=2000,Rhof=1.2,Croh=0.05),
					noiseLevel=list(Vp=0.1,Rhop=0.05),model="Barth-Muschelknautz",
					intervals=c(0,2,4,6,8,10,15,20,30)*1e-6,
					delta=c(0.0,0.02,0.03,0.05,0.1,0.3,0.3,0.2)
					){
	#--------------------------------------------------------------
	#Geometriedaten / geometric data
	#--------------------------------------------------------------
	#Da <- 1.260				      # Zyklondurchmesser [m]/ cyclone diameter [m] Parameter 8.1
	#H <- 2.500          		  # Zyklongesamthoehe  [m] / total cyclone height [m] Parameter 8.2
	#Dt <- 0.420          		  # Tauchrohrinnendurchmesser [m] / vortex finder diameter [m] Parameter 8.3
	#Ht <- 0.640          		  # Tauchrohreintauchtiefe / [m] Parameter 8.4
	#Xt <- 0            		  # Tauchrohrposition / position vortex finder
	#Phi <- 0             		# Tauchrohrposition / position vortex finder[degree] Parameter 8.5
	#Tt <- 0	           		  # Tauchrohrwanddicke [m] / thickness wall of vortex finder [m] Parameter 8.6
	#Epsilon <- 13.1340 		  # Konusneigungwinkel [grad] / angle of cyclone cone [degree] Parameter 8.7 old
	#Dx <- 0              		# Staubaustragsdurchmesser / discharge duct diameter [m] Parameter 8.8
	#He <- 0.600           		# Schlitzeinlaufhoehe [m] ]/inlet section height [m] Parameter 8.9
	#Be <- 0.200          		  # Schlitzeinlaufhoehe [m] / inlet section width [m] Parameter 8.10
	#Alpha <- 0         		# Einlaufwinkel [grad] / [degree] Parameter 8.11 
	#---------------------------------------------
	#Betriebsdaten / fluid parameters
	#-------------------------------------------
	#Mu <- 18.5*10^(-6)		  # Viskositaet / viscosity Zaehigkeit Luft
	######Vp <- 5000              # Volumenstrom in m/h / volumetric gas flow rate [m^3/h] !!!! Not used anymore, define Ve instead (inlet speed).
	#Ve <- (50/36)*0.12
	#Rhop <- 2000          # Partikeldichte in kg/m^3 / particle density [kg/m^3]
	#Rhof <- 1.2        	  # Gasdichte kg/m^3
	#Croh <- 0.05              # Rohgas Konzentration in [kg/m^3]
	#lambdag <- 1/200         # load free wall friction coefficient Eq.(2.7)
	#--------------------------------------------------------------
	if(is.na(x[[1]]))return(rep(NA,2))
	cyclone$Da<-x[1]
	if(length(x)>1)cyclone$H<-x[2]
	if(length(x)>2)cyclone$Dt<-x[3]
	if(length(x)>3)cyclone$Ht<-x[4]
	if(length(x)>4)cyclone$He<-x[5]
	if(length(x)>5)cyclone$Be<-x[6]
	
	fluid$Vp <- fluid$Ve * cyclone$He * cyclone$Be
	
	#--------------------------------------------------------------
	# Noise:
	fluid$Vp=if(deterministic[1]) fluid$Vp else fluid$Vp*(1-noiseLevel$Vp) + noiseLevel$Vp * fluid$Vp * 2 * runif(1)
	fluid$Rhop=if(deterministic[2]) fluid$Rhop else fluid$Rhop*(1-noiseLevel$Rhop) + noiseLevel$Rhop * fluid$Rhop * 2 * runif(1) #note: NEVER lower than Rhof, else sqrt(-) in model

	#--------------------------------------------------------------
	xmin <- intervals[-length(intervals)] 
	xmax <- intervals[-1]
	if(deterministic[3])
		xmean <- (xmax + xmin) / 2
	else
		xmean <- xmin + runif(length(xmin)) * (xmax - xmin)		
	#--------------------------------------------------------------
	if(model=="Barth-Muschelknautz"){
		C <- calculationBarthMuschelknautz(cyclone, fluid, xmean, delta)
		E <- C[1] #gesamtabscheidegrad
		Ew <- C[2] #abscheidegrad im wirbel
		PressureDrop <- C[3]
	}else{
		PressureDrop <- calculationMothes(1,cyclone, fluid)[2]
		E <- sum(matrix(unlist(lapply(xmean,calculationMothes,cyclone=cyclone,fluid=fluid)),2)[1,] * delta)
		Ew <- E #not included by Mothes
	}
	#--------------------------------------------------------------
	as.numeric(c(PressureDrop,-E,-Ew)) #first target: minimal pressure drop, second target: maximal col.eff. (inverted for minimization)
}
#funCyclone(c(1.26,2.5,.42,.65,.6,.2),fluid=list(Mu=1.85e-5,Ve=(50/36)/0.12,lambdag=1/200,Rhop=2000,Rhof=1.2,Croh=0.05))
#par(mfrow=c(2,1))
#fun <- function(x) funCyclone(c(1.260,2.500,.420,x))
#f <- function(x) sapply(x,fun)[2,]
#curve(f,from=0,to=1,main="Loeffler")
#fun2 <- function(x) funCyclone(c(1.260,2.500,.420,x),model="Mothes")
#f2 <- function(x) sapply(x,fun2)[2,]
#curve(f2,from=0,to=1,main="Mothes")
