!**************************************************************************************************
! Module Energy Use for Residential Appliances 
!
! Original Indian Version
! Author:   Bas van Ruijven
! Date:     May 2008
!
! Global Version By: Vassilis Daioglou
! Date:     May 2010
!**************************************************************************************************

!**************************************************************************************************
MODULE Appliance;
BEGIN

DOUBLE	af;		!new domain, annuity factor
DOUBLE	ctax;
DOUBLE  ppsf;		!Double used to determine to which extent the maximum energy savings are reached at a certain carbon tax level (declared as a timestep to allow automatic interpolation)

!**************************************************************************************************

#INCLUDE ../Data/DatAppliance_TIMER.m	! Read datafiles

! 1=fan, 2= Air Cooler, 3=Air Conditioner, 4=Refrigerator, 5=Microwave, 6=Washing Machine, 7=Clothes Dryer, 8=Dish Washer, 9=TV, 10=VCR/DVD, 11=PC/Other 
! Four major end-uses:
! 1) Cooling (1-3)
! 2) Food storage and processing (4-5)
! 3) Washing & Cleaning (6-8)
! 4) Entertainment (9-11)   

IMPORT REAL	
	CarbonTax[NR27](t),			! From cookingwater.m
	PriceSecFuel[NR27](t),	
	Fuel_Subsidy[NR27,TURQ,NECN](t),
    	convergence(t),
   	convergence2(t),
   	PCOpcT_ppp[NR27,TURQ](t),   	!HHExp per capita for TURQ in pppUSD2005
	CDR1[NR27,TURQ](t),
    	Households[NR27,TURQ](t),
    	FloorSpace[NR27,TURQ](t),   	! Floorspace per capita
    	HHFloor[NR27,TURQ](t),      	! Total Floorspace
    	Electrification[NR27,TURQ](t),  	! From drivers/electrification.m
    	CDD[NR27,NMT](t),      	! From CoolHeat.m
	InsulationCoolingCor[NR27,TURQ](t),	! Actual reduction in cooling demand (%) applied due to insulation of dwellings
    	Urbanization[NR27](t),
    	HHExpFac_ppp[NR27,TURQ](t), 	! From drivers.m
    	POP[NRCT](t),
    	POP_q[NR27,TURQ](t),
    	OutageFactor[NR27](t),
    	EffSmoothFactor(t),
    	PayBackTimeTC[NS];

REAL	AnnuityFactor[NR27,TURQ](t),
!	FuelPrice_Final[NR27,TURQ](t),
	CDR_extra	= 0,
	coe_kWh1[NR27](t),
	coe_kWh[NR27](t),
	EffCooling[NR27](t),
    	SatPrice[NR27](t),
    	Saturation[NR27,TURQ,7](t),
    	Availability[NR27,TURQ,3](t),  
    	Diffusion[NR27,TURQ,3](t),  	! Diffusion in TURQ for Fans, air coolers and air conditioners
    	Diffusion1[NR27,TUR,11](t), 	! Diffusion in TUR for rest
    	Diffusion3[NR27,TURQ,11](t),    	! Diffusion in Q for rest
    	DiffusionTOT[NR27,TURQ,11](t),
    	AdjDiffusion[NR27,TURQ,11](t),
    	AvgDiffusion[NR27,TUR,11](t),
    	ApplCapReq[NR27,TURQ,11](t),
    	ApplCap_ini[NR27,TURQ,11](t),
    	ApplCapTot[NR27,TURQ,11](t),
    	ApplCapTOTAL[NR27,TURQ,11](t),
    	ApplCapNew[NR27,TURQ,11](t),
    	ApplCapDepr[NR27,TURQ,11](t),
    	ApplNewCapSmth[NR27,11](t),
    	ApplCapTOTALphh[NR27,TURQ,11](t),
    	CoolApplCapphh[NR27,TURQ](t),
    	ClimateMaxSaturation[NR27](t),
    	UEC_min[NRTim,11](t),
    	UEC[NR27,TURQ,11](t),
	UEC_base1[NR27,11](t),
    	UEC_base[NR27,11](t),
	UEC_unit_new1[NR27,TURQ,11](t),		! kWh/L after coe effects accounted for, baseline coe
	UEC_unit_new2[NR27,TURQ,11](t),		! kWh/L after coe effects accounted for, changed coe
	UEC_diff[NR27,TURQ,11](t),		! difference
	UEC_unit_new[NR27,TURQ,11](t),
	UEC_new[NR27,TURQ,11](t),		! kWh/unit after coe effects accounted for
	StandByFac[NR27,11](t),		! Multiplier for energy savings assuming appliances have no stand-by
	EffApplFac[NR27,11](t),		! Multiplier for energy savings assuming changes in appliance use (more efficient use)
    	ApplEnReq[NR27,TURQ,11](t),     		! Energy Requirement of all appliences, tmin    
    	ApplEnUse_ini[NR27,TURQ,11](t),     
    	ApplEnUseTot[NR27,TURQ,11](t),      	! Energy requirement of all appliances, all other t
    	ApplEnUse_new[NR27,TURQ,11](t),    	! Energy use of all new (marginal) appliances at given t
    	ApplEnUse_Depr[NR27,TURQ,11](t),        	! Energy use of all depreciated appliances at given t
    	ApplEnUseTot_Cor[NR27,TURQ,11](t),      	! Energy requirement of all appliances, corrected for improvements in insulation (affects air conditioners/coorlers)
    	ApplEnUse_total[NR27,TURQ,11](t),       	! Energy requirement of all appliances
    	ApplUEC_aggr[NR27,TURQ,11](t),      	! Aggregate UEC when energy standards are introduced

    	Other[NR27,TURQ](t),
    	OtherShare[NR27,TURQ](t),   
    	ApplEnUseTOTALpc[NRC](t),
   ! 	ApplEnUse[NR27,TURQ,11](t),
    	ApplEnUse_tuss[NRC,TURQ,12](t),	!12 = other
    	ApplEnUse1[NRC,TURQ](t),
    	ApplEnAggr[NRC,TUR](t);
	   

REAL var1[NRC,11](t);		!These have to be set up in order to accomodate changing of UEC_beta, otherwise curve is discontinuous.
REAL var2[NRC,11](t);
REAL
ApplEnUse3[NRC,TURQ](t), !kWh/cap for cooling appliances
ApplEnUse4[NRC,TURQ](t); !kWh/cap for rest of appliances

EXPORT REAL 
ApplEnAggr2[NRC,TURQ,NECN](t), !GJ For ALL Appliances
ApplEnAggr3[NRC,TURQ,NECN](t), !GJ for cooling appliances
ApplEnAggr4[NRC,TURQ,NECN](t), !GJ for rest of appliances
UEC_Improv[NRC](t),		! Aggregate improvement in UEC compared to tscen
Appl_EnExp[NRC,TURQ](t); 

REAL
ApplEnAggr3pc[NR27,TURQ,NECN](t),
ApplEnAggr4pc[NR27,TURQ,NECN](t);

INTEGER ApplTechLT          		= 15;
EXPORT REAL    	
	ApplEnUseTOTAL[NR27,NECN](t),
	FuelPrice_Final[NR27,TURQ](t),
	ApplEnUse[NR27,TURQ,11](t);
	
!**************************************************************************************************
! Economic parameters
AnnuityFactor[R,i] 	= SWITCH( CDR1[R,i] > EPS ?
				(CDR1[R,i] + CDR_extra)/(1-(1+(CDR1[R,i] + CDR_extra))**-(ApplTechLT))
			  ELSE 0), R = 1 to NRTim, i = 1 to TURQ;

coe_kWh1[R]		= PriceSecFuel[R] * 1.252083 * 0.0036, R = 1 to NRTim;	! 2005$/kWh, baseline
coe_kWh[R]		= PriceSecFuel[R] * 1.252083 * 0.0036, R = 1 to NRTim;		! 2005$/kWh, future

FuelPrice_Final[R,i]	= PriceSecFuel[R] - Fuel_Subsidy[R,i,8], R = 1 to 26, i = 1 to TURQ;
FuelPrice_Final[NRTim,i]	= 0,i = 1 to TURQ;

!**************************************************************************************************
! Calcultations for Air Conditioning from Morna/Detlef
ClimateMaxSaturation[R] = 1.00 - 0.949 * EXP(-0.00187 * CDD[R,NMT]), R = 1 to NRTim;

!**************************************************************************************************
! Diffusion of Fans as function of Floorspace PC (Variation on McNeil/Letschert 2007)
! Historic saturation level is function of price development (Based on 1993/2002 difference for India)
SatPrice[R]     	= 1, R = 1 to 2;
SatPrice[R]     	= SatpriceREMI, R = 3 to 10;
SatPrice[R]     	= 1, R = 11 to 12;
SatPrice[R]     	= SatPriceREMI, R = 13 to 22;
Satprice[R]     	= 1, R = 23 to 24;
Satprice[R]     	= SatPriceREMI, R = 25 to NRTim;

Saturation[R,i,1]   	= SatPrice[R] * ClimateMaxSaturation[R], R = 1 to NRTim, i=1 to TURQ;           		!For fan only
Saturation[R,i,j]   	= AirCoolFrac[R] * SatPrice[R] * ClimateMaxSaturation[R], R = 1 to NRTim, i=1 to TURQ, j=2;      	!For air cooler only
Saturation[R,i,j]   	= (1-AirCoolFrac[R]) * SatPrice[R] * ClimateMaxSaturation[R], R = 1 to NRTim, i=1 to TURQ, j=3;  	!For air conditioner only   

! Fitted diffusion levels of Fans, Air Coolers and Air Conditioners (TURQ)
Availability[R,i,j]    	= EXP(-Beta[R,i,j]*EXP(-(Gamma[R,i,j]/1000)*PCOpcT_ppp[R,i])), R = 1 to NRTim, i=1 to TUR, j=1 to 3;
Availability[R,i,j]    	= EXP(-Beta[R,2,j]*EXP(-(Gamma[R,2,j]/1000)*PCOpcT_ppp[R,i])), R = 1 to NRTim, i=4 to 8, j=1 to 3;
Availability[R,i,j]    	= EXP(-Beta[R,3,j]*EXP(-(Gamma[R,3,j]/1000)*PCOpcT_ppp[R,i])), R = 1 to NRTim, i=9 to TURQ, j=1 to 3;

! Fitted diffusion levels of Fans, Air Coolers and Air Conditioners (TURQ)
Diffusion[R,i,j]    	= Saturation[R,i,j]*Availability[R,i,j], R = 1 to NRTim, i=1 to TURQ, j=1 to 3;

Diffusion1[R,i,j]   	= 0, R = 1 to NRTim, i = 1 to TUR, j = 1 to 3;
! DIFFUSION FOR REST OF APPLIANCES FOLLOW 3 BASIC PATTERNS
! In all cases Total averaged out from Urban and Rural.
! 1. Constant Growth: Clothes Dryer, Dish Washer
! 2. An 'income delay', income at which diffusion starts depends on time: Microwave, VCR, PC
! 3. 'Surface', rate of diffusion and saturation point vary: Refrigerator, Washing Machine, TV
! Appliances with income delay: !@1970 diffusion starts at 10000$(2005)ppp, !@2000 diffusion starts at 700$(2005)ppp.

! FOOD PREPARATION AND STORAGE APPLIANCES
Diffusion1[R,i,4]   	= phi1[R,4]*EXP(-phi2[R,i,4]*EXP(-(Phi3[R,i,4]/1000)*PCOpcT_ppp[R,i])), R = 1 to NRTim, i = 2 to TUR;   			! Refrigerators
Diffusion1[R,i,5]   	= phi1[R,5] * EXP(-phi2[R,i,5] * EXP(-(phi3[R,i,5]/1000) * (PCOpcT_ppp[R,i] - IncomeDelay))), R = 1 to NRTim, i = 2 to TUR; 	! Microwave
! CLEANING APPLIANCES
Diffusion1[R,i,6]   	= phi1[R,6]*EXP(-phi2[R,i,6]*EXP(-(Phi3[R,i,6]/1000)*PCOpcT_ppp[R,i])), R = 1 to NRTim, i = 2 to TUR;      		! Washing Machine
Diffusion1[R,i,7]   	= phi1[R,7] * EXP(-phi2[R,i,7] * EXP(-(phi3[R,i,7]/1000) * PCOpcT_ppp[R,i])), R = 1 to NRTim, i = 2 to TUR; 		! Clothes Dryer
Diffusion1[R,i,8]   	= phi1[R,8] * EXP(-phi2[R,i,8] * EXP(-(phi3[R,i,8]/1000)* PCOpcT_ppp[R,i])), R = 1 to NRTim, i = 2 to TUR;  		! Dish Washer
! ENTERTAINMENT APPLIANCES
Diffusion1[R,i,9]   	= phi1[R,9]*EXP(-phi2[R,i,9]*EXP(-(phi3[R,i,9]/1000)*PCOpcT_ppp[R,i])), R = 1 to NRTim, i = 2 to TUR;  			! Television
Diffusion1[R,i,10]  	= phi1[R,10] * EXP(-phi2[R,i,10]* EXP(-(phi3[R,i,10]/1000) * (PCOpcT_ppp[R,i] - IncomeDelay))), R = 1 to NRTim, i = 2 to TUR;   	! DVD/VCR
Diffusion1[R,i,11]  	= phi1[R,11] * EXP(-phi2[R,i,11] * EXP(-(phi3[R,i,11]/1000) * (PCOpcT_ppp[R,i] - IncomeDelay))), R = 1 to NRTim, i = 2 to TUR;  	! PC
!TOTAL allocation
Diffusion1[R,1,j]   	= (Diffusion1[R,2,j] * Urbanization[R]) + (Diffusion1[R,3,j] * (1-Urbanization[R])), R = 1 to NRTim, j = 4 to 11;   

!Quintile Allocation    
Diffusion3[R,i,j] 	= SWITCH( HHExpFac_ppp[R,i] > EPS ?
	            		(Appliance_a[j]*LOG(ABS(HHExpFac_ppp[R,i]))+1) * Diffusion1[R,2,j]
             		  ELSE 0), R = 1 to NRTim, i = 4 to 8, j = 4 to 11;
Diffusion3[R,i,j] 	= SWITCH( HHExpFac_ppp[R,i] > EPS ? 
            		 	(Appliance_a[j]*LOG(ABS(HHExpFac_ppp[R,i]))+1) * Diffusion1[R,3,j] 
            		  ELSE 0), R = 1 to NRTim, i = 9 to TURQ, j = 4 to 11;

DiffusionTOT[R,i,j] 	= Diffusion[R,i,j], R = 1 to NRTim, i = 1 to TURQ, j = 1 to 3;  ! Only for fans and AC's
DiffusionTOT[R,1,j] 	= Diffusion1[R,1,j], R = 1 to NRTim, j = 4 to 11;
DiffusionTOT[R,i,j] 	= Diffusion1[R,i,j], R = 1 to NRTim, i = 2 to 3, j = 4 to 11;
DiffusionTOT[R,i,j] 	= Diffusion3[R,i,j], R = 1 to NRTim, i = 4 to TURQ, j = 4 to 11;

! Adjust for Electrification
AdjDiffusion[R,i,j] 	= DiffusionTOT[R,i,j] * Electrification[R,i], R = 1 TO NRTim, i=1 to TURQ, j=1 to 11;  !It has to be electrified

AvgDiffusion[R,1,j] 	= SWITCH( households[R,1] > EPS ?
           		 	LSUM(i = 4 to TURQ,AdjDiffusion[R,i,j]*households[R,i]) / households[R,1]
           		  ELSE 0), R = 1 to NRTim, j=1 to 11;
AvgDiffusion[R,2,j] 	= SWITCH( households[R,2] > EPS ?
            			LSUM(i=4 to TURQ-5,AdjDiffusion[R,i,j]*households[R,i]) / households[R,2]
           		  ELSE 0), R = 1 to NRTim, j=1 to 11;
AvgDiffusion[R,3,j] 	= SWITCH( households[R,3] > EPS ?
            			LSUM(i=9 to TURQ,AdjDiffusion[R,i,j]*households[R,i]) / households[R,3]
           		  ELSE 0), R = 1 to NRTim, j=1 to 11;
!**************************************************************************************************
! Vintage capital stock model

ApplCapReq[R,i,j]       = households[R,i]*AdjDiffusion[R,i,j]*(0.04*HHFloor[R,i]), R = 1 to NRTim, i=1 to TURQ, j=1; !Fans only
ApplCapReq[R,i,j]       = households[R,i]*AdjDiffusion[R,i,j], R = 1 to NRTim, i=1 to TURQ, j=2 to 11;              !Rest
ApplCap_ini[R,i,j]  	= LAST(ApplCap_ini[R,i,j],ApplCapReq[R,i,j]), R = 1 to NRTim, i=1 to TURQ, j=1 to 11;
ApplCapTot[R,i,j]       = LAST(ApplCapTot[R,i,j],ApplCap_ini[R,i,j]) - ApplCapDepr[R,i,j] + ApplCapNew[R,i,j], R = 1 to NRTim, i=1 to TURQ, j=1 to 11;
ApplCapNew[R,i,j]   	= MAX(0,ApplCapReq[R,i,j] + ApplCapDepr[R,i,j] - LAST(ApplCapTot[R,i,j],ApplCap_ini[R,i,j])), R = 1 to NRTim, i=1 to TURQ, j=1 to 11;
ApplCapDepr[R,i,j]      = SWITCH( (T - T.MIN) < (ApplTechLT-2) AND ApplTechLT>EPS ? ApplCap_ini[R,i,j]/ApplTechLT,
                  		(T - T.MIN) < (ApplTechLT+2) ? LSUM(k = 0 to 4, 1/5*NLAST(ApplCapNew[R,i,j],ApplTechLT + (k-2),0)) + MAX(0,(t.min + ApplTechLT + 2 - t)/5) * ApplCap_ini[R,i,j]/ApplTechLT
             	 	  ELSE    LSUM(k = 0 to 4, 1/5*NLAST(ApplCapNew[R,i,j],ApplTechLT + (k-2),0))), R = 1 to NRTim, i=1 to TURQ, j=1 to 11;
ApplNewCapSmth[R,j] 	= LSUM(i=4 to TURQ,LAVG(k=1 to 5,NLAST(ApplCapNew[R,i,j],k,ApplCapNew[R,i,j]))), R = 1 to NRTim, j=1 to 11;
ApplCapTOTAL[R,1,j] 	= LSUM(i=4 to TURQ,ApplCapTot[R,i,j]), R = 1 to NRTim, j=1 to 11;
ApplCapTOTAL[R,2,j] 	= LSUM(i=4 to TURQ-5,ApplCapTot[R,i,j]), R = 1 to NRTim, j=1 to 11;
ApplCapTOTAL[R,3,j] 	= LSUM(i=9 to TURQ,ApplCapTot[R,i,j]), R = 1 to NRTim, j=1 to 11;
ApplCapTOTAL[R,i,j] 	= ApplCapTot[R,i,j], R = 1 to NRTim, i=4 to TURQ, j=1 to 11;
ApplCapTOTALphh[R,i,j]  = SWITCH( households[R,i] > EPS ?
            			ApplCapTOTAL[R,i,j]/households[R,i]
            		  ELSE 0), R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11;
CoolApplCapphh[R,i] 	= LSUM ( j = 2 to 3, ApplCapTOTALphh[R,i,j]), R = 1 to NRTim, i = 1 to TURQ;    ! # of air coolers and conditioners per household, for calibration purposes.
!**************************************************************************************************
! Energy Use (in GJ/yr), UEC in kWh
! Unit Energy Consumption for Fan, Air Cooler and Air Conditioner

EffCooling[R]		= SWITCH( t > tscen ?
				MAX(Effcooling_fut[R],EffCooling_ctax[R](CarbonTax[R]))
			  ELSE EffCooling1[R]), R = 1 to NRTim;



! Unit energy consumption for rest of appliances, autonomous decline over time (determine dby UEC_beta), described in normalised units:
! UEC_beta is scenario dependent.
! 4. Refrigerator: kWh/L (total volume)
! 5. Microwave: kWh/unit
! 6. Washing Machine: kWh/L (washing volume)
! 7. Clothes Dryer: kWh/kg (load)
! 8. Dish Washer: kWh/washing cycles 	(not cycles per year, but cycles per wash!!!)
! 9/10/11. TV/DVD-VCR/PC: kWh/unit

UEC_min[R,i]		= SWITCH( t <= tscen ? UEC_min_in[R,i]
			  ELSE (Convergence2 * UEC_min_scen[R,i]) + ((1-convergence2) * UEC_min_in[R,i])), R = 1 TO NRTim, i = 4 TO 11;

UEC_base1[R,i]		= UEC_alpha[R,i] * UEC_beta[R,i]**(t-1970) + UEC_min[R,i], R = 1 to NRTim, i = 4 to 11;

var1[R,i]		= UEC_base1[R,i](tscen) - UEC_min[R,i], R = 1 to NRTim, i = 4 to 11;
var2[R,i]		= SWITCH( t > tscen ?
				(UEC_beta_fut[R,i]*ApplCtaxEff(CarbonTax[R]))**(t-tscen)
			  ELSE 0), R = 1 to NRTim, i = 4 to 11;

UEC_base[R,i]		= SWITCH ( t <= tscen ?
				UEC_base1[R,i]
			  ELSE var1[R,i] * var2[R,i] + UEC_min[R,i]), R = 1 to NRTim, i = 4 to 11;
! Deviation from baseline due to coe
	! Below variables irrelevant for cooling appliances
	UEC_unit_new1[R,i,j] 	= 0.0, R = 1 TO NRTim, i = 1 TO TURQ, j = 1 TO 3;
	UEC_unit_new2[R,i,j] 	= 0.0, R = 1 TO NRTim, i = 1 TO TURQ, j = 1 TO 3;
	UEC_diff[R,i,j] 	= 0, R = 1 TO NRTim, i = 1 TO TURQ, j = 1 TO 3;
	UEC_unit_new[R,i,j] 	= 0, R = 1 TO NRTim, i = 1 TO TURQ, j = 1 TO 3;
	UEC_new[R,i,j] 	= 0, R = 1 TO NRTim, i = 1 TO TURQ, j = 1 TO 3;

UEC_unit_new1[R,i,j] 	= elast_coeff[1,j](AnnuityFactor[R,i]) * LOG(coe_kWh1[R]) + elast_coeff[2,j](AnnuityFactor[R,i]), R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11; 	! kWh/unit, at baseline coe
UEC_unit_new2[R,i,j]	= elast_coeff[1,j](AnnuityFactor[R,i]) * LOG(coe_kWh[R]) + elast_coeff[2,j](AnnuityFactor[R,i]), R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11;	! kWh/unit at new coe
UEC_diff[R,i,j]		= UEC_unit_new1[R,i,j] - UEC_unit_new2[R,i,j], R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11; ! this is done because we want the curve to go through UECbase in t=2010
UEC_unit_new[R,i,j]	= SWITCH( UEC_base[R,j] - UEC_diff[R,i,j] > UEC_min[R,j] ?
				UEC_base[R,j] - UEC_diff[R,i,j]
			  ELSE UEC_min[R,j]), R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11;
UEC_new[R,i,j]		= UEC_unit_new[R,i,j] * UEC_unit[R,j], R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11;	! kWh/yr, of new refrigerators purchased

! Standyby and efficient appliance factors = 1 unless overridden from scenario
StandByFac[R,i]		= SWITCH (FlagStandByMode = 1 AND t > tscen ?
				StandByFac_in[R,i]
			ELSE 1), R = 1 TO NRTim, i = 1 TO 11;	

EffApplFac[R, i]		= SWITCH (FlagEffAppl = 1 AND t > tscen ?
				EffApplFac_in[R,i]
			ELSE 1), R = 1 TO NRTim, i = 1 TO 11;

! Final Unit Energy Consumption values:
UEC[R,i,j]      	= 0.0401 * CDD[R,NMT] + 22.282, R = 1 to NRTim, i = 1 to TURQ, j = 1;   ! Fans

UEC[R,i,j]      	= SWITCH(PCOpcT_ppp[R,i] > EPS ?
               			MAX(2.5/EffCooling[R]*CDD[R,NMT]*(0.6053*LOG(PCOpcT_ppp[R,i])-3.1897),400)*300/2160
               		  ELSE 0), R = 1 to NRTim, i=1 to TURQ, j=2;  !Air cooler
UEC[R,i,j]     		= SWITCH(PCOpcT_ppp[R,i] > EPS ?
                		MAX(2.5/EffCooling[R]*CDD[R,NMT]*(0.6053*LOG(PCOpcT_ppp[R,i])-3.1897),400)
               	 	  ELSE 0), R = 1 to NRTim, i=1 to TURQ, j=3;  !Air Conditioner

UEC[R,i,j]		= UEC_new[R,i,j] * StandByFac[R,j] * EffApplFac[R,j], R = 1 to NRTim, i = 1 to TURQ, j = 4 to 11;

! Determine improvements in efficiency, to be used by service sector
UEC_Improv[R]	= SWITCH( t < tscen ? 1
		ELSE	((UEC[R,2,4] + UEC[R,2,6])/2) / ((UEC[R,2,4](tscen) + UEC[R,2,6](tscen))/2)), R = 1 TO NRC;

!**************************************************************************************************
! Other unnacounted for appliance energy use, result from calibration to IEA data (30 years of Energy Use in IEA countries).
! Related to income.
! 'Other' Calculated on a per capita basis, since IEA data is on a per capita basis.
! The coefficients of the "other" are based on regions which seem to have very high, when compared to default REMG, miscellaneous applies (Canada, USA, Japan)
! and regions with low miscelaneous appliances (Western Europe, australia, and ROW by default).
Other[R,i]      	= SWITCH( PCOpcT_ppp[R,i] > EPS ?	
				SWITCH( Other_coeff[R,1] * LOG(ABS(PCOpcT_ppp[R,i])) - Other_coeff[R,2] > 0 ?
            				Other_coeff[R,1] * LOG(ABS(PCOpcT_ppp[R,i])) - Other_coeff[R,2]
            			ELSE 0)
			  ELSE 0), R = 1 to NRTim, i = 1 to TURQ;     ! Extra kwh/cap per region, FOR TOTAL

! Since UEC has been 	calculated on a marginal basis, and energy stock model is necessary to get aggregate UEC's over appliance lifetime
! NOT for Fans, Air Coolers and Air Conditioners which are climate based
ApplEnReq[R,i,j]    	= ApplCapReq[R,i,j] * UEC[R,i,j], R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11;
ApplEnUse_ini[R,i,j]    = LAST(ApplEnUse_ini[R,i,j],ApplEnReq[R,i,j]), R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11;
ApplEnUseTot[R,i,j] 	= (LAST(ApplEnUseTot[R,i,j],ApplEnUse_ini[R,i,j]) - ApplEnUse_Depr[R,i,j] + ApplEnUse_new[R,i,j]), R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11;
ApplEnUse_new[R,i,j]    = ApplCapNew[R,i,j] * UEC[R,i,j], R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11;
ApplEnUse_Depr[R,i,j]   = SWITCH( (t - t.min) < (ApplTechLT-2) ?   
					ApplCapDepr[R,i,j] * UEC_ini[R,i,j],
        			(t - t.min) < (ApplTechLT) + 2 ?
            				LSUM(k = 0 to 4, 1/5*NLAST(ApplEnUse_new[R,i,j],ApplTechLT + (k-2),0)) + MAX(0,(t.min + ApplTechLT + 2 - t)/5)*ApplCapDepr[R,i,j] * UEC_ini[R,i,j],
       		 	  ELSE LSUM(k = 0 to 4, 1/5*NLAST(ApplEnUse_new[R,i,j],ApplTechLT + (k-2),0))), R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11;

! Make corrections of cooling appliances when there is insulation (via renovation). This has to be done separately in order to not interfere with the energy stock model
ApplEnUseTot_Cor[R,i,1]	= ApplEnUseTot[R,i,1], R = 1 TO NRTim, i =  1 TO TURQ;
ApplEnUseTot_Cor[R,i,j]	= InsulationCoolingCor[R,i] * ApplEnUseTot[R,i,j], R = 1 TO NRTim, i =  1 TO TURQ, j = 2 TO 3;;
ApplEnUseTot_Cor[R,i,j]	= ApplEnUseTot[R,i,j], R = 1 TO NRTim, i =  1 TO TURQ, j = 4 TO 11;

ApplEnUse_total[R,1,j]  = LSUM(i = 4 to TURQ, ApplEnUseTot_Cor[R,i,j]), R = 1 to NRTim, j = 1 to 11;    ! kWh, OUTAGE FACTOR NOT ACCOUNTED FOR YET
ApplEnUse_total[R,2,j]  = LSUM(i = 4 to TURQ-5, ApplEnUseTot_Cor[R,i,j]), R = 1 to NRTim, j = 1 to 11;
ApplEnUse_total[R,3,j]  = LSUM(i = 9 to TURQ, ApplEnUseTot_Cor[R,i,j]), R = 1 to NRTim, j = 1 to 11;
ApplEnUSe_total[R,i,j]  = ApplEnUSeTot_Cor[R,i,j], R = 1 to NRTim, i = 4 to TURQ, j = 1 to 11;

ApplEnUse[R,i,j]    	= (ApplEnUse_total[R,i,j] * OutageFactor[R])/POP_q[R,i], R = 1 to 26, i = 1 to TURQ, j = 1 to 11;   ! kWh/cap
ApplEnUse[NRTim,i,j]    = 0, i = 1 to TURQ, j = 1 to 11;

ApplUEC_aggr[R,i,j] 	= SWITCH ( ApplCapTotal[R,i,j] > EPS ? 
				ApplEnUse_total[R,i,j]/ApplCapTOTAL[R,i,j]
		 	  ELSE 0), R = 1 to 26, i = 1 to TURQ, j = 1 to 11;
ApplUEC_aggr[NRTim,i,j] = 0, i = 1 to TURQ, j = 1 to 11;

! Electricity use for ALL appliances NOTE IT IS PER CAPITA!!!!!!!!!!!
ApplEnUse1[R,i] 	= (LSUM(j=1 to 11,ApplEnUse[R,i,j])) + Other[R,i], R = 1 to NRTim, i = 1 to TURQ;   ! kWh/cap, all appliances

ApplEnUse3[R,i]		= (LSUM(j=1 to 3,ApplEnUse[R,i,j])), R = 1 to NRTim, i = 1 to TURQ;
ApplEnUse4[R,i]		= (LSUM(j=4 to 11,ApplEnUse[R,i,j])) + Other[R,i], R = 1 to NRTim, i = 1 to TURQ;

! output for TUSS
ApplEnUse_tuss[R,i,j]	= ApplEnUse_total[R,i,j] * OutageFactor[R] * MJperkWh/1000, 	R = 1 to NRTim, i = 1 to TURQ, j = 1 to 11; 	! GJ
ApplEnUse_tuss[R,1,12]  = LSUM(i=4 to TURQ, Other[R,i] * POP_q[R,i] * OutageFactor[R]) * MJperkWh/1000, R = 1 to NRTim;			! GJ
ApplEnUse_tuss[R,2,12]  = LSUM(i=4 to TURQ-5, Other[R,i] * POP_q[R,i] * OutageFactor[R]) * MJperkWh/1000, R = 1 to NRTim;		! GJ
ApplEnUse_tuss[R,3,12]  = LSUM(i=9 to TURQ, Other[R,i] * POP_q[R,i] * OutageFactor[R]) * MJperkWh/1000, R = 1 to NRTim;			! GJ
ApplEnUse_tuss[R,i,12]  = Other[R,i] * POP_q[R,i] * OutageFactor[R] * MJperkWh/1000, R = 1 to NRTim, i = 4 to TURQ;		! GJ

OtherShare[R,i] 	= SWITCH( ApplEnUse1[R,i] > 0 AND Other[R,i] > 0 ?
            			(Other[R,i]/ApplEnUse1[R,i])*100
           		  ELSE 0), R = 1 to NRTim, i = 1 to TURQ;
! GJ for appliances TUR
ApplEnAggr[R,1] 	= LSUM(i=4 to TURQ,ApplEnUse1[R,i] * POP_q[R,i]) * MJperkWh/1000, R = 1 to NRTim;       ! GJ
ApplEnAggr[R,2] 	= LSUM(i=4 to TURQ-5,ApplEnUse1[R,i] * POP_q[R,i]) * MJperkWh/1000, R = 1 to NRTim; 	! GJ
ApplEnAggr[R,3] 	= LSUM(i=9 to TURQ,ApplEnUse1[R,i] * POP_q[R,i]) * MJperkWh/1000, R = 1 to NRTim;       ! GJ

!ALL appliances
ApplEnAggr2[R,i,j]	= 0, R = 1 to NRTim, i = 1 to TURQ, j = 1 to 7;
ApplEnAggr2[R,i,NECN]	= ApplEnUse1[R,i] * POP_q[R,i] * MJperkWh/1000, R = 1 to 26, i = 4 to TURQ;			! GJ

ApplEnAggr2[R,1,NECN]	= LSUM(i=4 to TURQ, ApplEnAggr2[R,i,NECN]), R = 1 to 26;
ApplEnAggr2[R,2,NECN]	= LSUM(i=4 to TURQ-5, ApplEnAggr2[R,i,NECN]), R = 1 to 26;
ApplEnAggr2[R,3,NECN]	= LSUM(i =9 to TURQ, ApplEnAggr2[R,i,NECN]), R = 1 to 26;
ApplEnAggr2[27,i,NECN]	= 0, i = 1 to TURQ;

! Cooling appliances
ApplEnAggr3[R,i,j]	= 0.0, R = 1 to NRC2, i = 1 to TURQ, j = 1 to 7;
ApplEnAggr3[R,i,NECN]	= ApplEnUse3[R,i] * POP_q[R,i] * MJperkWh/1000, R = 1 to NRC2, i = 4 to TURQ;			! GJ

ApplEnAggr3[R,1,NECN]	= LSUM(i=4 to TURQ, ApplEnAggr3[R,i,NECN]), R = 1 to NRC2;
ApplEnAggr3[R,2,NECN]	= LSUM(i=4 to TURQ-5, ApplEnAggr3[R,i,NECN]), R = 1 to NRC2;
ApplEnAggr3[R,3,NECN]	= LSUM(i =9 to TURQ, ApplEnAggr3[R,i,NECN]), R = 1 to NRC2;
ApplEnAggr3[NRC,i,j]	= LSUM(R =1 to NRC2, ApplEnAggr3[R,i,j]), i = 1 to TURQ, j = 1 TO NECN;

ApplEnAggr3pc[R,i,j]	= SWITCH( POP_q[R,i] > EPS ?
				ApplEnAggr3[R,i,j]/POP_q[R,i]
			  ELSE 0), R = 1 to NRTim, i  = 1 to TURQ, j = 1 to NECN;

! Rest of appliances
ApplEnAggr4[R,i,j]	= 0, R = 1 to NRTim, i = 1 to TURQ, j = 1 to 7;
ApplEnAggr4[R,i,NECN]	= ApplEnUse4[R,i] * POP_q[R,i] * MJperkWh/1000, R = 1 to 26, i = 4 to TURQ;			! GJ

ApplEnAggr4[R,1,NECN]	= LSUM(i=4 to TURQ, ApplEnAggr4[R,i,NECN]), R = 1 to 26;
ApplEnAggr4[R,2,NECN]	= LSUM(i=4 to TURQ-5, ApplEnAggr4[R,i,NECN]), R = 1 to 26;
ApplEnAggr4[R,3,NECN]	= LSUM(i =9 to TURQ, ApplEnAggr4[R,i,NECN]), R = 1 to 26;
ApplEnAggr4[27,i,NECN]	= 0, i = 1 to TURQ;

ApplEnAggr4pc[R,i,j]	= SWITCH( POP_q[R,i] > EPS ?
				ApplEnAggr4[R,i,j]/POP_q[R,i]
			  ELSE 0), R = 1 to NRTim, i = 1 to TURQ, j = 1 to NECN;
! In NECN for export
ApplEnUseTOTAL[R,EC]   	= 0, R = 1 to NRTim, EC=1 to NECN-1;        ! GJ
ApplEnUseTOTAL[R,NECN]  = ApplEnAggr[R,1], R = 1 to NRTim;      ! GJ

! Per cap for comparison
ApplEnUseTOTALpc[R] 	= SWITCH( POP_q[R,1] > 0 ?          ! GJ/cap
            			ApplEnAggr[R,1]/POP_q[R,1]
            		  ELSE 0), R = 1 to NRTim;

! Expenditures for Appliance Energy
Appl_EnExp[R,i]	= FuelPrice_Final[R,i] * ApplEnAggr2[R,i,8], R = 1 to 26, i = 1 to TURQ;
Appl_EnExp[NRTim,i]	= 0, i = 1 to TURQ;

END;