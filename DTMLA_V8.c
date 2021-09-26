#include <windows.h>
#include <math.h>
#include <string.h>
#include "usersurf.h"
#define ITMAX 1.0e+12 
#define Center_location_resolution 1.0e-12
#define SnellSign(U) U*U/(Dz*Dz + U * U) - (FD->n1)*(FD->n1)*(D - U)*(D - U) / (RT*RT + (D - U)*(D - U))

/*

Original code written by James Sutter
Sept 21, 2001

Modified by Kenneth Moore to fix bug in plano case, and add paraxial power due to r^2 term
Feb 24, 2005

Modified by Mitch Harter to use FIXED_DATA5 structure
August 12, 2016

 Modified by Omer Korech to fix some bug in plano case, and add a feature.
 The added feature: Tilt each microlens to face a point (Dx,Dy,Dz), where Dx,Dy,Dz are LDE parameters
 More accurately, the tilt is towards a ray the emanates from (Dx,Dy,Dz) and refracts from the previous surface.
 Jan 3rd, 2021

 Modified by Omer Korech to enable shifting the apexes according to the chief ray of previous RMLA (Regular Multi Lens Array aka EvenArray)
 Jan 15th, 2021

 Modified by Omer Korech to enable variable curvature and conic constant
 March 2nd, 2021

This surface breaks up the beam into numerous separate
beams, and so most OpticStudio features will fail to work with this surface.

However, the spot diagrams, image analysis, etc, all work okay.

*/

/* Gets the ray x and y as inputs, and should provide the nearest apex coordinates cx, cy as output s*/
int GetCellCenter(double x, double y, double *Rcx, double *Rcy, double *DTcx, double *DTcy, FIXED_DATA5* FD);

/* Sag of untilted micro lens as a function of its local x and y (x and y relative to the micro lens optical axis) */
double SagOfUntiltedAsphere(double x, double y, FIXED_DATA5* FD, double LensDistance);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);

/*Get surface normal coordinates (l,m,n) of untilted micro lens as a funtion of its local x and y */
int CalculateSurfaceNormalOfUntiltedMicroLens(double x, double y, double* l, double* m, double* n, FIXED_DATA5* FD, double LensDistance);

/* Propogates a ray. pV is a pointer to V. V may stand for = x, y, z, l, m, n. V is in the microlens coordinate frame of reference */
int PropogateThroughUntiltedAsphere(double* px, double* py, double* pz, double* pl, double* pm, double* pn, double* path, FIXED_DATA5* FD, double LensDistance);

/*GtoL converts global coordiates, where the lens is tilted, to lens coordinates, where the lens is rotationally symetric */
int GtoL(double * x, double * y, double * z, double Rcx, double Rcy, FIXED_DATA5* FD);

/*LtoG converts lens coordinates, where the lens is rotationally symetric, to global coordiates, where the lens is tilted */
int LtoG(double* x, double* y, double* z, double Rcx, double Rcy, FIXED_DATA5* FD);

int RoundToInt(double num);

int  __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD)
{
	int i, error;
	double power, ApexDistance;
	double Rcx, Rcy, DTcx, DTcy, x, y, z, l, m, n; // R stands for RMLA (Regular MicroLens Array), DT stands for Dectentered & Tilted MLA
		
	switch (FD->type)
	{
	case 0:
		/* ZEMAX is requesting general information about the surface */
		switch (FD->numb)
		{
		case 0:
			/* ZEMAX wants to know the name of the surface */
			/* do not exceed 12 characters */
			strcpy(UD->string, "DTMLA");
			break;
		case 1:
			/* ZEMAX wants to know if this surface is rotationally symmetric */
			/* it is not, so return a null string */
			UD->string[0] = '\0';
			break;
		case 2:
			/* ZEMAX wants to know if this surface is a gradient index media */
			/* it is not, so return a null string */
			UD->string[0] = '\0';
			break;
		}
		break;
	case 1:
		/* ZEMAX is requesting the names of the parameter columns */
		/* the value FD->numb will indicate which value ZEMAX wants. */
		/* they are all "Unused" for this surface type */
		/* returning a null string indicates that the parameter is unused. */
		switch (FD->numb)
		{
		case 1:
			strcpy(UD->string, "2nd Order Term");
			break;
		case 2:
			strcpy(UD->string, "4th Order Term");
			break;
		case 3:
			strcpy(UD->string, "6th Order Term");
			break;
		case 4:
			strcpy(UD->string, "8th Order Term");
			break;
		case 5:
			strcpy(UD->string, "10th Order Term");
			break;
		case 6:
			strcpy(UD->string, "RMLA Thickness");
			break;
		case 7:
			strcpy(UD->string, "Ar");
			break;
		case 8:
			strcpy(UD->string, "Br");
			break;
		case 9:
			strcpy(UD->string, "RMLA Wx");
			break;
		case 10:
			strcpy(UD->string, "RMLA Wy");
			break;
		case 11:
			strcpy(UD->string, "Dx");
			break;
		case 12:
			strcpy(UD->string, "Dy");
			break;
		case 13:
			strcpy(UD->string, "Dz");
			break;
		case 14:
			strcpy(UD->string, "Ac");
			break;
		case 15:
			strcpy(UD->string, "Bc");
			break;
		default:
			UD->string[0] = '\0';
			break;
		}
		break;
	case 2:
		/*
		This case is now deprecated. It used to be necessary when inputs for a surface
		needed to be separately defined as "parameter" data and "extra" data. This is no
		longer the case, as OpticStudio now supports the ability to define all inputs as
		"parameter" data, as long as the FIXED_DATA5 structure is being used by the DLL.
		*/
		break;
	case 3:
		/* ZEMAX wants to know the sag of the surface */
		/* if there is an alternate sag, return it as well */
		/* otherwise, set the alternate sag identical to the sag */
		/* The sag is sag1, alternate is sag2. */

		UD->sag1 = 0.0;
		UD->sag2 = 0.0;

		/* figure out the center coordinates of which "cell" we are in */
				
		x = UD->x;
		y = UD->y;
		z = UD->z;

		l = 0;
		m = 0;
		n = 1;

		/* make sure nx and ny are both odd, otherwise the chief ray is a problem...*/
		
		if (GetCellCenter(x, y, &Rcx, &Rcy, &DTcx, &DTcy, FD) ){
			return -1;
		}
		
		/* offset the coordinates */
		x -= DTcx;
		y -= DTcy;	

		/*For V8 (Variable radius and conic) Calculate lens distance from rotation axis*/
		ApexDistance = sqrt((DTcx - FD->param[11]) * (DTcx - FD->param[11]) + (DTcy - FD->param[12]) * (DTcy - FD->param[12]));

		//	Dx = FD->param[11];
		//	Dy = FD->param[12];

		error = GtoL(&x, &y, &z, Rcx, Rcy, FD); //Rotates x, y, z to a coordinate system where the lens is untilted 
		error = error || GtoL(&l, &m, &n, Rcx, Rcy, FD); //Rotates l, m, n to a coordinate system where the lens is untilted 
		PropogateThroughUntiltedAsphere(&x, &y, &z, &l, &m, &n, &power, FD, ApexDistance); // TIR error are irrelevant for case 3. power is a dummy variable
		error = error || LtoG(&x, &y, &z, Rcx, Rcy, FD); //Rotates x, y, z back to a coordinate system where the lens is tilted 		

		if (error) {
			return -1;
		}
		else {
			UD->sag1 = z;
			UD->sag2 = UD->sag1;
		}
		break;

	case 4:
		/* ZEMAX wants a paraxial ray trace to this surface */
		/* x, y, z, and the optical path are unaffected, at least for this surface type */
		/* for paraxial ray tracing, the return z coordinate should always be zero. */
		/* paraxial surfaces are always planes with the following normals */

		/* for a lens array, only consider the single lens on axis */
		UD->ln = 0.0;
		UD->mn = 0.0;
		UD->nn = -1.0;
		power = (FD->n2 - FD->n1)*(FD->cv + 2.0*FD->param[1]);

		if ((UD->n) != 0.0)
		{
			(UD->l) = (UD->l) / (UD->n);
			(UD->m) = (UD->m) / (UD->n);

			(UD->l) = (FD->n1*(UD->l) - (UD->x)*power) / (FD->n2);
			(UD->m) = (FD->n1*(UD->m) - (UD->y)*power) / (FD->n2);

			/* normalize */
			(UD->n) = sqrt(1 / (1 + (UD->l)*(UD->l) + (UD->m)*(UD->m)));
			/* de-paraxialize */
			(UD->l) = (UD->l)*(UD->n);
			(UD->m) = (UD->m)*(UD->n);
		}
		break;
	case 5:
		/* ZEMAX wants a real ray trace to this surface */

		x = UD->x;
		y = UD->y;
		z = UD->z;
		
		if (GetCellCenter(x, y, &Rcx, &Rcy, &DTcx, &DTcy, FD)) {
			return -1;
		}

		/* offset the coordinates */
		x -= DTcx;
		y -= DTcy;

		/*For V8 (Variable radius and conic) Calculate lens distance from rotation axis*/
		ApexDistance = sqrt((DTcx - FD->param[11]) * (DTcx - FD->param[11]) + (DTcy - FD->param[12]) * (DTcy - FD->param[12]));
		
		error = GtoL(&x, &y, &z, Rcx, Rcy, FD); //Rotates x, y, z to a coordinate system where the lens is untilted 
		error = error || GtoL(&UD->l, &UD->m, &UD->n, Rcx, Rcy, FD); //Rotates l, m, n to a coordinate system where the lens is untilted 
		if (error)
		{
			return (-1);
		}
		error =PropogateThroughUntiltedAsphere(&x, &y, &z, &UD->l, &UD->m, &UD->n, &UD->path, FD, ApexDistance);
		if (error)
		{
			return (error);
		}

		error = LtoG(&x, &y, &z, Rcx, Rcy, FD); //Rotates x, y, z back to a coordinate system where the lens is tilted 
		error = error || LtoG(&UD->l, &UD->m, &UD->n, Rcx, Rcy, FD); //Rotates l, m, n back to a coordinate system where the lens is tilted

		if (error)
		{
			return (-1);
		}
		
		// CalculateSurfaceNormalOfUntiltedMicroLens: // Might not be nessecery
	
		/* restore coordinates */
		UD->x = x + DTcx;
		UD->y = y + DTcy;
		UD->z = z;

		break;

	case 6:
		/* ZEMAX wants the index, dn/dx, dn/dy, and dn/dz at the given x, y, z. */

		/* This is only required for gradient index surfaces, so return dummy values */
		UD->index = FD->n2;
		UD->dndx = 0.0;
		UD->dndy = 0.0;
		UD->dndz = 0.0;
		break;
	case 7:
		/* ZEMAX wants the "safe" data. */
		/* this is used by ZEMAX to set the initial values for all parameters and extra data */
		/* when the user first changes to this surface type. */
		/* this is the only time the DLL should modify the data in the FIXED_DATA FD structure */
		for (i = 1; i < 16; i++) FD->param[i] = 0.0;
		
		FD->param[9] = 0.01; //RMLA Wx
		FD->param[10] = 0.01; //RMLA Wy
		FD->param[11] = 0; //Dx
		FD->param[12] = 0; //Dy
		FD->param[13] = 350.0; //Dz			
		break;
	case 8:
		/* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
		break;
	case 9:
		/* ZEMAX is calling the DLL for the last time, do any memory release here. */
		break;
	}
	return 0;
}

int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn)
{
	double nr, cosi, cosi2, rad, cosr, gamma;
	if (thisn != nextn)
	{
		nr = thisn / nextn;
		cosi = fabs((*l) * ln + (*m) * mn + (*n) * nn);
		cosi2 = cosi * cosi;
		if (cosi2 > 1) cosi2 = 1;
		rad = 1 - ((1 - cosi2) * (nr * nr));
		if (rad < 0) {
			return(-1);
		}
		cosr = sqrt(rad);
		gamma = nr * cosi - cosr;
		(*l) = (nr * (*l)) + (gamma * ln);
		(*m) = (nr * (*m)) + (gamma * mn);
		(*n) = (nr * (*n)) + (gamma * nn);
	}
	return 0;
}

int GetCellCenter(double x, double y, double *Rcx, double *Rcy, double *DTcx, double *DTcy, FIXED_DATA5* FD)
{
	double RXr, RYr, RZr, Magnitue;
	double RT, RWx, RWy, Dx, Dy, Dz;
	double D, IntersectOfStraightLine, InitialGuess;
	double Umin = 0;
	double l, m, n;
	double LowerBound, UpperBound, middle;
	int iter = 0;

	RT = FD->param[6];
	RWx = FD->param[9];
	RWy = FD->param[10];
	if (RWx <= 0.0 || RWy <= 0.0) {
		return(-1);
	}

	Dx = FD->param[11];
	Dy = FD->param[12];
	Dz = FD->param[13];
	if (Dz < RT) {
		return (-1);
	}

	/* Finds the interseption point at the RMLA plane
	of a ray that emenates in the pupil center,
	refracts by a plane in the RMLA location (disregarding the RMLA lens structure)
	and ends at FD->x, FD-y, FD-z
	The location is reported as the distance ( Umin ),
	along the 1D projection on the RMLA plane
	of the segment of the line from the start point (pupil center) to the end point (UD->x,y,z)*/

	D = sqrt((x - Dx)*(x - Dx) + (y - Dy)*(y - Dy)); // The length of the projection on the RMLA plane of the segment of the line from the start point(pupil center) to the end point(UD->x, y, z)

	if (D > Center_location_resolution / 100.0) { //This is usually true, on axis rays is a common expample for an exception.

		IntersectOfStraightLine = Dz * D / (RT + Dz); // Intersection of the ray if there was no refraction

		UpperBound = D;
		LowerBound = IntersectOfStraightLine;

		while (UpperBound - LowerBound > Center_location_resolution) { // 1.0e-6 --> 1 nm accuracy
			Umin = (UpperBound + LowerBound) / 2.0;

			if (SnellSign(Umin) > 0) {
				UpperBound = Umin;
			}
			else
			{
				LowerBound = Umin;
			}

			if (iter++ > ITMAX) break;
		}

		//Refine the evaluation:

		if (UpperBound < LowerBound) {
			middle = UpperBound;
			UpperBound = LowerBound;
			LowerBound = middle;
		}

		while (fabs(SnellSign(Umin) > Center_location_resolution / 10.0)) {

			middle = (UpperBound + LowerBound) / 2.0;

			Umin = middle - 0.5 * ((middle - LowerBound)*(middle - LowerBound) * (SnellSign(middle) - SnellSign(UpperBound)) - (middle - UpperBound)*(middle - UpperBound)*(SnellSign(middle) - SnellSign(LowerBound))) / (middle - LowerBound)* (SnellSign(middle) - SnellSign(UpperBound)) - (middle - UpperBound)*(SnellSign(middle) - SnellSign(LowerBound));

			if (Umin > middle) {
				UpperBound = Umin;
			}
			else
			{
				LowerBound = Umin;
			}

			if (iter++ > 1000) break;
		}


		// Sainity check	
		if (Umin > D || Umin < 0) {
			return (-1);
		}

		// Check that the result satisfy snell's law
		if (fabs(Umin / sqrt(Dz*Dz + Umin * Umin) - (FD->n1)*(D - Umin) / sqrt(RT*RT + (D - Umin)*(D - Umin))) > 0.0001) {
			return (-1);
		}


		/*Back to 2D
		find the vector coordinates of the intersection point, in a global coordinate frame
		(apriori the equations seem more elegant in the pupil frame of reference.
		However, the discretization in the lens location is conducted of the global frame of reference)
		Begin at Dx,Dy and propoaget Umin/D of the distance towards x,y*/
		x = Dx + Umin / D * (x - Dx);
		y = Dy + Umin / D * (y - Dy);


	}

	//Find the center of nearest RMLA micro lens
	*Rcx = RWx * RoundToInt(x / RWx);
	*Rcy = RWy * RoundToInt(y / RWy);
	
	// Find the center of nearest RMLA micro lens, relative to rotation point 
	RXr = *Rcx - Dx;
	RYr = *Rcy - Dy;
	RZr = Dz;		
	
	//Calculate the directional cosines of a ray from the rotation point to the RMLA apex
	Magnitue = sqrt(RXr * RXr + RYr * RYr + RZr * RZr);
	l = RXr / Magnitue;
	m = RYr / Magnitue;
	n = RZr / Magnitue;
	Refract(1, FD->n1, &l, &m, &n, 0, 0, -1);
	
	// DTMLA uLens coordinate
	*DTcx = *Rcx + l * RT / n;	
	*DTcy = *Rcy + m * RT / n;	

	return 0;
}

double SagOfUntiltedAsphere(double x, double y, FIXED_DATA5* FD, double LensDistance)
{
	double r2, alpha, Sag, a[6];
	int i;

	double Conic, Curvature;
	Curvature = FD->param[7] * LensDistance * LensDistance + FD->param[8] * LensDistance + FD->cv;
	Conic = FD->param[14] * LensDistance * LensDistance + FD->param[15] * LensDistance + FD->k;

	r2 = x * x + y * y;
	alpha = 1 - (1 + Conic) * Curvature * Curvature * r2;
	if (alpha < 0)
	{
		/* lens is not big enough to fill the cell, radius is too short */
		/* assume it is plane between lenses */
		alpha = 0.0;
		//r2 = 0.0; // Flatten the pillars between the lens
	}
	Sag = (Curvature * r2) / (1 + sqrt(alpha));

	/* now the aspheric terms */

	/* aspheric terms*/
	for (i = 1; i <= 5; i++)
	{
		a[i] = FD->param[i]; // even terms
	}

	Sag += a[1] * r2 + a[2] * r2 * r2 + a[3] * r2 * r2 * r2 + a[4] * r2 * r2 * r2 * r2 + a[5] * r2 * r2 * r2 * r2 * r2;

	return Sag;
}

int CalculateSurfaceNormalOfUntiltedMicroLens(double x, double y, double* l, double* m, double* n, FIXED_DATA5* FD, double LensDistance)
{
	
	double r2, alpha, mm, mx, my, a[6];
	int i;

	double Conic, Curvature;
	Curvature = FD->param[7] * LensDistance * LensDistance + FD->param[8] * LensDistance + FD->cv;
	Conic = FD->param[14] * LensDistance * LensDistance + FD->param[15] * LensDistance + FD->k;

	r2 = x * x + y * y;
	if (r2 < Center_location_resolution)
	{
		(*l) = 0;
		(*m) = 0;
		(*n) = -1;
	}
	else
	{
		alpha = 1.0 - (1.0 + Conic) * Curvature * Curvature * r2;
		if (alpha < 0) {	
			//return(FD->surf); // Ray misses the surface
			// The surfaces is approxiamated as flat
			(*l) = 0;
			(*m) = 0;
			(*n) = -1;
			return 0;			
		}
		alpha = sqrt(alpha);

		// mm * r = dz/dr
		mm = (Curvature / (1.0 + alpha)) * (2.0 + (Curvature * Curvature * r2 * (1.0 + Conic)) / (alpha * (1.0 + alpha)));

		/* now the aspheric terms */			
		for (i = 1; i <= 5; i++)
		{
			a[i] = FD->param[i]; // even terms
		}
		mm += 2.0 * a[1] + 4.0 * a[2] * r2 + 6.0 * a[3] * r2 * r2  + 8.0 * a[4] * r2 * r2 * r2 + 10.0 * a[5] * r2 * r2 * r2 * r2;

		mx = x * mm;  // = (x/r) * ( dz/dr )
		my = y * mm; // = (y/r) * ( dz/dr )
		
		(*n) = -sqrt(1 / (1 + (mx * mx) + (my * my)));
		(*l) = -mx * (*n);
		(*m) = -my * (*n);
	}
	return 0;
}

int PropogateThroughUntiltedAsphere(double * px, double * py, double * pz, double * pl, double * pm, double * pn, double * path, FIXED_DATA5* FD, double LensDistance)
{
	int loop;
	double t, dz, sag;
	double ln, mn,nn;

	/* make sure we do at least 1 loop */
	t = 100.0;
	loop = 0;


	while (fabs(t) > 1e-10)
	{
		sag = SagOfUntiltedAsphere(*px, *py, FD, LensDistance);

		/* okay, now with sag in hand, how far are we away in z? */
		dz = sag - *pz;

		/* now compute how far along the beam trajectory this is */
		/* note this will crash if n == 0!! */
		t = dz / (*pn);

		/* propagate the additional "t" distance */
		*px = *px + (*pl) * t;
		*py = *py + (*pm) * t;
		*pz = *pz + (*pn) * t;

		/* add in the optical path */
		*path= *path + t;

		/* prevent infinte loop if no convergence */
		loop++;
		if (loop > 5000) {  // Extending this number 5000 --> more than 5000, may increase prescision
			if( fabs(t) > 1e-6)//More than 1 nm between beam and surface
			{ 
				return(-1);
			}
			else {
				t = 0.0; // Round the less than 1 nm path difference to 0.
			}
			
		}
	}

	/* now do the normals */	
	if (CalculateSurfaceNormalOfUntiltedMicroLens(*px, *py, &ln, &mn, &nn, FD, LensDistance)) {
		return(-1);
	}

	if (Refract(FD->n1, FD->n2, pl, pm, pn, ln, mn, nn))
	{
		return(-FD->surf);							
	}
	return 0;
}

int GtoL(double* x, double* y, double* z, double Rcx, double Rcy, FIXED_DATA5* FD)
{
	double Xg, Yg, Zg; 
	double Dx, Dy, Dz;	
	double Xr, Yr, Zr, r, R;
	double SinTheta, CosTheta, SinPhi, SinTilt, CosTilt;

	Dx = FD->param[11];
	Dy = FD->param[12];
	Dz = FD->param[13];

	Xr = Rcx - Dx;
	Yr = Rcy - Dy;
	Zr = Dz;
	r = sqrt(Xr * Xr + Yr * Yr + Zr * Zr);
	R = sqrt(Xr * Xr + Yr * Yr);

	if (fabs(Zr) < 1e-6) {
		return(-1);
	}
	
	if (R < Center_location_resolution) {
		return(0);
	}
	
	SinTheta = Yr / R;
	CosTheta = Xr / R;	
	SinPhi = R / r;
	SinTilt = SinPhi / (FD->n1);
	CosTilt = sqrt(1 - SinTilt * SinTilt);

	Xg = *x;
	Yg = *y;
	Zg = *z;

	*x = CosTilt * (Xr * Xg + Yr * Yg) / R - Zg * SinTilt;
	*y = (Xr * Yg - Yr * Xg) / R;
	*z = SinTilt * (Xr * Xg + Yr * Yg) / R + Zg * CosTilt;

	return 0;
}

int LtoG(double* x, double* y, double* z, double Rcx, double Rcy, FIXED_DATA5 * FD)
{
	double Xl, Yl, Zl; 
	double Dx, Dy, Dz;	
	double Xr, Yr, Zr, r, R;
	double SinTheta, CosTheta, SinPhi, SinTilt, CosTilt;

	Dx = FD->param[11];
	Dy = FD->param[12];
	Dz = FD->param[13];

	Xr = Rcx - Dx;
	Yr = Rcy - Dy;
	Zr = Dz;
	r = sqrt(Xr * Xr + Yr * Yr + Zr * Zr);
	R = sqrt(Xr * Xr + Yr * Yr);

	if (fabs(Zr) < 1e-6) {
		return(-1);
	}

	if (R < Center_location_resolution) {
		return(0);
	}

	SinTheta = Yr / R;
	CosTheta = Xr / R;
	SinPhi = R / r;
	SinTilt = SinPhi / (FD->n1);
	CosTilt = sqrt(1 - SinTilt * SinTilt);

	Xl = *x;
	Yl = *y;
	Zl = *z;

	*x = ( Xr * (Xl * CosTilt + Zl * SinTilt) - Yr * Yl ) / R;
	*y = ( Yr * (Xl * CosTilt + Zl * SinTilt) + Xr * Yl ) / R;
	*z = Zl * CosTilt - Xl * SinTilt;

	return 0;
}

int RoundToInt(double num)
{
	return (num >= 0) ? (int)(num + 0.5) : (int)(num - 0.5);
}









