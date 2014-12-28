#include "../gsop.cuh"
/*
 * indentation_kernel.cu
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */

struct IndentationConstant {

	/*
	 * Cantilever tip parameters
	 */
	float3 tipCoord;  // Coordinates of the cantilever tip
	float3 chipCoord; // Coordinates of the cantilever base
	float3 chipCoord0; // Initial coordinates of the cantilever base
	float3 direction; // Direction of the cantilever base movement
	float V; // Velocity of the cantilever base (in A/indentation update period)
	float dx; // Current displacement of the cantilever base
	int moveSurface; // If the surface should be moved rather than cantilever base
	int fixTransversal; // If the cantilever tip should be fixed to the vector of the base movement (all perpendicular movement is supressed)
	float cantileverKs; // Spring constant for the spring connecting cantilever base with the tip
	float tipZeta;

	/*
	 * The interaction potentials between cantilever tip and amino acid and between the surface and amino acid are computed as
	 * U=epsilon*(A*(sigma/r)^12 + B*(sigma/r)^6),
	 * where epsilon (in kcal/mol) is the repulsion strength, sigma (in A) is the repulsion radii,
	 * A and B are dimentionless potential coefficients. Usually:
	 * A = 0, B = 1 for pure repulsion (U~1/r^6)
	 * A = 1, B = -2 for full LJ (U=epsilon*((sigma/r)^12 - 2*(sigma/r)^6)
	 * For the interaction between the tip and amino acid, r=r_ti-R, with r_ti is the distance between the tip center and the amino-acid, R is the tip Radii
	 * For the interaction between the surface and amino acid, r is the distance between the surface and the amino-acid
	 * For the computational purposes, the A'=12.0*epsilon*sigma^12*A and B'=6.0*epsilon*sigma^6*B are used
	 */
	float tipRadius; // Radius of the cantilever tip (in A) for the interaction between the tip and aminos
	float tipAprime; // A' for the tip-amino interactions
	float tipBprime; // B' for the tip-amino interactions
	float surfAprime; // A' for the surface-amino interactions
	float surfBprime; // B' for the surface-amino interactions

	/*
	 * Mica parameters
	 */
	float3 micaN; // Normal vector
	float3 micaR0; // A sample point on the surface to define its position at the start of simulation
	float3 micaR; // A sample point on the surface to define its position at the current timestep
	float4* h_tipForces; // Forces, acting on each amino-acid due to the interaction with the cantilever tip (on CPU)
	float4* d_tipForces; // Forces, acting on each amino-acid due to the interaction with the cantilever tip (on GPU)
	int surfaceSize; // Number of beads to add for the surface representation (along each axis, the total number would be squared)
	int surfaceBeadsCount; // Total number (surfaceSize^2)
	float surfaceStep; // Distance between beads that represent surface
	float3* surfacePointsR0; // Positions of all the beads representing the surface at the start
	float4* h_surfacePointsCoord; // Current position of all the beads representing the surface (on CPU)
	float4* d_surfacePointsCoord; // Current position of all the beads representing the surface (on GPU)

	/*
	 * Pairlist for the interactions of the amino-acids with the surface beads
	 */
	float pairsCutoff2;
	int* h_micaListCounts;
	int* d_micaListCounts;
	int* h_micaList;
	int* d_micaList;

};

IndentationConstant hc_indentation;
__device__ __constant__ IndentationConstant c_indentation;

__global__ void indentation_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i];
		float4 f = c_gsop.d_forces[d_i];
		float4 f_in = c_indentation.d_tipForces[d_i];//make_float4(0.0, 0.0, 0.0, 0.0);
		float4 df;
		float4 dr;
		dr.x =  coord.x - c_indentation.tipCoord.x;
		dr.y =  coord.y - c_indentation.tipCoord.y;
		dr.z =  coord.z - c_indentation.tipCoord.z;
		dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		float r2 = dr.w - c_indentation.tipRadius;
		r2 = 1.0f/r2;
		r2 = r2*r2;
		float r6 = r2*r2*r2;
		df.w = c_indentation.tipAprime*r6 + c_indentation.tipBprime;
		df.w = df.w*r6/(dr.w*(dr.w - c_indentation.tipRadius));
		df.x = dr.x*df.w;
		df.y = dr.y*df.w;
		df.z = dr.z*df.w;
		/*dr.w = sqrt(dr.x*dr.x+dr.y*dr.y+dr.z*dr.z);
		df.w = dr.w - c_indentation.tipRadius;
		df.w = powf(df.w, 7);
		//if(dr.w < c_pairs.pairsCutoff2){
			df.w = c_indentation.a6/df.w;
			df.w = -6.0*c_indentation.el*df.w/dr.w;
			df.x = dr.x*df.w;
			df.y = dr.y*df.w;
			df.z = dr.z*df.w;
		//}*/
		f_in.x += df.x;
		f_in.y += df.y;
		f_in.z += df.z;
		c_indentation.d_tipForces[d_i] = f_in;

		dr.x = coord.x - c_indentation.micaR.x;
		dr.y = coord.y - c_indentation.micaR.y;
		dr.z = coord.z - c_indentation.micaR.z;
		dr.w = dr.x*c_indentation.micaN.x + dr.y*c_indentation.micaN.y + dr.z*c_indentation.micaN.z;

		dr.w = 1.0f/dr.w;
		dr.w = dr.w*dr.w;
		r6 = dr.w*dr.w*dr.w;
		df.w = c_indentation.surfAprime*r6 + c_indentation.surfBprime;
		df.w = df.w*r6*dr.w;
		f.x += c_indentation.micaN.x*df.w + df.x;
		f.y += c_indentation.micaN.y*df.w + df.y;
		f.z += c_indentation.micaN.z*df.w + df.z;

		/*dr.w = powf(dr.w, 8);
		dr.w = 6.0*c_indentation.el*c_indentation.a6/dr.w;
		f.x += c_indentation.micaN.x*dr.w + df.x;
		f.y += c_indentation.micaN.y*dr.w + df.y;
		f.z += c_indentation.micaN.z*dr.w + df.z;*/
		c_gsop.d_forces[d_i] = f;
	}
}

__global__ void indentationDiscreteSurf_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i];
		float4 f = c_gsop.d_forces[d_i];
		float4 dr;
		dr.x =  coord.x - c_indentation.tipCoord.x;
		dr.y =  coord.y - c_indentation.tipCoord.y;
		dr.z =  coord.z - c_indentation.tipCoord.z;
		dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		coord.w = dr.w - c_indentation.tipRadius;
		coord.w = 1.0f/coord.w;
		coord.w = coord.w*coord.w;
		float r6 = coord.w*coord.w*coord.w;
		coord.w = c_indentation.tipAprime*r6 + c_indentation.tipBprime;
		coord.w = coord.w*r6/(dr.w*(dr.w - c_indentation.tipRadius));
		dr.x = dr.x*coord.w;
		dr.y = dr.y*coord.w;
		dr.z = dr.z*coord.w;
		f.x += dr.x;
		f.y += dr.y;
		f.z += dr.z;

		coord = c_indentation.d_tipForces[d_i];
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;
		c_indentation.d_tipForces[d_i] = coord;
		coord = c_gsop.d_coord[d_i];

		int i;
		for(i = 0; i < c_indentation.d_micaListCounts[d_i]; i++){
			dr = c_indentation.d_surfacePointsCoord[c_indentation.d_micaList[i*c_gsop.width + d_i]];

			dr.x = coord.x - dr.x;
			dr.y = coord.y - dr.y;
			dr.z = coord.z - dr.z;
			dr.w = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

			dr.w = 1.0f/dr.w;
			r6 = dr.w*dr.w*dr.w;
			coord.w = c_indentation.surfAprime*r6 + c_indentation.surfBprime;
			coord.w = coord.w*r6*dr.w;
			f.x += dr.x*coord.w;
			f.y += dr.y*coord.w;
			f.z += dr.z*coord.w;
		}

		c_gsop.d_forces[d_i] = f;
	}
}

__global__ void generateMicaList_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int i;
		int count = 0;
		float4 coord = c_gsop.d_coord[d_i];
		for(i = 0; i < c_indentation.surfaceBeadsCount; i++){
			float4 r2 = c_indentation.d_surfacePointsCoord[i];
			r2.x = coord.x - r2.x;
			r2.y = coord.y - r2.y;
			r2.z = coord.z - r2.z;
			r2.w = r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
			if(r2.w < c_indentation.pairsCutoff2){
				c_indentation.d_micaList[count*c_gsop.width + d_i] = i;
				count ++;
			}
		}
		c_indentation.d_micaListCounts[d_i] = count;
	}
}
