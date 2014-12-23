#include "../gsop.cuh"
/*
 * indentation_kernel.cu
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */

struct IndentationConstant {
	float3 tipCoord;
	float3 chipCoord;
	float3 chipCoord0;
	float3 direction;
	float V;
	float dx;
	int moveSurface;
	int fixTransversal;
	float tipRadius;
	float cantileverKs;
	float tipa6;
	float tipel;
	float surfa6;
	float surfel;
	float tipAprime;
	float tipBprime;
	float surfAprime;
	float surfBprime;
	float tipZeta;
	float3 micaN;
	float3 micaR0;
	float3 micaR;
	float4* h_tipForces;
	float4* d_tipForces;
	float4 tipForce;
	FILE* out;
	long int retractionStep;
	int showTipSurface;
	int surfaceSize;
	int surfaceBeadsCount;
	float surfaceStep;
	float3* surfacePointsR0;
	float4* h_surfacePointsCoord;
	float4* d_surfacePointsCoord;

	float pairsCutoff2;
	int* h_micaListCounts;
	int* d_micaListCounts;
	int* h_micaList;
	int* d_micaList;

	int outputFreq;
	float3 cantileverVector;
	float4 fav;
	float3 tipCoordAv;
	float3 chipCoordAv;
	float kDeltaXAv;
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
