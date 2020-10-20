texture<float2, 2, cudaReadModeElementType> gDAdapt;
texture<float, 2, cudaReadModeElementType> gAdapt;

texture<float, 2, cudaReadModeElementType> tex_x;
texture<float, 2, cudaReadModeElementType> tex_y;

__global__ void SetZero(float *oData, int iSize)
{
	int x = 256*blockIdx.x+threadIdx.x;
	if (x<iSize)
		oData[x]=0;
}

__device__ void DDACore(float iX,
						float iY,
						float *data,
						float *resolution,
						int cmpIndex,
						size_t size,
						float DIS,
						int RESOLUTION)
{
	float dx=iX-data[cmpIndex*2];
	float dy=iY-data[cmpIndex*2+1];
	dx = dx-floorf(dx+0.5f);
	dy = dy-floorf(dy+0.5f);

	if((fabs(dx)<DIS)&&(fabs(dy)<DIS))
	{
		int index_x=floor(dx*RESOLUTION/2.0f/DIS+RESOLUTION/2);
		int index_y=floor(dy*RESOLUTION/2.0f/DIS+RESOLUTION/2);
		float addVal = (float)1/(float)size;
		atomicAdd(resolution+index_y*RESOLUTION+index_x, addVal);

	}
}

__global__ void DDA(float *data, float *resolution,size_t size,float DIS,int RESOLUTION)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<size)
	{
		float x = data[valindex*2];
		float y = data[valindex*2+1];
		for(int i=0;i<size;i++)
		if(valindex!=i)
		{
			DDACore(x, y, data, resolution, i, size, DIS, RESOLUTION);
		}
	}
//	__syncthreads();

}

__global__ void DDAAdapt(float *data, float *resolution,size_t size,float DIS,int RESOLUTION, float iAspect)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<size)
	{
		float2 pst0 = make_float2(data[valindex*2], data[valindex*2+1]);
		float r0 = tex2D(gAdapt, pst0.x, pst0.y);
		for(int i=0;i<size;i++)
		if(valindex!=i)
		{
			float2 pst1 = make_float2(data[i*2], data[i*2+1]);
			float r1 = tex2D(gAdapt, pst1.x, pst1.y);
			float dx=pst0.x-pst1.x;
			float dy=pst0.y-pst1.y;
#ifndef _NON_WRAP
			dx = dx-floorf(dx+0.5f);
			dy = dy-floorf(dy+0.5f);
#endif

			dy *= iAspect;

			float A = 2.0f/(r0+r1);
			dx *= A;
			dy *= A;

			if((dx*dx<DIS*DIS)&&(dy*dy<DIS*DIS))
			{
				int index_x=dx*RESOLUTION/2.0f/DIS+RESOLUTION/2.0f;
				int index_y=dy*RESOLUTION/2.0f/DIS+RESOLUTION/2.0f;
				float addVal = (float)1/(float)size;
				atomicAdd(resolution+index_y*RESOLUTION+index_x, addVal);
			}
		}
	}
//	__syncthreads();

}

__device__ void CalcMoveLenCore(float *iData,
								float iX,
								float iY,
								int iCmpIndex,
								float iDDARange,
								float& oAddX,
								float& oAddY
								)
{
	float dx_temp=iX-iData[iCmpIndex*2];
	float dy_temp=iY-iData[iCmpIndex*2+1];
	dx_temp = dx_temp-floorf(dx_temp+0.5f);
	dy_temp = dy_temp-floorf(dy_temp+0.5f);

	if((dy_temp*dy_temp+dx_temp*dx_temp)<iDDARange*iDDARange*0.81)
	{
		float u=dx_temp/2/iDDARange+0.5f;
		float v=dy_temp/2/iDDARange+0.5f;

		float mx = tex2D(tex_x, u, v);
		float my = tex2D(tex_y, u, v);

		oAddX += mx;
		oAddY += my;
	}
}

__device__ void CalcMoveLen(float *data,
							int valindex,
							size_t size,
							float DIS,
							int RESOLUTION,
							float& oMoveX,
							float& oMoveY)
{
	float addx = 0;
	float addy = 0;
	float x = data[valindex*2];
	float y = data[valindex*2+1];
	for(int i=0;i<size;i++)
	{   
		if(valindex!=i)
		{
			CalcMoveLenCore(data, x, y, i, DIS, addx, addy);
		}
	}
	oMoveX = addx;
	oMoveY = addy;
}

__global__ void CalcForce(float *data,
						  float2 *oForce,
						  float *oPower,
						  size_t size,
						  float DIS,
						  int RESOLUTION)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;

	if(valindex<size)
	{
		float addx=0;
		float addy=0;
		CalcMoveLen(data, valindex, size, DIS, RESOLUTION, addx, addy);
		oForce[valindex] = make_float2(addx, addy);
		oPower[valindex] = addx*addx+addy*addy;
	}

}

__global__ void MoveForce(float *data,
						  float *newData,
						  int size,
						  float2 *force,
						  float scale,
						  float disturbScale,
						  float *disturb)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;

	if(valindex<size)
	{
#ifdef _GLOBAL_DISTURB
		newData[valindex*2] = data[valindex*2]+force[valindex].x*scale+(disturb[valindex*2]-0.5f)*disturbScale;
		newData[valindex*2+1] = data[valindex*2+1]+force[valindex].y*scale+(disturb[valindex*2+1]-0.5f)*disturbScale;
#else
		newData[valindex*2] = data[valindex*2]+force[valindex].x*scale;
		newData[valindex*2+1] = data[valindex*2+1]+force[valindex].y*scale;
#endif
	}
}

__global__ void MoveAdapt(float *data,float *data_new, float* randNum,size_t size,float DIS,int RESOLUTION,float CONSTANT, int iRows, int iCols, float iAspect)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;

	if(valindex<size)
	{
		float mx = 0;
		float my = 0;
		float2 pst0 = make_float2(data[valindex*2], data[valindex*2+1]);
		float r0 = tex2D(gAdapt, pst0.x, pst0.y);
#ifdef _DER_ADJ
		float2 dR0 = tex2D(gDAdapt, pst0.x+0.5f/iCols, pst0.y+0.5f/iRows);
#else // _DER_ADJ
		float2 dR0 = tex2D(gDAdapt, pst0.x, pst0.y);
#endif // _DER_ADJ
		for(int i=0;i<size;i++)
		{ 
			if(valindex!=i)
			{
				float2 pst1 = make_float2(data[i*2], data[i*2+1]);
				float r1 = tex2D(gAdapt, pst1.x, pst1.y);
				float2 dpst = make_float2(pst1.x-pst0.x, pst1.y-pst0.y);
#ifndef _NON_WRAP
				dpst.x = dpst.x-floorf(dpst.x+0.5f);
				dpst.y = dpst.y-floorf(dpst.y+0.5f);
#endif
				dpst.y *= iAspect;

				float A = 2.0f/(r0+r1);
				float2 scaleDpst = make_float2(dpst.x*A, dpst.y*A);
				float len = sqrt(scaleDpst.x*scaleDpst.x+scaleDpst.y*scaleDpst.y);
				if (len<DIS*0.9f)
				{
					float tmp0 = -2.0f/(r0+r1)/(r0+r1);
					float2 dA = make_float2(tmp0*dR0.x, tmp0*dR0.y);
					float daxx = dA.x*dpst.x-A;
					float dayx = dA.x*dpst.y;
					float dayy = dA.y*dpst.y-A;
					float daxy = dA.y*dpst.x;
					float fx = tex2D(tex_x, scaleDpst.x/DIS/2+0.5f, scaleDpst.y/DIS/2+0.5f);
					float fy = tex2D(tex_y, scaleDpst.x/DIS/2+0.5f, scaleDpst.y/DIS/2+0.5f);

					mx += (fx*daxx+fy*dayx);
					my += (fy*dayy+fx*daxy);
				}

			}
		}
		float moveX;
		float moveY;
#ifdef _GLOBAL_DISTURB
		float randx = randNum[valindex*2]-0.5f;
		float randy = randNum[valindex*2+1]-0.5f;

		moveX = CONSTANT*(mx+randx*_DISTURB_FACTOR/sqrt((float)size))*r0;
		moveY = CONSTANT*(my+randy*_DISTURB_FACTOR/sqrt((float)size))*r0;
#else // _GLOBAL_DISTURB
		moveX = CONSTANT*mx*r0;
		moveY = CONSTANT*my*r0;
#endif // _GLOBAL_DISTURB
		moveY /= iAspect;

		data_new[valindex*2] += moveX;
		data_new[valindex*2+1] += moveY;
	}
}

__global__ void Change(float *data,float *data_new,size_t size)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<size)
	{
		data_new[valindex*2]-=floor(data_new[valindex*2]);
		data_new[valindex*2+1]-=floor(data_new[valindex*2+1]);
		data[valindex*2]=data_new[valindex*2];
		data[valindex*2+1]=data_new[valindex*2+1];

	}
}

__global__ void Gaussian(float* res, float* target,float *cov, float *oDiff,size_t size,int RESOLUTION,int Kernel_size,float Deviation, float* iGaussCoe)
{
    unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<size)
	{
		int kernelAlign = Kernel_size*2+1;
		cov[valindex]=0;
		int index_x=valindex%RESOLUTION;
		int index_y=valindex/RESOLUTION;
		float weight=0;
		float cmp = RESOLUTION*0.9f/2;
		cmp = cmp*cmp;
		float dx=index_x-RESOLUTION/2.0f+0.5f;
		float dy=index_y-RESOLUTION/2.0f+0.5f;
		if (dx*dx+dy*dy<cmp)
			oDiff[valindex]=res[valindex]-target[valindex];
		else
			oDiff[valindex]=0;
		for(int i=index_x-Kernel_size;i<=index_x+Kernel_size;i++)
		for(int j=index_y-Kernel_size;j<=index_y+Kernel_size;j++)
		{
			int index_new=j*RESOLUTION+i;
			float value=0;
			float x = i-RESOLUTION/2.0f+0.5f;
			float y = j-RESOLUTION/2.0f+0.5f;

			int gaussIndex = (i-index_x+Kernel_size)+(j-index_y+Kernel_size)*kernelAlign;
			float gaussVal = iGaussCoe[gaussIndex];
			if(x*x+y*y<cmp)
			   value=(res[index_new]-target[index_new])*gaussVal;
			cov[valindex]+=value;
			weight+=gaussVal;
		}
		cov[valindex]=cov[valindex]/weight;
	}
}


__global__ void  Gaussian2( float *cov1, float *cov_x, float *cov_y,size_t size,int RESOLUTION,int Kernel_size,float Deviation, float* iGaussDerCoe)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<size)
	{
		int kernelAlign = Kernel_size*2+1;
		cov_x[valindex]=0;
		cov_y[valindex]=0;
		int index_x=valindex%RESOLUTION;
		int index_y=valindex/RESOLUTION;
		for(int i=index_x-Kernel_size;i<=index_x+Kernel_size;i++)
			for(int j=index_y-Kernel_size;j<=index_y+Kernel_size;j++)

			{
				if(i>-1&&j>-1&&i<RESOLUTION&&j<RESOLUTION)
				{
					int index_new=j*RESOLUTION+i;
					int gaussIndex = (i-index_x+Kernel_size)+(j-index_y+Kernel_size)*kernelAlign;
					cov_x[valindex]+=cov1[index_new]*iGaussDerCoe[gaussIndex*2];
					cov_y[valindex]+=cov1[index_new]*iGaussDerCoe[gaussIndex*2+1];
				}

			}
	}
//	__syncthreads();
}

template <int BlockSize>
__global__ void InitGaussCoe(float* oGaussCoe, float* oGaussDerCoe, int iSize, float iSigma)
{
	int i = blockIdx.x*BlockSize+threadIdx.x;
	int j = blockIdx.y*BlockSize+threadIdx.y;
	if (i<iSize && j<iSize)
	{
		float scale = 1.0f/(2*3.141592654f*iSigma*iSigma);

		float x = i-iSize/2.0f+0.5f;
		float y = j-iSize/2.0f+0.5f;

		float gaussVal = scale * exp(-(x*x+y*y)/2.0f/iSigma/iSigma);

		int index = j*iSize+i;
		oGaussCoe[index]=gaussVal;
		oGaussDerCoe[index*2]=-gaussVal*gaussVal*x/iSigma/iSigma;
		oGaussDerCoe[index*2+1]=-gaussVal*gaussVal*y/iSigma/iSigma;
	}
}

