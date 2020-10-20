int gCellReso = 1;
int* gCellBeg = NULL;
int* gCellEnd = NULL;
int* gSampleKeys = NULL;

__global__ void UpdateKeys(float2* iData, int iSize, int* oKeys, int iCellReso);

__global__ 
void FindCellStartEnd(int* iKeys,
					  int iKeyNum,
					  int* oCellStart,
					  int* oCellEnd)
{
	extern __shared__ int sharedHash[];    // blockSize + 1 elements
    int index = blockIdx.x*blockDim.x + threadIdx.x;
	
    int key;
    // handle case when no. of particles not multiple of block size
    if (index < iKeyNum) {
        key = iKeys[index];

        // Load hash data into shared memory so that we can look 
        // at neighboring particle's hash value without loading
        // two hash values per thread
	    sharedHash[threadIdx.x+1] = key;

	    if (index > 0 && threadIdx.x == 0)
	    {
		    // first thread in block must load neighbor particle hash
		    sharedHash[0] = iKeys[index-1];
	    }
	}

	__syncthreads();

	if (index < iKeyNum) {
		// If this particle has a different cell index to the previous
		// particle then it must be the first particle in the cell,
		// so store the index of this particle in the cell. 
		// As it isn't the first particle, it must also be the cell end of
		// the previous particle's cell

	    if (index == 0 || key != sharedHash[threadIdx.x])
	    { 
		    oCellStart[key] = index;
            if (index > 0)
                oCellEnd[sharedHash[threadIdx.x]] = index;
	    }

        if (index == iKeyNum - 1)
        {
            oCellEnd[key] = index + 1;
        }
	}
}

int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

void computeGridSize(int n, int blockSize, int &numBlocks, int &numThreads)
{
    numThreads = min(blockSize, n);
    numBlocks = iDivUp(n, numThreads);
}

void InitGrid(int iSampleSize)
{
	cudaMalloc((void**)&gCellBeg, gCellReso*gCellReso*sizeof(int));
	cudaMalloc((void**)&gCellEnd, gCellReso*gCellReso*sizeof(int));
	cudaMalloc((void**)&gSampleKeys, iSampleSize*sizeof(int));
}

__device__ __host__
int CalcCellIDInt(int iCellX, int iCellY, int iCellReso)
{
	return iCellX+iCellY*iCellReso;
}

__device__ __host__
int CalcCellID(float ix, float iy, int iCellReso)
{
	//ix = floor(ix*iCellReso);
	//iy = floor(iy*iCellReso);
	//int x = ix;
	//int y = iy;
	//return CalcCellIDInt(x, y, iCellReso);
	return CalcCellIDInt(int(floor(ix*iCellReso)), int(floor(iy*iCellReso)), iCellReso);
}

void BuildGrid(float2* ioData,
			   int* ioKeys,
			   int iDataNum,
			   int* oCellStart,
			   int* oCellEnd
			   )
{
    int numThreads, numBlocks;
    computeGridSize(iDataNum, 256, numBlocks, numThreads);
	
	UpdateKeys<<<numBlocks, numThreads>>>(ioData, iDataNum, ioKeys, gCellReso);

	int celNum = gCellReso*gCellReso;
	thrust::sort_by_key(
		thrust::device_ptr<int>(ioKeys),
		thrust::device_ptr<int>(ioKeys+iDataNum),
		thrust::device_ptr<float2>(ioData));
	cudaMemset(oCellStart, 0, celNum*sizeof(int));
	cudaMemset(oCellEnd, 0, celNum*sizeof(int));

    int smemSize = sizeof(int)*(numThreads+1);
	FindCellStartEnd<<<numBlocks, numThreads, smemSize>>>
		(ioKeys, iDataNum, oCellStart, oCellEnd);

	//thrust::host_vector<int> cellStart(thrust::device_ptr<int>(oCellStart),
	//	thrust::device_ptr<int>(oCellStart+gCellReso*gCellReso));
	//thrust::host_vector<int> cellEnd(thrust::device_ptr<int>(oCellEnd),
	//	thrust::device_ptr<int>(oCellEnd+gCellReso*gCellReso));
	//for (int checkI = 0; checkI < gCellReso*gCellReso; checkI++)
	//{
	//	thrust::host_vector<float2> data(thrust::device_ptr<float2>(ioData+cellStart[checkI]),
	//		thrust::device_ptr<float2>(ioData+cellEnd[checkI]-1));
	//	for (int i=0; i<data.size(); i++)
	//		if (CalcCellID(data[i].x, data[i].y, gCellReso)!=checkI)
	//			printf("!");
	//}
	////for (int i=0; i<cellStart.size(); i++)
	////	printf("(%d,%d)\t", cellStart[i], cellEnd[i]);
	////printf("\n");
	//printf("\n");
}

__device__ __host__
int CalcCellBegY(float iy, float iDDARange, int iCellReso)
{
	//return 0;
	return int(floor((iy-iDDARange)*iCellReso));
}

__device__ __host__
int CalcCellEndY(float iy, float iDDARange, int iCellReso)
{
	return int(ceil((iy+iDDARange)*iCellReso));
}

__device__
int CalcCellBegX(float ix, float iy, int iCellY, float iDDARange, int iCellReso)
{
	//return 0;
	float yup = iy-float(iCellY)/iCellReso;
	float ydown = iy-float(iCellY+1)/iCellReso;
	float ylen = min(fabs(yup), fabs(ydown));
	if (yup*ydown<=0)
		ylen = 0;
	float dx = sqrt(iDDARange*iDDARange-ylen*ylen);
	return int( floor((ix-dx)*iCellReso) );
}

__device__
int CalcCellEndX(float ix, float iy, int iCellY, float iDDARange, int iCellReso)
{
	//return iCellReso-1;
	float yup = iy-float(iCellY)/iCellReso;
	float ydown = iy-float(iCellY+1)/iCellReso;
	float ylen = min(fabs(yup), fabs(ydown));
	if (yup*ydown<=0)
		ylen = 0;
	float dx = sqrt(iDDARange*iDDARange-ylen*ylen);
	return int( ceil((ix+dx)*iCellReso) );
}

__device__
int wrap(int iVal, int iHighBound)
{
	int retVal = iVal - iHighBound*(iVal/iHighBound);
	retVal += iHighBound;
	return retVal%iHighBound;
}

__global__ void DDAGrid(float2* iData,
					  int iDataSize,
					  int iCellReso,
					  int* iCellIndexBeg,
					  int* iCellIndexEnd,
					  float* oHist,
					  int iHistReso,
					  float iDDARange)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<iDataSize)
	{
		float2 pst = iData[valindex];
		float dataX = pst.x;
		float dataY = pst.y;
		int cellYBeg = CalcCellBegY(dataY, iDDARange, iCellReso);
		int cellYEnd = CalcCellEndY(dataY, iDDARange, iCellReso);
		for (int scellY = cellYBeg; scellY <= cellYEnd; scellY++)
		{
			int cellY = wrap(scellY, iCellReso);
			int cellXBeg = CalcCellBegX(dataX, dataY, scellY, iDDARange, iCellReso);
			int cellXEnd = CalcCellEndX(dataX, dataY, scellY, iDDARange, iCellReso);
			for (int scellX = cellXBeg; scellX <= cellXEnd; scellX++)
			{
				int cellX = wrap(scellX, iCellReso);
				int cellID = CalcCellIDInt(cellX, cellY, iCellReso);
				int begIdx = iCellIndexBeg[cellID];
				int endIdx = iCellIndexEnd[cellID];
				for (int i=begIdx; i<endIdx; i++)
				if (i!=valindex)
				{
					DDACore(dataX, dataY, (float*)iData, oHist, i, iDataSize, iDDARange, iHistReso);
				}
			}
		}
	}
}

__global__ void ChangeAndCount(float *data,float *data_new,size_t size,int iCellReso, int* oCount)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<size)
	{
		float dataOldX = data[valindex*2];
		float dataOldY = data[valindex*2+1];
		float dataNewX = data_new[valindex*2];
		float dataNewY = data_new[valindex*2+1];
		dataNewX -= floor(dataNewX);
		dataNewY -= floor(dataNewY);

		oCount[valindex] =
			(CalcCellID(dataNewX, dataNewY, iCellReso)!=CalcCellID(dataOldX, dataOldY, iCellReso));

		data[valindex*2]=data_new[valindex*2]=dataNewX;
		data[valindex*2+1]=data_new[valindex*2+1]=dataNewY;

	}
}

__global__ void UpdateKeys(float2* iData, int iSize, int* oKeys, int iCellReso)
{
	int valindex = blockIdx.x * 256 + threadIdx.x;
	if (valindex<iSize)
	{
		float2 pst = iData[valindex];
		oKeys[valindex] = CalcCellID(pst.x,
							pst.y,
							iCellReso);
	}
}

__global__ void CalcForceGrid(float2 *iData,
							int iDataSize,
							float2* oForce,
							float* oPower,
							float iDDARange,
							int* iCellIndexBeg,
							int* iCellIndexEnd,
							int iCellReso
								)
{
	unsigned long valindex = blockIdx.x * 256 + threadIdx.x;
	if(valindex<iDataSize)
	{
		float addX = 0;
		float addY = 0;
		float2 pst = iData[valindex];
		float dataX = pst.x;
		float dataY = pst.y;
		int cellYBeg = CalcCellBegY(dataY, iDDARange, iCellReso);
		int cellYEnd = CalcCellEndY(dataY, iDDARange, iCellReso);
		for (int scellY = cellYBeg; scellY <= cellYEnd; scellY++)
		{
			int cellY = wrap(scellY, iCellReso);
			int cellXBeg = CalcCellBegX(dataX, dataY, scellY, iDDARange, iCellReso);
			int cellXEnd = CalcCellEndX(dataX, dataY, scellY, iDDARange, iCellReso);
			for (int scellX = cellXBeg; scellX <= cellXEnd; scellX++)
			{
				int cellX = wrap(scellX, iCellReso);
				int cellID = CalcCellIDInt(cellX, cellY, iCellReso);
				int begIdx = iCellIndexBeg[cellID];
				int endIdx = iCellIndexEnd[cellID];
				for (int i=begIdx; i<endIdx; i++)
				if (i!=valindex)
				{
					CalcMoveLenCore(
						(float*)iData,
						dataX,
						dataY,
						i,
						iDDARange,
						addX,
						addY);
				}
			}
		}
		oForce[valindex] = make_float2(addX, addY);
		oPower[valindex] = addX*addX+addY*addY;
	}
}
