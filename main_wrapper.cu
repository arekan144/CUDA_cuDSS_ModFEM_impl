#include "main_wrapper.cuh"

#define cleanUpFun() \
{ \
    cudssMatrixDestroy(A_h); \
    cudssMatrixDestroy(x_h); \
    cudssMatrixDestroy(b_h); \
    cudssDataDestroy(handle, solverData); \
    cudssConfigDestroy(solverConfig); \
    cudssDestroy(handle); \
    cudaFree(colptr_d); \
    cudaFree(rowptr_d); \
    cudaFree(values_d); \
    cudaFree(b_d); \
    cudaFree(x_d); \
} 

#ifdef _DEBUG
#include <iostream>
#define interpretCudaStatus(call, status, msg) \
status = call; \
if (status != cudaSuccess) { \
    std::cout << "Error: CUDA API returned on msg { " << msg << " }, details:\n"; \
    std::cout << cudaGetErrorString(cudaGetLastError()); \
    std::cout << std::endl; \
    cleanUpFun(); \
    return -1; \
}

static std::string giveCudssString(cudssStatus_t status) { 
    switch (status)
    {
    case CUDSS_STATUS_SUCCESS: return "SUCCESS";
    case CUDSS_STATUS_NOT_INITIALIZED: return "NOT INITIALIZED";
    case CUDSS_STATUS_ALLOC_FAILED: return "ALLOC FAILED";
    case CUDSS_STATUS_INVALID_VALUE: return "INVALID VALUE";
    case CUDSS_STATUS_NOT_SUPPORTED: return "NOT SUPPORTED";
    case CUDSS_STATUS_EXECUTION_FAILED: return "EXECUTION FAILED";
    case CUDSS_STATUS_INTERNAL_ERROR: return "INTERNAL ERROR";
    default: return "UKNOWN STATUS";
    }
}

#define interpretCudssStatus(call, status, msg) \
status = call; \
if (status != CUDSS_STATUS_SUCCESS) { \
    std::cout << "Example FAILED: CUDA API returned error { " << msg << " }\n"; \
    std::cout << "Msg: " << giveCudssString(status) << std::endl; \
    cleanUpFun(); \
    return -1; \
} 

#else
#define interpretCudaStatus(call, status, msg) \
if(call != cudaSuccess){ \
    cleanUpFun(); \
    return -1; \
}

#define interpretCudssStatus(call, status, msg) \
if(call != CUDSS_STATUS_SUCCESS){ \
    cleanUpFun(); \
    return -1; \
}
#endif


int doDecomposition(SparseStructures::CSR& matrix, double* b, double** x)
{
    cudaError_t cudaStatus = cudaSuccess;
    cudssStatus_t cudssStatus = CUDSS_STATUS_SUCCESS;

    int* colptr_d = nullptr;
    int* rowptr_d = nullptr;
    double* values_d = nullptr;
    double* b_d = nullptr;
    double* x_d = nullptr;

    cudssHandle_t handle = nullptr;
    cudssConfig_t config = nullptr;
    cudssData_t data = nullptr;

    cudssMatrix_t A_h = nullptr;
    cudssMatrixType_t mtype = CUDSS_MTYPE_GENERAL;
    cudssMatrixViewType_t mview = CUDSS_MVIEW_FULL;
    cudssIndexBase_t base = CUDSS_BASE_ZERO;

    cudssMatrix_t x_h = nullptr;
    cudssMatrix_t b_h = nullptr;

    cudssConfig_t solverConfig = nullptr;
    cudssData_t solverData = nullptr;

    interpretCudaStatus(
        cudaSetDevice(0), cudaStatus, "cudaSetDevice"
    );

    interpretCudaStatus(
        cudaMalloc((&colptr_d), static_cast<size_t>(matrix.getNNZ()) * sizeof(int)),
        cudaStatus, "cudaMalloc colptr");
    interpretCudaStatus(
        cudaMemcpy(colptr_d,matrix.getColInd(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy colptr");

    interpretCudaStatus(
        cudaMalloc(&rowptr_d, static_cast<size_t>(matrix.getN()+1u) * sizeof(int)),
        cudaStatus, "cudaMalloc rowptr");
    interpretCudaStatus(
        cudaMemcpy(rowptr_d, matrix.getRowPtr(),
            static_cast<size_t>(matrix.getN() + 1u) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy rowptr");

    interpretCudaStatus(
        cudaMalloc(&values_d, static_cast<size_t>(matrix.getNNZ()) * sizeof(double)),
        cudaStatus, "cudaMalloc values");
    interpretCudaStatus(
        cudaMemcpy(values_d, matrix.getAcsr(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy values");
    
    interpretCudaStatus(
        cudaMalloc(&b_d, static_cast<size_t>(matrix.getN())*sizeof(double)),
        cudaStatus, "cudaMalloc b");
    interpretCudaStatus(
        cudaMemcpy(b_d, b,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy b");

    interpretCudaStatus(
        cudaMalloc(&x_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc x");

    interpretCudssStatus(
        cudssCreate(&handle),
        cudssStatus, "cudssCreate");

    /*
    * (7.1) ONLY valid IndexTypes: R_32I, R_64I !
    */
    interpretCudssStatus(
        cudssMatrixCreateCsr(&A_h,   // [out] Matrix handle
            matrix.getN(),           // Number of rows
            matrix.getN(),           // Number of columns
            matrix.getNNZ(),         // Number of non-zeros
            rowptr_d,                // Row start offsets
            nullptr,                 // Row end offsets [Usefull when batch-loading, otherwise null]
            colptr_d,                // Column indices of the matrix
            values_d,                // Values of the dense matrix
            CUDA_R_32I,              // Index type of the matrix
            CUDA_R_64F,              // Data type of the matrix
            mtype,                   // Matrix type of the matrix
            mview,                   // Matrix view of the matrix
            base),                   // Indexing base
        cudssStatus, "cudssMatrixCreateCsr");
    
    interpretCudssStatus(
        cudssMatrixCreateDn(&x_h, matrix.getN(), 1ll, matrix.getN(), x_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn x");
    interpretCudssStatus(
        cudssMatrixCreateDn(&b_h, matrix.getN(), 1ll, matrix.getN(), b_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn b");

    interpretCudssStatus(
        cudssConfigCreate(&solverConfig),
        cudssStatus, "cudssConfigCreate");
    interpretCudssStatus(
        cudssDataCreate(handle, &solverData),
        cudssStatus, "cudssDataCreate");

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_ANALYSIS,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute ANALYSIS");
    
    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_FACTORIZATION,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute FACTORIZATION");

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_FACTORIZATION,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute FACTORIZATION");
   
    interpretCudaStatus(cudaDeviceSynchronize(), cudaStatus, "cudaDeviceSynchronize");

    (*x) = new double[matrix.getN()];
    interpretCudaStatus(
        cudaMemcpy((*x), x_d,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyDeviceToHost),
        cudaStatus, "cudaMemcpy x");

    cleanUpFun();

	return 0;
}
