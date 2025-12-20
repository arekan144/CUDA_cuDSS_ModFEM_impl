#include "cudss_wrapper.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cudss.h>
#include <cusparse.h>

//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"

#define cleanCUDAFun() { \
    cudaFree(colptr_d); \
    cudaFree(rowptr_d); \
    cudaFree(values_d); \
    cudaFree(b_d); \
    cudaFree(x_d); \
}

#define cleanCUDSSFun() { \
    cudssMatrixDestroy(A_h); \
    cudssMatrixDestroy(x_h); \
    cudssMatrixDestroy(b_h); \
    cudssDataDestroy(handle, solverData); \
    cudssConfigDestroy(solverConfig); \
    cudssDestroy(handle); \
}

#define cleanImageFun() { \

#ifdef _DEBUG
#include <iostream>
#define interpretCudaStatus(call, status, msg, cleanFun1, cleanFun2) \
status = call; \
if (status != cudaSuccess) { \
    std::cout << "Error: CUDA API returned on msg { " << msg << " }, details:\n"; \
    std::cout << cudaGetErrorString(cudaGetLastError()); \
    std::cout << std::endl; \
    cleanFun1(); \
    cleanFun2(); \
    return 1; \
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

#define interpretCudssStatus(call, status, msg, cleanFun1, cleanFun2) \
status = call; \
if (status != CUDSS_STATUS_SUCCESS) { \
    std::cout << "Error: cuDSS API returned error { " << msg << " }\n"; \
    std::cout << "Msg: " << giveCudssString(status) << std::endl; \
    cleanFun1(); \
    cleanFun2(); \
    return 2; \
} 

#else
#define interpretCudaStatus(call, status, msg, cleanFun1, cleanFun2) \
if(call != cudaSuccess){ \
    cleanFun1(); \
    cleanFun2(); \
    return 1; \
}

#define interpretCudssStatus(call, status, msg, cleanFun1, cleanFun2) \
if(call != CUDSS_STATUS_SUCCESS){ \
    cleanFun1(); \
    cleanFun2(); \
    return 2; \
}
#endif

__global__ void __fillImage(int* ) {

}

#define nothing()

#define cleanImageData() { \
 cudaFree(data_d);\
};

#include <iostream>
int saveSparcityImage(int n, int nnz, int* rowptr_post, int* colind_post) {

    size_t image_dimn = 1024;
    if (n < 1024)
        image_dimn = static_cast<size_t>(n);

    typedef unsigned char byte;
    
    byte* data_d = nullptr;
    byte* data_h = nullptr;
    

    /*interpretCudaStatus(
        cudaMalloc((&data_d), sizeof(byte) * image_dimn * image_dimn),
        cudaStatus,
        "cudaMalloc :img: data_d",
        cleanImageData,
        nothing);*/
    /*
    interpretCudaStatus(
        cudaMalloc((&rowptr_cudss), row_size),
        cudaStatus,
        "cudaMalloc :img: data_d",
        cleanImageData,
        nothing);
    interpretCudaStatus(
        cudaMalloc((&colptr_cudss), col_size),
        cudaStatus,
        "cudaMalloc :img: data_d",
        cleanImageData,
        nothing);
        */

    /*size_t verify;
    interpretCudssStatus(
        cudssDataGet(handle,
            matrix_data,
            CUDSS_DATA_PERM_REORDER_ROW,
            &rowptr_cudss,
            row_size,
            &verify),
        cudssStatus,
        "cudssDataGet :img: perm_reorder_row",
        cleanImageData,
        nothing
    );

    if (verify != row_size)
    {
        cleanImageData();
        return -1;
    }
    
    interpretCudssStatus(
        cudssDataGet(handle,
            matrix_data,
            CUDSS_DATA_PERM_REORDER_COL,
            &colptr_cudss,
            col_size,
            &verify),
        cudssStatus,
        "cudssDataGet :img: perm_reorder_col",
        cleanImageData,
        nothing
    );
    if (verify != col_size)
    {
        cleanImageData();
        return -1;
    }
    
    int* rowptr_h = new int[static_cast<size_t>(n + 1)];
    int* colptr_h = new int[static_cast<size_t>(nnz)];

    cudaMemcpy(rowptr_h, rowptr_cudss, row_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(colptr_h, colptr_cudss, col_size, cudaMemcpyDeviceToHost);

    for (auto i = 0llu; i < static_cast<size_t>(5); i++)
    {
        auto ky = static_cast<size_t>(rowptr_h[i + 1llu] - rowptr_h[i]);
        for (auto j = 0llu; j < ky; j++)
        {
            std::cout << colptr_h[static_cast<size_t>(rowptr_h[i])] << ' ';
        }
        std::cout << '\n';
    }

    delete[] rowptr_h;
    delete[] colptr_h;
*/

    cleanImageData();
    return 0;
}
#include <stdio.h>
int cuDSSOnlyAnalisysAndSpPattern(SparseStructures::CSR& matrix,
    double* b, double** x, short matrix_type, short view_type, short index_base) {
    cudaError_t cudaStatus = cudaSuccess;
    cudssStatus_t cudssStatus = CUDSS_STATUS_SUCCESS;

    int* colptr_d = nullptr;
    int* rowptr_d = nullptr;
    double* values_d = nullptr;
    double* b_d = nullptr;
    double* x_d = nullptr;

    cudssHandle_t handle = nullptr;

    cudssMatrix_t A_h = nullptr;
    cudssMatrixType_t mtype = static_cast<cudssMatrixType_t>(matrix_type); // info for types can be found cudss.h, or cudss_wrapper.h
    cudssMatrixViewType_t mview = static_cast<cudssMatrixViewType_t>(view_type);
    cudssIndexBase_t base = static_cast<cudssIndexBase_t>(index_base);

    cudssMatrix_t x_h = nullptr;
    cudssMatrix_t b_h = nullptr;

    cudssConfig_t solverConfig = nullptr;
    cudssData_t solverData = nullptr;

    interpretCudaStatus(
        cudaSetDevice(0),
        cudaStatus, "cudaSetDevice",
        cleanCUDAFun, cleanCUDSSFun);

    cudaStream_t stream;
    interpretCudaStatus(
        cudaStreamCreate(&stream),
        cudaStatus, "cudaStreamCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssCreate(&handle),
        cudssStatus, "cudssCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssSetStream(handle, stream),
        cudssStatus, "cudssSetStream",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc((&colptr_d), static_cast<size_t>(matrix.getNNZ()) * sizeof(int)),
        cudaStatus, "cudaMalloc colptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&rowptr_d, static_cast<size_t>(matrix.getN() + 1u) * sizeof(int)),
        cudaStatus, "cudaMalloc rowptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&values_d, static_cast<size_t>(matrix.getNNZ()) * sizeof(double)),
        cudaStatus, "cudaMalloc values",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&b_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc b",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&x_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc x",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMemcpy(colptr_d, matrix.getColInd(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy colptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(rowptr_d, matrix.getRowOff(),
            static_cast<size_t>(matrix.getN() + 1u) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy rowptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(values_d, matrix.getAcsr(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy values",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(b_d, b,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy b",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssMatrixCreateCsr(&A_h,   // [out] Matrix handle
            matrix.getN(),           // Number of rows [host]
            matrix.getN(),           // Number of columns [host]
            matrix.getNNZ(),         // Number of non-zeros [host]
            rowptr_d,                // Row start offsets [device]
            nullptr,                 // Row end offsets [Usefull when batch-loading, otherwise null] 
            colptr_d,                // Column indices of the matrix [device]
            values_d,                // Values of the dense matrix (nonzeros) [device]
            CUDA_R_32I,              // Index type of the matrix [(7.1) ONLY valid IndexTypes: R_32I, R_64I !] [host]
            CUDA_R_64F,              // Data type of the matrix [host]
            mtype,                   // Type of the matrix: HERMITIAN, GENERAL, SYMMETRIC [host]
            mview,                   // View of the matrix: FULL, UPPER, LOWER [host]
            base),                   // Indexing base: ZERO, ONE [host]
        cudssStatus, "cudssMatrixCreateCsr",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssMatrixCreateDn(&x_h, matrix.getN(), 1ll, matrix.getN(), x_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn x",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudssStatus(
        cudssMatrixCreateDn(&b_h, matrix.getN(), 1ll, matrix.getN(), b_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn b",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssConfigCreate(&solverConfig),
        cudssStatus, "cudssConfigCreate",
        cleanCUDAFun, cleanCUDSSFun);
    cudssAlgType_t reorder_alg = CUDSS_ALG_DEFAULT;
    interpretCudssStatus(
    cudssConfigSet(solverConfig, CUDSS_CONFIG_REORDERING_ALG,
        &reorder_alg, sizeof(cudssAlgType_t)),
        cudssStatus, "cudssConfigSet REORDERING_ALG",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssDataCreate(handle, &solverData),
        cudssStatus, "cudssDataCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_ANALYSIS,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute ANALYSIS",
        cleanCUDAFun, cleanCUDSSFun);

    cudaDeviceSynchronize();

    size_t sizeWritten;
    std::cout << matrix.getN() << "\n";
    size_t perm_s = matrix.getN() * sizeof(int);
    int* perm_col = nullptr;
 
    std::cout << "COL\n";
    interpretCudaStatus(
        cudaMalloc(&perm_col, perm_s),
        cudaStatus, "cudaMalloc x",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(cudssDataGet(handle, solverData, CUDSS_DATA_PERM_REORDER_COL, perm_col,
        perm_s, &sizeWritten), cudssStatus, "cudssDataGet for reorder row perm",
        cleanCUDSSFun, cleanCUDAFun);

    /*int* perm_row = nullptr;
    std::cout << "ROW\n";
    interpretCudaStatus(
        cudaMalloc(&perm_col, perm_s),
        cudaStatus, "cudaMalloc x",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(cudssDataGet(handle, solverData, CUDSS_DATA_PERM_REORDER_ROW, perm_row,
        perm_s, &sizeWritten), cudssStatus, "cudssDataGet for reorder row perm",
        cleanCUDSSFun, cleanCUDAFun); // not needed, sizeWritten returns 0*/

    //saveSparcityImage(matrix.getN(), matrix.getNNZ(), handle, solverData);

    cudaFree(perm_col);

    cleanCUDSSFun();
    cleanCUDAFun();
    return 0;
}

int cuSparseDecomposition(SparseStructures::CSR& matrix,
    double* b, double** x,  short matrix_type, short view_type, short index_base)
{

    return 0;
}

int cuDSSDecomposition(SparseStructures::CSR& matrix,
    double* b, double** x,  short matrix_type,  short view_type, short index_base) {
    cudaError_t cudaStatus = cudaSuccess;
    cudssStatus_t cudssStatus = CUDSS_STATUS_SUCCESS;

    int* colptr_d = nullptr;
    int* rowptr_d = nullptr;
    double* values_d = nullptr;
    double* b_d = nullptr;
    double* x_d = nullptr;

    cudssHandle_t handle = nullptr;

    cudssMatrix_t A_h = nullptr;
    cudssMatrixType_t mtype = static_cast<cudssMatrixType_t>(matrix_type); // info for types can be found cudss.h, or cudss_wrapper.h
    cudssMatrixViewType_t mview = static_cast<cudssMatrixViewType_t>(view_type);
    cudssIndexBase_t base = static_cast<cudssIndexBase_t>(index_base);

    cudssMatrix_t x_h = nullptr;
    cudssMatrix_t b_h = nullptr;

    cudssConfig_t solverConfig = nullptr;
    cudssData_t solverData = nullptr;

    interpretCudaStatus(
        cudaSetDevice(0), 
        cudaStatus, "cudaSetDevice",
        cleanCUDAFun, cleanCUDSSFun);

    cudaStream_t stream;
    interpretCudaStatus(
        cudaStreamCreate(&stream), 
        cudaStatus, "cudaStreamCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssCreate(&handle),
        cudssStatus, "cudssCreate",
        cleanCUDAFun, cleanCUDSSFun);
    
    interpretCudssStatus(
        cudssSetStream(handle, stream),
        cudssStatus, "cudssSetStream",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc((&colptr_d), static_cast<size_t>(matrix.getNNZ()) * sizeof(int)),
        cudaStatus, "cudaMalloc colptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&rowptr_d, static_cast<size_t>(matrix.getN() + 1u) * sizeof(int)),
        cudaStatus, "cudaMalloc rowptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&values_d, static_cast<size_t>(matrix.getNNZ()) * sizeof(double)),
        cudaStatus, "cudaMalloc values",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&b_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc b",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMalloc(&x_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc x",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMemcpy(colptr_d, matrix.getColInd(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy colptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(rowptr_d, matrix.getRowOff(),
            static_cast<size_t>(matrix.getN() + 1u) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy rowptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(values_d, matrix.getAcsr(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy values",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(b_d, b,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy b",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssMatrixCreateCsr(&A_h,   // [out] Matrix handle
            matrix.getN(),           // Number of rows
            matrix.getN(),           // Number of columns
            matrix.getNNZ(),         // Number of non-zeros
            rowptr_d,                // Row start offsets
            nullptr,                 // Row end offsets [Usefull when batch-loading, otherwise null]
            colptr_d,                // Column indices of the matrix
            values_d,                // Values of the dense matrix (nonzeros)
            CUDA_R_32I,              // Index type of the matrix [(7.1) ONLY valid IndexTypes: R_32I, R_64I !]
            CUDA_R_64F,              // Data type of the matrix
            mtype,                   // Type of the matrix: HERMITIAN, GENERAL, SYMMETRIC
            mview,                   // View of the matrix: FULL, UPPER, LOWER
            base),                   // Indexing base: ZERO, ONE
        cudssStatus, "cudssMatrixCreateCsr",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssMatrixCreateDn(&x_h, matrix.getN(), 1ll, matrix.getN(), x_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn x",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudssStatus(
        cudssMatrixCreateDn(&b_h, matrix.getN(), 1ll, matrix.getN(), b_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn b",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssConfigCreate(&solverConfig),
        cudssStatus, "cudssConfigCreate",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudssStatus(
        cudssDataCreate(handle, &solverData),
        cudssStatus, "cudssDataCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_ANALYSIS,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute ANALYSIS",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_FACTORIZATION,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute FACTORIZATION",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_SOLVE,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute SOLVE",
        cleanCUDAFun, cleanCUDSSFun);
    cleanCUDSSFun();
    
    interpretCudaStatus(
        cudaStreamSynchronize(stream),
        cudaStatus, "cudaStreamSynchronize",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaStreamDestroy(stream),
        cudaStatus, "cudaStreamDestroy",
        cleanCUDAFun, cleanCUDSSFun);

    (*x) = new double[matrix.getN()];
    interpretCudaStatus(
        cudaMemcpy((*x), x_d,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyDeviceToHost),
        cudaStatus, "cudaMemcpy x",
        cleanCUDAFun, cleanCUDSSFun);
    
    cleanCUDAFun();
    return 0;
}

int cuDSSDecompositionWithSynchronization(ChronoTimer& timer, SparseStructures::CSR& matrix,
    double* b, double** x, short matrix_type,  short view_type, short index_base)
{
    timer.setStartTime(); // time for starting library, setting device, allocating memory and copying
    cudaError_t cudaStatus = cudaSuccess;
    cudssStatus_t cudssStatus = CUDSS_STATUS_SUCCESS;

    int* colptr_d = nullptr;
    int* rowptr_d = nullptr;
    double* values_d = nullptr;
    double* b_d = nullptr;
    double* x_d = nullptr;

    cudssHandle_t handle = nullptr;

    cudssMatrix_t A_h = nullptr;
    cudssMatrixType_t mtype = static_cast<cudssMatrixType_t>(matrix_type); // info for types can be found cudss.h, or cudss_wrapper.h
    cudssMatrixViewType_t mview = static_cast<cudssMatrixViewType_t>(view_type);
    cudssIndexBase_t base = static_cast<cudssIndexBase_t>(index_base);

    cudssMatrix_t x_h = nullptr;
    cudssMatrix_t b_h = nullptr;

    cudssConfig_t solverConfig = nullptr;
    cudssData_t solverData = nullptr;

    interpretCudaStatus(
        cudaSetDevice(0),
        cudaStatus, "cudaSetDevice",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssCreate(&handle),
        cudssStatus, "cudssCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc((&colptr_d), static_cast<size_t>(matrix.getNNZ()) * sizeof(int)),
        cudaStatus, "cudaMalloc colptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(colptr_d, matrix.getColInd(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy colptr",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc(&rowptr_d, static_cast<size_t>(matrix.getN() + 1u) * sizeof(int)),
        cudaStatus, "cudaMalloc rowptr",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(rowptr_d, matrix.getRowOff(),
            static_cast<size_t>(matrix.getN() + 1u) * sizeof(int), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy rowptr",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc(&values_d, static_cast<size_t>(matrix.getNNZ()) * sizeof(double)),
        cudaStatus, "cudaMalloc values",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(values_d, matrix.getAcsr(),
            static_cast<size_t>(matrix.getNNZ()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy values",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc(&b_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc b",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudaStatus(
        cudaMemcpy(b_d, b,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyHostToDevice),
        cudaStatus, "cudaMemcpy b",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(
        cudaMalloc(&x_d, static_cast<size_t>(matrix.getN()) * sizeof(double)),
        cudaStatus, "cudaMalloc x",
        cleanCUDAFun, cleanCUDSSFun);

    timer.saveTimeNow();

    timer.setStartTime(); // time for creation of csr matrix, dense X and B matrixes

    interpretCudssStatus(
        cudssMatrixCreateCsr(&A_h,   // [out] Matrix handle
            matrix.getN(),           // Number of rows
            matrix.getN(),           // Number of columns
            matrix.getNNZ(),         // Number of non-zeros
            rowptr_d,                // Row start offsets
            nullptr,                 // Row end offsets [Usefull when batch-loading, otherwise null]
            colptr_d,                // Column indices of the matrix
            values_d,                // Values of the dense matrix (nonzeros)
            CUDA_R_32I,              // Index type of the matrix [(7.1) ONLY valid IndexTypes: R_32I, R_64I !]
            CUDA_R_64F,              // Data type of the matrix
            mtype,                   // Type of the matrix: HERMITIAN, GENERAL, SYMMETRIC
            mview,                   // View of the matrix: FULL, UPPER, LOWER
            base),                   // Indexing base: ZERO, ONE
        cudssStatus, "cudssMatrixCreateCsr",
        cleanCUDAFun, cleanCUDSSFun);


    interpretCudssStatus(
        cudssMatrixCreateDn(&x_h, matrix.getN(), 1ll, matrix.getN(), x_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn x",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudssStatus(
        cudssMatrixCreateDn(&b_h, matrix.getN(), 1ll, matrix.getN(), b_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
        cudssStatus, "cudssMatrixCreateDn b",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudssStatus(
        cudssConfigCreate(&solverConfig),
        cudssStatus, "cudssConfigCreate",
        cleanCUDAFun, cleanCUDSSFun);
    interpretCudssStatus(
        cudssDataCreate(handle, &solverData),
        cudssStatus, "cudssDataCreate",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(cudaDeviceSynchronize(), cudaStatus, "cudaDeviceSynchronize",
        cleanCUDAFun, cleanCUDSSFun);
    timer.saveTimeNow();
    timer.setStartTime(); // time for ANALYSIS

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_ANALYSIS,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute ANALYSIS",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(cudaDeviceSynchronize(), cudaStatus, "cudaDeviceSynchronize",
        cleanCUDAFun, cleanCUDSSFun);
    timer.saveTimeNow();
    timer.setStartTime(); // time for FACTORIZATION

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_FACTORIZATION,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute FACTORIZATION",
        cleanCUDAFun, cleanCUDSSFun);

    interpretCudaStatus(cudaDeviceSynchronize(), cudaStatus, "cudaDeviceSynchronize",
        cleanCUDAFun, cleanCUDSSFun);
    timer.saveTimeNow();
    timer.setStartTime(); // time for SOLVE

    interpretCudssStatus(
        cudssExecute(handle,
            CUDSS_PHASE_SOLVE,
            solverConfig,
            solverData,
            A_h,
            x_h,
            b_h),
        cudssStatus, "cudssExecute SOLVE",
        cleanCUDAFun, cleanCUDSSFun);


    interpretCudaStatus(cudaDeviceSynchronize(), cudaStatus, "cudaDeviceSynchronize",
        cleanCUDAFun, cleanCUDSSFun);
    timer.saveTimeNow();
    timer.setStartTime(); // time for copying the solution


    (*x) = new double[matrix.getN()];
    interpretCudaStatus(
        cudaMemcpy((*x), x_d,
            static_cast<size_t>(matrix.getN()) * sizeof(double), cudaMemcpyDeviceToHost),
        cudaStatus, "cudaMemcpy x",
        cleanCUDAFun, cleanCUDSSFun);
    timer.saveTimeNow();
    cleanCUDSSFun();
    cleanCUDAFun();

    return 0;
}
