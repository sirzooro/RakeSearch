// Search for pairs of diagonal Latin squares by the method of rows permutation

#include "MovePairSearch.h"
#include "boinc_api.h"

#ifdef __SSE2__
#include "immintrin.h"
#endif
#ifdef __ARM_NEON
#include "arm_neon.h"
#endif

#ifndef _MSC_VER
#define ffs __builtin_ffs
#else
#include <intrin.h>
inline int ffs(int i)
{
  unsigned long index;
  if (_BitScanForward(&index, i))
    return index + 1;
  else
    return 0;
}
#endif

#define SetBit(bitfield, bitno) ((bitfield) |= (1u << (bitno)))
#define ClearBit(bitfield, bitno) ((bitfield) &= ~(1u << (bitno)))

#define AllBitsMask(numbits) ((1u << (numbits)) - 1)

inline void MovePairSearch::CopyRow(int* __restrict dst, int* __restrict src)
{
  int n = 0;
#ifdef __AVX__
  for (; n < Rank - 7; n += 8)
  {
    __m256i v = _mm256_loadu_si256((const __m256i*)&src[n]);
    _mm256_storeu_si256((__m256i*)&dst[n], v);
  }
#endif
#ifdef __SSE2__
  for (; n < Rank - 3; n += 4)
  {
    __m128i v = _mm_loadu_si128((const __m128i*)&src[n]);
    _mm_storeu_si128((__m128i*)&dst[n], v);
  }
#endif
  for (; n < Rank; n++)
  {
    dst[n] = src[n];
  }
}

inline void MovePairSearch::SetRow(int* dst, int val)
{
  int n = 0;
#ifdef __AVX__
  __m256i v_avx = _mm256_set1_epi32(val);
  for (; n < Rank - 7; n += 8)
  {
    _mm256_storeu_si256((__m256i*)&dst[n], v_avx);
  }
#endif
#ifdef __SSE2__
#ifdef __AVX__
  __m128i v_sse = _mm256_castsi256_si128(v_avx);
#else
  __m128i v_sse = _mm_set1_epi32(val);
#endif
  for (; n < Rank - 3; n += 4)
  {
    _mm_storeu_si128((__m128i*)&dst[n], v_sse);
  }
#endif
  for (; n < Rank; n++)
  {
    dst[n] = val;
  }
}

// Default constructor
MovePairSearch::MovePairSearch()
{
  Reset();
}

// Reset search settings
void MovePairSearch::Reset()
{
  // Reset the DLS generator settings
  squareAGenerator.Reset();

  // Reset the next search settings
  ClearBeforeNextSearch();

  // Reset the values of global counters
  totalPairsCount = 0;
  totalSquaresWithPairs = 0;
  totalProcessedSquaresSmall = 0;
  totalProcessedSquaresLarge = 0;

  startParametersFileName = "start_parameters.txt";
  resultFileName = "result.txt";
  checkpointFileName = "checkpoint.txt";
  tempCheckpointFileName = "tmp_checkpoint.txt";

  // Set the header constant for the file of parameters or checkpoint
  moveSearchGlobalHeader = "# Move search of pairs OLDS status";
  moveSearchComponentHeader = "# Move search component status";

  // Reset the initialization flag
  isInitialized = 0;
}


// Reset the variables before the next search
void MovePairSearch::ClearBeforeNextSearch()
{
  // Reset the values of matrices of squares A and B,
  // also the matrix of rows usage when forming square B
  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      squareA[i][j] = -1;
      squareB[i][j] = -1;
    }
    rowsHistory[i] = AllBitsMask(Rank);
  }

  // Reset the values in the vectors of rows usage in the next permutation
  // and the rows numbers used for the current square
  for (int i = 0; i < Rank; i++)
  {
    currentSquareRows[i] = -1;
  }

  // Reset the counter of pairs found for the given DLS
  pairsCount = 0;
}


// Search initialization
void MovePairSearch::InitializeMoveSearch(string start, string result,
                                       string checkpoint, string temp)
{
  fstream startFile;
  fstream checkpointFile;

  // Read the filenames
  startParametersFileName = start;
  resultFileName = result;
  checkpointFileName = checkpoint;
  tempCheckpointFileName = temp;

  // Read the generator state and the search state from the checkpoint file or initial values
    // Open files with initial parameters and the checkpoint file
    startFile.open(startParametersFileName.c_str(), std::ios_base::in);
    checkpointFile.open(checkpointFileName.c_str(), std::ios_base::in);

    // Read the state
	if (checkpointFile.is_open())
	{
		// Read the state from the checkpoint file
		try
		{
			isStartFromCheckpoint = 1;
			Read(checkpointFile);
		}
		catch (...)
		{
			cerr << "Error opening checkpoint file! Starting with workunit start parameters." << endl;
			isStartFromCheckpoint = 0;
		}
    }
    else
    {
        isStartFromCheckpoint = 0;
    }

	if (isStartFromCheckpoint != 1)
    {
		// Read the state from the initial parameters file
		Read(startFile);
		isStartFromCheckpoint = 0;
    }

    // Close the files
    startFile.close();
    checkpointFile.close();
}


// Read the search state from stream
void MovePairSearch::Read(istream& is)
{
  string marker;

  // Reset the initialization flag
  isInitialized = 0;

  // Read the search state
    // Find the start marker of the state
    do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
    }
    while (marker != moveSearchGlobalHeader);

    // Read the state of the DLS generator
    is >> squareAGenerator;

    // Find the marker of permutation component
	do
	{
		std::getline(is, marker);

		if (is.eof())
		{
			throw ("Expected start marker, but EOF found.");
		}
	}
    while (marker != moveSearchComponentHeader);

    // Read the search variables by permutation (actually, the variables containing statistics)
    is >> pairsCount;
    is >> totalPairsCount;
    is >> totalSquaresWithPairs;
    is >> totalProcessedSquaresLarge;
    is >> totalProcessedSquaresSmall;

  // Set the initialization flag
  isInitialized = 1;
}


// Write the search state into stream
void MovePairSearch::Write(ostream& os)
{
  // Write the search state
    // Writing the header
    os << moveSearchGlobalHeader << endl;
    os << endl;

    // Write the state of the DLS generator
    os << squareAGenerator;
    os << endl;

    // Write the header of permutation block
    os << moveSearchComponentHeader << endl;
    os << endl;

    // Writing statistical values
    os << pairsCount << " " << totalPairsCount << " " << totalSquaresWithPairs << endl;
    os << totalProcessedSquaresLarge << " " << totalProcessedSquaresSmall << endl;
    os << endl;
}


// Create a checkpoint
void MovePairSearch::CreateCheckpoint()
{
  ofstream checkpointFile;

  checkpointFile.open(tempCheckpointFileName.c_str(), std::ios_base::out);
  if (checkpointFile.is_open())
  {
    Write(checkpointFile);
    checkpointFile.close();
    remove(checkpointFileName.c_str());
    rename(tempCheckpointFileName.c_str(), checkpointFileName.c_str());
  }
  else
    cerr << "Error opening checkpoint file!" << endl;
}


// Start the search for orthogonal squares by the method of rows permutation
void MovePairSearch::StartMoveSearch()
{
  // Subscribe to the event of finding a DLS
  squareAGenerator.Subscribe(this);

  // Start DLS generation
  squareAGenerator.Start();

  // Unsubscribe from the event of finding a DLS
  squareAGenerator.Unsubscribe();

  // Output the search results
  ShowSearchTotals();
}


#if defined(__ARM_NEON) && !defined(__aarch64__)
__attribute__((always_inline))
inline void MovePairSearch::transposeMatrix4x4(int srcRow, int srcCol, int destRow, int destCol)
{
  uint16x4_t v1, v2;
  v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+1][srcCol+0]));
  v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+1][srcCol+2]));
  uint16x4_t v1_1 = vuzp_u16(v1, v2).val[0];
  v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+2][srcCol+0]));
  v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+2][srcCol+2]));
  uint16x4_t v2_1 = vuzp_u16(v1, v2).val[0];
  v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+3][srcCol+0]));
  v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+3][srcCol+2]));
  uint16x4_t v3_1 = vuzp_u16(v1, v2).val[0];
  v1 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+4][srcCol+0]));
  v2 = vld1_u16((uint16_t*)(&squareA_Mask[srcRow+4][srcCol+2]));
  uint16x4_t v4_1 = vuzp_u16(v1, v2).val[0];
  
  uint16x4x2_t v12_2 = vtrn_u16(v1_1, v2_1);
  uint16x4x2_t v34_2 = vtrn_u16(v3_1, v4_1);

  uint32x2x2_t v13_3 = vtrn_u32(vreinterpret_u32_u16(v12_2.val[0]), vreinterpret_u32_u16(v34_2.val[0]));
  uint32x2x2_t v24_3 = vtrn_u32(vreinterpret_u32_u16(v12_2.val[1]), vreinterpret_u32_u16(v34_2.val[1]));
  
  vst1_u32((uint32_t*)(&squareA_MaskT[destRow+0][destCol+0]), v13_3.val[0]);
  vst1_u32((uint32_t*)(&squareA_MaskT[destRow+1][destCol+0]), v24_3.val[0]);
  vst1_u32((uint32_t*)(&squareA_MaskT[destRow+2][destCol+0]), v13_3.val[1]);
  vst1_u32((uint32_t*)(&squareA_MaskT[destRow+3][destCol+0]), v24_3.val[1]);
}
#endif


// Event processor of DLS generation, will start the search for its pair
void MovePairSearch::OnSquareGenerated(const Square& newSquare)
{
  // Reset before the search for orthogonal squares
  ClearBeforeNextSearch();

  // Write the found square

  // Copy square and generate masks first
#ifdef __AVX2__
  // AVX2 has "shift by vector" instruction, use it here
  // Note: AVX512 instructions which use ZMM registers cause too big
  // CPU frequency throttling. It does not make sense to use them in this
  // one place only, as everything else will be slowed down too.
  int n = 0;
  for (; n < Rank*Rank-7; n += 8)
  {
    // Copy data
    __m256i v = _mm256_load_si256 ((__m256i*)(&newSquare.Matrix[0][0] + n));
    _mm256_store_si256((__m256i*)(&squareA[0][0] + n), v);

    // Calculate bitmasks
    v = _mm256_sllv_epi32(_mm256_set1_epi32(1), v);
    _mm256_store_si256((__m256i*)(&squareA_Mask[0][0] + n), v);
  }
  // Process remaining elements
  for (; n < Rank*Rank; n++)
  {
    int x = *(&newSquare.Matrix[0][0] + n);
    *((&squareA[0][0] + n)) = x;
    *((&squareA_Mask[0][0] + n)) = 1 << x;
  }
#elif defined(__SSSE3__)
  // SSSE3 added shuffle instruction, which can be used to build small lookup table.
  // Maximum val (256) needs 9 bytes. Fortunately all values can be decreased by one,
  // so all of them will fit in 8 bits here.
  const __m128i vcLut = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 255, 127, 63, 31, 15, 7, 3, 1, 0);
  const __m128i vc1 = _mm_set1_epi16(1);
  const __m128i vc0 = _mm_setzero_si128();
  int n = 0;
  for (; n < Rank*Rank-7; n += 8)
  {
    // Copy data
    __m128i v1 = _mm_load_si128((__m128i*)(&newSquare.Matrix[0][0] + n));
    __m128i v2 = _mm_load_si128((__m128i*)(&newSquare.Matrix[0][0] + n + 4));
    _mm_store_si128((__m128i*)(&squareA[0][0] + n), v1);
    _mm_store_si128((__m128i*)(&squareA[0][0] + n + 4), v2);

    // Pack two 32x4 vectors into one 8x16
    v1 = _mm_packs_epi32(v1, v2);
    v1 = _mm_packs_epi16(v1, vc0);
    // Get maak fom LUT
    v1 = _mm_shuffle_epi8(vcLut, v1);
    // Convert vector 8x16 to 16x8, and increment values by one
    v1 = _mm_unpacklo_epi8(v1, vc0);
    v1 = _mm_add_epi16(v1, vc1);

    // Convert vector 16x8 into two 32x4, and store results
    v2 = _mm_unpackhi_epi16(v1, vc0);
    v1 = _mm_unpacklo_epi16(v1, vc0);

    _mm_store_si128((__m128i*)(&squareA_Mask[0][0] + n), v1);
    _mm_store_si128((__m128i*)(&squareA_Mask[0][0] + n + 4), v2);
  }
  // Process remaining elements
  for (; n < Rank*Rank; n++)
  {
    int x = *(&newSquare.Matrix[0][0] + n);
    *((&squareA[0][0] + n)) = x;
    *((&squareA_Mask[0][0] + n)) = 1 << x;
  }
#else
  // Default non-SIMD code
  // Note: this will be autovectorized for ARM NEON. gcc has limit how many times
  // it can unroll loop, so single loop with manual vectorization would be slower.
  // Two nested loops are below limit, so autovectorization creates expected
  // machine code here.
  for (int i = 0; i < Rank; i++)
  {
    for (int j = 0; j < Rank; j++)
    {
      squareA[i][j] = newSquare.Matrix[i][j];
      squareA_Mask[i][j] = 1u << newSquare.Matrix[i][j];
    }
  }
#endif

  // Create transposed copy of squareA_Mask if needed
#ifdef __SSE2__
  __m128i v1, v2;
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[1][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[1][4]));
  __m128i v1_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[2][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[2][4]));
  __m128i v2_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[3][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[3][4]));
  __m128i v3_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[4][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[4][4]));
  __m128i v4_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[5][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[5][4]));
  __m128i v5_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[6][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[6][4]));
  __m128i v6_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[7][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[7][4]));
  __m128i v7_1 = _mm_packs_epi32(v1, v2);
  v1 = _mm_loadu_si128((__m128i*)(&squareA_Mask[8][0]));
  v2 = _mm_loadu_si128((__m128i*)(&squareA_Mask[8][4]));
  __m128i v8_1 = _mm_packs_epi32(v1, v2);

  __m128i v1_2 = _mm_unpacklo_epi16(v1_1, v2_1);
  __m128i v2_2 = _mm_unpackhi_epi16(v1_1, v2_1);
  __m128i v3_2 = _mm_unpacklo_epi16(v3_1, v4_1);
  __m128i v4_2 = _mm_unpackhi_epi16(v3_1, v4_1);
  __m128i v5_2 = _mm_unpacklo_epi16(v5_1, v6_1);
  __m128i v6_2 = _mm_unpackhi_epi16(v5_1, v6_1);
  __m128i v7_2 = _mm_unpacklo_epi16(v7_1, v8_1);
  __m128i v8_2 = _mm_unpackhi_epi16(v7_1, v8_1);

  __m128i v1_3 = _mm_unpacklo_epi32(v1_2, v3_2);
  __m128i v2_3 = _mm_unpackhi_epi32(v1_2, v3_2);
  __m128i v3_3 = _mm_unpacklo_epi32(v2_2, v4_2);
  __m128i v4_3 = _mm_unpackhi_epi32(v2_2, v4_2);
  __m128i v5_3 = _mm_unpacklo_epi32(v5_2, v7_2);
  __m128i v6_3 = _mm_unpackhi_epi32(v5_2, v7_2);
  __m128i v7_3 = _mm_unpacklo_epi32(v6_2, v8_2);
  __m128i v8_3 = _mm_unpackhi_epi32(v6_2, v8_2);

  __m128i v1_4 = _mm_unpacklo_epi64(v1_3, v5_3);
  __m128i v2_4 = _mm_unpackhi_epi64(v1_3, v5_3);
  __m128i v3_4 = _mm_unpacklo_epi64(v2_3, v6_3);
  __m128i v4_4 = _mm_unpackhi_epi64(v2_3, v6_3);
  __m128i v5_4 = _mm_unpacklo_epi64(v3_3, v7_3);
  __m128i v6_4 = _mm_unpackhi_epi64(v3_3, v7_3);
  __m128i v7_4 = _mm_unpacklo_epi64(v4_3, v8_3);
  __m128i v8_4 = _mm_unpackhi_epi64(v4_3, v8_3);

  _mm_store_si128((__m128i*)(&squareA_MaskT[0][0]), v1_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[1][0]), v2_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[2][0]), v3_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[3][0]), v4_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[4][0]), v5_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[5][0]), v6_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[6][0]), v7_4);
  _mm_store_si128((__m128i*)(&squareA_MaskT[7][0]), v8_4);

  for (int i = 1; i < Rank; i++)
  {
    for (int j = 8; j < Rank; j++)
    {
      squareA_MaskT[j][i-1] = squareA_Mask[i][j];
    }
  }
#elif defined(__ARM_NEON)
#ifdef __aarch64__
  uint16x8_t v1, v2;
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[1][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[1][4]));
  uint16x8_t v1_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[2][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[2][4]));
  uint16x8_t v2_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[3][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[3][4]));
  uint16x8_t v3_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[4][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[4][4]));
  uint16x8_t v4_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[5][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[5][4]));
  uint16x8_t v5_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[6][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[6][4]));
  uint16x8_t v6_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[7][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[7][4]));
  uint16x8_t v7_1 = vuzp1q_u16(v1, v2);
  v1 = vld1q_u16((uint16_t*)(&squareA_Mask[8][0]));
  v2 = vld1q_u16((uint16_t*)(&squareA_Mask[8][4]));
  uint16x8_t v8_1 = vuzp1q_u16(v1, v2);
  
  uint16x8_t v1_2 = vtrn1q_u16(v1_1, v2_1);
  uint16x8_t v2_2 = vtrn2q_u16(v1_1, v2_1);
  uint16x8_t v3_2 = vtrn1q_u16(v3_1, v4_1);
  uint16x8_t v4_2 = vtrn2q_u16(v3_1, v4_1);
  uint16x8_t v5_2 = vtrn1q_u16(v5_1, v6_1);
  uint16x8_t v6_2 = vtrn2q_u16(v5_1, v6_1);
  uint16x8_t v7_2 = vtrn1q_u16(v7_1, v8_1);
  uint16x8_t v8_2 = vtrn2q_u16(v7_1, v8_1);

  uint32x4_t v1_3 = vtrn1q_u32(vreinterpretq_u32_u16(v1_2), vreinterpretq_u32_u16(v3_2));
  uint32x4_t v2_3 = vtrn1q_u32(vreinterpretq_u32_u16(v2_2), vreinterpretq_u32_u16(v4_2));
  uint32x4_t v3_3 = vtrn2q_u32(vreinterpretq_u32_u16(v1_2), vreinterpretq_u32_u16(v3_2));
  uint32x4_t v4_3 = vtrn2q_u32(vreinterpretq_u32_u16(v2_2), vreinterpretq_u32_u16(v4_2));
  uint32x4_t v5_3 = vtrn1q_u32(vreinterpretq_u32_u16(v5_2), vreinterpretq_u32_u16(v7_2));
  uint32x4_t v6_3 = vtrn1q_u32(vreinterpretq_u32_u16(v6_2), vreinterpretq_u32_u16(v8_2));
  uint32x4_t v7_3 = vtrn2q_u32(vreinterpretq_u32_u16(v5_2), vreinterpretq_u32_u16(v7_2));
  uint32x4_t v8_3 = vtrn2q_u32(vreinterpretq_u32_u16(v6_2), vreinterpretq_u32_u16(v8_2));
  
  uint64x2_t v1_4 = vtrn1q_u64(vreinterpretq_u64_u32(v1_3), vreinterpretq_u64_u32(v5_3));
  uint64x2_t v2_4 = vtrn1q_u64(vreinterpretq_u64_u32(v2_3), vreinterpretq_u64_u32(v6_3));
  uint64x2_t v3_4 = vtrn1q_u64(vreinterpretq_u64_u32(v3_3), vreinterpretq_u64_u32(v7_3));
  uint64x2_t v4_4 = vtrn1q_u64(vreinterpretq_u64_u32(v4_3), vreinterpretq_u64_u32(v8_3));
  uint64x2_t v5_4 = vtrn2q_u64(vreinterpretq_u64_u32(v1_3), vreinterpretq_u64_u32(v5_3));
  uint64x2_t v6_4 = vtrn2q_u64(vreinterpretq_u64_u32(v2_3), vreinterpretq_u64_u32(v6_3));
  uint64x2_t v7_4 = vtrn2q_u64(vreinterpretq_u64_u32(v3_3), vreinterpretq_u64_u32(v7_3));
  uint64x2_t v8_4 = vtrn2q_u64(vreinterpretq_u64_u32(v4_3), vreinterpretq_u64_u32(v8_3));
  
  vst1q_u64((uint64_t*)(&squareA_MaskT[0][0]), v1_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[1][0]), v2_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[2][0]), v3_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[3][0]), v4_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[4][0]), v5_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[5][0]), v6_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[6][0]), v7_4);
  vst1q_u64((uint64_t*)(&squareA_MaskT[7][0]), v8_4);
  
  for (int i = 1; i < Rank; i++)
  {
    for (int j = 8; j < Rank; j++)
    {
      squareA_MaskT[j][i-1] = squareA_Mask[i][j];
    }
  }
#else // !__aarch64__
  transposeMatrix4x4(0, 0, 0, 0);
  transposeMatrix4x4(0, 4, 4, 0);
  transposeMatrix4x4(4, 0, 0, 4);
  transposeMatrix4x4(4, 4, 4, 4);
  
  for (int i = 1; i < Rank; i++)
  {
    for (int j = 8; j < Rank; j++)
    {
      squareA_MaskT[j][i-1] = squareA_Mask[i][j];
    }
  }
#endif // !__aarch64__
#else
  // Non-SIMD code does not use squareA_MaskT, so nothing here
#endif // __SSE2__

  // Start the rows permutation
  MoveRows();

  // Check the mutual orthogonality of the squares
  if (pairsCount > 0)
  {
    CheckMutualOrthogonality();
  }

  // Gather statistics on the processed squares
  totalProcessedSquaresSmall++;

  // Fix the information about state of processing
  if (totalProcessedSquaresSmall % CheckpointInterval == 0)
  {
    // Update the completion progress for the BOINC client
    double fraction_done;
    if(Rank == 8)
      fraction_done = (double)(totalProcessedSquaresSmall)/5000000.0;
    else
    {
      if(Rank == 9)
        fraction_done = (double)(totalProcessedSquaresSmall)/275000000.0;
      else
        fraction_done = 0.999999999;
    }

    if(fraction_done >=1) fraction_done = (double)(totalProcessedSquaresSmall)/300000000.0;
    if(fraction_done >=1) fraction_done = 0.999999999;

    boinc_fraction_done(fraction_done); // Tell the BOINC client the completion fraction

    // Check if the BOINC client is able to create a checkpoint,
    // and if yes, start the function of its creation
    if (boinc_time_to_checkpoint()) {
      CreateCheckpoint();
      boinc_checkpoint_completed(); // BOINC client knows the checkpoint has been created
    }

    if(isDebug)
    {
      cout << "# ------------------------" << endl;
      cout << "# Processed " << totalProcessedSquaresLarge << " milliards and " << totalProcessedSquaresSmall << " squares." << endl;
      cout << "# Last processed square:" << endl;
      cout << endl;
      cout << newSquare;
      cout << "# ------------------------" << endl;
    }
  }
}

#define DBG_UP()
#define DBG_DOWN()

// Permute the rows of the given DLS, trying to find ODLS for it
void MovePairSearch::MoveRows()
{
  int currentRowId;
  int gettingRowId = -1;

  int diagonalValues1, diagonalValues2;

  int diagonalValuesHistory[Rank][2];

  int rowsUsage; // Flags of the rows usage at the current moment; rowsUsage[number of the row] = 0 | 1, where 0 means the row is already used, 1 - not.

  int rowCandidates; // Rows which still have to be checked

  // Write the 1st row of square A into square B for the search of normalized squares
  CopyRow(&squareB[0][0], &squareA[0][0]);

  // Mark the usage of the 1st row, because it is fixed
  rowsUsage = AllBitsMask(Rank) & ~1u;
  ClearBit(rowsHistory[0], 0);
  currentSquareRows[0] = 0;

#if !defined (__SSE2__) && !defined(__ARM_NEON)
  // For non-vectorized builds start from 1st row, and mark all rows except 0th as candidates
  currentRowId = 1;
  rowCandidates = AllBitsMask(Rank) & ~1u;
#else
  // SSE2/AVX2 version performs duplicate check when it steps forward to next row
  // and saves result as a candidates. So we have to start from 0th row and
  // mark it as the only candidate in order to perform duplicate check for 1st row.
  currentRowId = 0;
  rowCandidates = 1;
#endif

  // Set bits for diagonal values in 1st row
  diagonalValuesHistory[0][0] = 1u << squareB[0][0];
  diagonalValuesHistory[0][1] = 1u << squareB[0][Rank - 1];

  diagonalValues1 = diagonalValuesHistory[0][0];
  diagonalValues2 = diagonalValuesHistory[0][1];

#ifdef __ARM_NEON
  // Set the powers of 2
  const uint16_t powersOf2[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };
#ifdef __aarch64__
  const uint16x8_t vPowersOf2 = vld1q_u16(powersOf2);
#else
  const uint16x4_t vPowersOf2Lo = vld1_u16(powersOf2);
  const uint16x4_t vPowersOf2Hi = vld1_u16(powersOf2+4);
#endif
#endif

  // Note: code below assumes that rowCandidates is non-zero at beginning
  // so 1st nested while loop should execute first. If it may not be the case,
  // change code to handle this.
  // 1st loop (used to be "if (rowCandidates)" part) - handle case when at least one row candidate is present
  while (1)
  {
    // Select a row from the initial square for the position currentRowId of the generated square
    // Process the search result
    while(1)
    {
      gettingRowId = __builtin_ctz(rowCandidates);
      // Process the new found row

      // Check diagonality of the generated part of the square
      // Check the main diagonal and secondary diagonal
      // Get bits for current row
      int bit1 = squareA_Mask[gettingRowId][currentRowId];
      int bit2 = squareA_Mask[gettingRowId][Rank - 1 - currentRowId];

#if !defined(__SSE2__) && !defined(__ARM_NEON)
      // Duplicate check
      int duplicationDetected = (diagonalValues1 & bit1) | (diagonalValues2 & bit2);

        // Process the results of checking the square for diagonality
        if (duplicationDetected)
        {
          // Mark the row in the history of the used rows
          ClearBit(rowCandidates, gettingRowId);
          
          if (0 == rowCandidates)
            break;
          else
            continue;
        }
#endif

        // Mark the row in the history of the used rows
        ClearBit(rowCandidates, gettingRowId);

        {
          // Write the row into the array of the current rows
          currentSquareRows[currentRowId] = gettingRowId;

          // Step forward depending on the current position
          if (currentRowId == Rank - 1)
          {
            // Write rows into the square in correct order
            for (int n = 1; n < Rank; ++n)
            {
              CopyRow(&squareB[n][0], &squareA[currentSquareRows[n]][0]);
            }

            // Process the found square
            ProcessOrthoSquare();
            break;
          }
          else
          {
            // Save new bitmasks for diagonal values for further use
            diagonalValues1 |= bit1;
            diagonalValues2 |= bit2;
            diagonalValuesHistory[currentRowId][0] = diagonalValues1;
            diagonalValuesHistory[currentRowId][1] = diagonalValues2;

            // Save remaining candidates in row history
            rowsHistory[currentRowId] = rowCandidates;

            // Mark the row in the array of the used rows
            ClearBit(rowsUsage, gettingRowId);

#if !defined(__SSE2__) && !defined(__ARM_NEON)
            // Set new row candidates
            rowCandidates = rowsUsage;
#endif

            // Step forward
            currentRowId++;
            DBG_UP();

#ifdef __SSE2__
            // load bitmasks for columns which will be on diagonals
            // for performance reasons load this as a row from transposed square
            // also excluse 0th element, row 0 has fixed position in square
            __m128i vCol1 = _mm_load_si128((const __m128i*)&squareA_MaskT[currentRowId][0]);
            __m128i vCol2 = _mm_load_si128((const __m128i*)&squareA_MaskT[Rank - 1 - currentRowId][0]);

            // AND loaded values with diagnonal masks
            __m128i vDiagMask1 = _mm_set1_epi16(diagonalValues1);
            __m128i vDiagMask2 = _mm_set1_epi16(diagonalValues2);

            vCol1 = _mm_and_si128(vCol1, vDiagMask1);
            vCol2 = _mm_and_si128(vCol2, vDiagMask2);

            // non-zero means that number is duplicated, zero means that it is unique
            // OR these values together first
            vCol1 = _mm_or_si128(vCol1, vCol2);

#if defined(__AVX512F__) && defined(__AVX512VL__)
            // check if result is zero and get result as a bitmask
            __mmask8 resultMask = _mm_testn_epi16_mask(vCol1, vCol1);

            // add one bit for 0th row, and AND result with rowsUsage
            rowCandidates = (resultMask << 1) & rowsUsage;
#else // !AVX512
            // check if result is zero
            vCol1 = _mm_cmpeq_epi16(vCol1, _mm_setzero_si128());

            // create mask from vector
            // there are 2 bits per result, so we need to pack int16 to int8 first
            vCol1 = _mm_packs_epi16(vCol1, _mm_setzero_si128());
            unsigned int mask = _mm_movemask_epi8(vCol1);

            // add one bit for 0th row, and AND result with rowsUsage
            rowCandidates = (mask << 1) & rowsUsage;
#endif // !AVX512
#elif defined(__ARM_NEON)
#ifdef __aarch64__
            // load bitmasks for columns which will be on diagonals
            // for performance reasons load this as a row from transposed square
            // also excluse 0th element, row 0 has fixed position in square
            uint16x8_t vCol1 = vld1q_u16((const uint16_t*)&squareA_MaskT[currentRowId][0]);
            uint16x8_t vCol2 = vld1q_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][0]);

            // AND loaded values with diagnonal masks
            uint16x8_t vDiagMask1 = vdupq_n_u16(diagonalValues1);
            uint16x8_t vDiagMask2 = vdupq_n_u16(diagonalValues2);

            vCol1 = vandq_u16(vCol1, vDiagMask1);
            vCol2 = vandq_u16(vCol2, vDiagMask2);

            // non-zero means that number is duplicated, zero means that it is unique
            // OR these values together first
            vCol1 = vorrq_u16(vCol1, vCol2);

            // check if result is zero
            vCol1 = vceqq_u16(vCol1, vdupq_n_u16(0));

            // create mask from vector
            uint16x8_t v = vandq_u16(vCol1, vPowersOf2);
            uint32_t mask = vaddvq_u64(vpaddlq_u32(vpaddlq_u16(v)));

            // add one bit for 0th row, and AND result with rowsUsage
            rowCandidates = (mask << 1) & rowsUsage;
#else // !__aarch64__
            // load bitmasks for columns which will be on diagonals
            // for performance reasons load this as a row from transposed square
            // also excluse 0th element, row 0 has fixed position in square
            uint16x4_t vCol1a = vld1_u16((const uint16_t*)&squareA_MaskT[currentRowId][0]);
            uint16x4_t vCol1b = vld1_u16((const uint16_t*)&squareA_MaskT[currentRowId][4]);
            
            uint16x4_t vCol2a = vld1_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][0]);
            uint16x4_t vCol2b = vld1_u16((const uint16_t*)&squareA_MaskT[Rank - 1 - currentRowId][4]);

            // AND loaded values with diagnonal masks
            uint16x4_t vDiagMask1 = vdup_n_u16(diagonalValues1);
            uint16x4_t vDiagMask2 = vdup_n_u16(diagonalValues2);

            vCol1a = vand_u16(vCol1a, vDiagMask1);
            vCol1b = vand_u16(vCol1b, vDiagMask1);
            
            vCol2a = vand_u16(vCol2a, vDiagMask2);
            vCol2b = vand_u16(vCol2b, vDiagMask2);

            // non-zero means that number is duplicated, zero means that it is unique
            // OR these values together first
            vCol1a = vorr_u16(vCol1a, vCol2a);
            vCol1b = vorr_u16(vCol1b, vCol2b);

            // check if result is zero
            vCol1a = vceq_u16(vCol1a, vdup_n_u16(0));
            vCol1b = vceq_u16(vCol1b, vdup_n_u16(0));

            // create mask from vector
            vCol1a = vand_u16(vCol1a, vPowersOf2Lo);
            vCol1b = vand_u16(vCol1b, vPowersOf2Hi);
            
            vCol1a = vorr_u16(vCol1a, vCol1b);
            
            uint64x1_t v = vpaddl_u32(vpaddl_u16(vCol1a));
            uint32_t mask = vget_lane_u32(vreinterpret_u32_u64(v), 0);

            // add one bit for 0th row, and AND result with rowsUsage
            rowCandidates = (mask << 1) & rowsUsage;
#endif
#endif // AVX2/SSE2
            if (!rowCandidates)
                break;
          }
        }
    }

    // 2nd loop (used to be "else" part) - handle case when there are no row candidates
    while(1)
    {
      // Process not-founding of the new row: step backward, clear the flags of usage,
      // the history of usage, the list of current rows and clear the square itself

        // Step backward
        currentRowId--;
        DBG_DOWN();
        // Check if we are done
        if (0 == currentRowId)
          return;
        // Get saved values for previous row
        diagonalValues1 = diagonalValuesHistory[currentRowId-1][0];
        diagonalValues2 = diagonalValuesHistory[currentRowId-1][1];
        // Clear the flag of row usage
        SetBit(rowsUsage, currentSquareRows[currentRowId]);
        // Get saved candidates
        rowCandidates = rowsHistory[currentRowId];
        if (rowCandidates)
            break;
    }
  }
}


// Process the found orthogonal square
void MovePairSearch::ProcessOrthoSquare()
{
  int isDifferent = 0;      // The number of differences in the rows with the initial square (to avoid generating its copy)

  Square a(squareA);        // Square A as an object
  Square b(squareB);        // Square B as an object

  int orthoMetric = Rank*Rank;  // The value of orthogonality metric saying that the squares are fully orthogonal

  // Check the square for being a copy of the initial one
  isDifferent = 0;

  for (int i = 0; i < Rank; i++)
  {
    if (currentSquareRows[i] != i)
    {
      isDifferent = 1;
      break;
    }
  }

  // Process the found square
  if (isDifferent && Square::OrthoDegree(a, b) == orthoMetric
      && b.IsDiagonal() && b.IsLatin() && a.IsDiagonal() && b.IsLatin())
  {
    // Write the information about the found square
      // Increase the squares counter
      pairsCount++;
      totalPairsCount++;

      // Remember the basic square
      if (pairsCount == 1)
      {
        orthoSquares[pairsCount - 1] = a;
        totalSquaresWithPairs++;
      }

      // Remember the pair square
      if (pairsCount < OrhoSquaresCacheSize)
      {
        orthoSquares[pairsCount] = b;
      }

      // The stream for output into the results file
      // It must be here, creation and destruction of iostream is costly!
      ofstream resultFile;
      resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
      if (!resultFile.is_open())
      {
        std::cerr << "Error opening file!";
      }

      // Output the header
      if (pairsCount == 1)
      {
        if(isDebug)
        {
          // Output the information about the 1st square in a pair as a header
          cout << "{" << endl;
          cout << "# ------------------------" << endl;
          cout << "# Detected pair for the square: " << endl;
          cout << "# ------------------------" << endl;
          cout << a;
          cout << "# ------------------------" << endl;
        }
        // Write the information into file
        if (resultFile.is_open())
        {
          resultFile << "{" << endl;
          resultFile << "# ------------------------" << endl;
          resultFile << "# Detected pair for the square: " << endl;
          resultFile << "# ------------------------" << endl;
          resultFile << a;
          resultFile << "# ------------------------" << endl;
        }
      }

      // Output the information about the found pair
        if(isDebug)
        {
          // Output the information into the console
          cout << b << endl;
        }

        // Output the information into the file
        if (resultFile.is_open())
        {
          resultFile << b << endl;
          resultFile.close();
        }
  }
}


// Check the mutual orthogonality of a set of squares found in the current search
void MovePairSearch::CheckMutualOrthogonality()
{
  int orthoMetric = Rank*Rank;
  int maxSquareId;
  ofstream resultFile;

  // Detect the upper bound of the squares to process
  if (pairsCount < OrhoSquaresCacheSize)
  {
    maxSquareId = pairsCount;
  }
  else
  {
    maxSquareId = OrhoSquaresCacheSize - 1;
  }

  // Open the file with results
  resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
  if (!resultFile.is_open()) { cout << "Error opening file!"; return; }

  // Check the mutual orthogonality of a set of squares
  for (int i = 0; i <= maxSquareId; i++)
  {
    for (int j = i + 1; j <= maxSquareId; j++)
    {
      if (Square::OrthoDegree(orthoSquares[i], orthoSquares[j]) == orthoMetric)
      {
        if(isDebug) cout << "# Square " << i << " # " << j << endl;
        resultFile << "# Square " << i << " # " << j << endl;
      }
    }
  }
  if(isDebug) cout << endl;
  resultFile << endl;

  // Write the total number of the found DLS pairs
  if(isDebug) cout << "# Pairs found: " << pairsCount << endl;
  resultFile << "# Pairs found: " << pairsCount << endl;

  // Set the end marker of the results section
  if(isDebug) cout << "}" << endl;
  resultFile << "}" << endl;

  // Close the file with results
  resultFile.close();
}

// Display the total results of the search
void MovePairSearch::ShowSearchTotals()
{
  ofstream resultFile;

  if(isDebug)
  {
    // Write the results into the console
    cout << "# ------------------------" << endl;
    cout << "# Total pairs found: " << totalPairsCount << endl;
    cout << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
    cout << "# ------------------------" << endl;
  }

  // Write the results into the file
  resultFile.open(resultFileName.c_str(), std::ios_base::binary | std::ios_base::app);
  if (resultFile.is_open())
  {
    resultFile << "# ------------------------" << endl;
    resultFile << "# Total pairs found: " << totalPairsCount << endl;
    resultFile << "# Total squares with pairs: " << totalSquaresWithPairs << endl;
    resultFile << "# Processes " << totalProcessedSquaresLarge << " milliards " << totalProcessedSquaresSmall << " squares" << endl;
    resultFile << "# ------------------------" << endl;
    resultFile.close();
  }
  else cerr << "Error opening file!" << endl;
}
