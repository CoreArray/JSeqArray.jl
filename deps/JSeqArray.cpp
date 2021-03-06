// ===========================================================
//
// SeqArray.cpp: the C/C++ codes for the JSeqArray package
//
// Copyright (C) 2017    Xiuwen Zheng
//
// This file is part of JSeqArray.
//
// JSeqArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// JSeqArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with JSeqArray.
// If not, see <http://www.gnu.org/licenses/>.

#include "Index.h"

#include <set>
#include <algorithm>

#include "ReadByVariant.h"
// #include "ReadBySample.h"
#include <ctype.h>



// ===========================================================
// Library Functions
// ===========================================================

extern "C"
{

using namespace CoreArray;
using namespace JSeqArray;


// ===========================================================
// Open a GDS file
// ===========================================================

/// initialize a SeqArray file
JL_DLLEXPORT void SEQ_File_Init(int file_id)
{
	COREARRAY_TRY
		CFileInfo &file = GetFileInfo(file_id);
		file.Selection();  // force to initialize selection
	COREARRAY_CATCH
}

/// finalize a SeqArray file
JL_DLLEXPORT void SEQ_File_Done(int file_id)
{
	COREARRAY_TRY
		map<int, CFileInfo>::iterator p = GDSFile_ID_Info.find(file_id);
		if (p != GDSFile_ID_Info.end())
			GDSFile_ID_Info.erase(p);
	COREARRAY_CATCH
}



// ===========================================================
// Progress bar
// ===========================================================

JL_DLLEXPORT void *SEQ_ProgressInit(C_Int64 count, C_BOOL verbose)
{
	return new CProgressStdOut(count, verbose);
}

JL_DLLEXPORT void SEQ_ProgressForward(void *ptr)
{
	if (ptr)
	{
		CProgressStdOut *obj = (CProgressStdOut*)ptr;
		obj->Forward();
	}
}

JL_DLLEXPORT void SEQ_ProgressDone(void *ptr)
{
	if (ptr)
	{
		CProgressStdOut *obj = (CProgressStdOut*)ptr;
		delete obj;
	}
}




// ===========================================================
// Set a working space
// ===========================================================

/// push the current filter to the stack
JL_DLLEXPORT void SEQ_FilterPush(int file_id, C_BOOL new_flag)
{
	COREARRAY_TRY
		map<int, CFileInfo>::iterator it = GDSFile_ID_Info.find(file_id);
		if (it != GDSFile_ID_Info.end())
		{
			if (new_flag || it->second.SelList.empty())
				it->second.SelList.push_back(TSelection());
			else
				it->second.SelList.push_back(it->second.SelList.back());
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// pop up the previous filter from the stack
JL_DLLEXPORT void SEQ_FilterPop(int file_id)
{
	COREARRAY_TRY
		map<int, CFileInfo>::iterator it = GDSFile_ID_Info.find(file_id);
		if (it != GDSFile_ID_Info.end())
		{
			if (it->second.SelList.size() <= 1)
				throw ErrSeqArray("No filter can be pop up.");
			it->second.SelList.pop_back();
		} else
			throw ErrSeqArray("The GDS file is closed or invalid.");
	COREARRAY_CATCH
}


/// set a working space with selected sample id with a logical vector
JL_DLLEXPORT void SEQ_SetSample(int file_id,  jl_array_t *samp_sel,
	C_BOOL intersect, C_BOOL verbose)
{
	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(file_id);
		TSelection &Sel = File.Selection();
		C_BOOL *pArray = Sel.pSample();
		size_t Count = File.SampleNum();

		if (!jl_is_nothing(samp_sel))
		{
			if (!intersect)
			{
				if (jl_array_len(samp_sel) != Count)
					throw ErrSeqArray("Invalid length of 'sample'.");
				memcpy(pArray, (void*)jl_array_data(samp_sel), Count);
			} else {
				if (jl_array_len(samp_sel) != File.SampleSelNum())
				{
					throw ErrSeqArray(
						"Invalid length of 'sample' "
						"(should be equal to the number of selected samples).");
				}
				// set selection
				C_BOOL *p = (C_BOOL*)jl_array_data(samp_sel);
				for (size_t i=0; i < Count; i++)
				{
					if (*pArray)
						*pArray = ((*p++) != 0);
					pArray ++;
				}
			}
		} else {
			// reset the filter
			memset(pArray, 1, Count);
		}

		int n = File.SampleSelNum();
		if (verbose)
			jl_printf(JL_STDOUT, "Number of selected samples: %s\n", PrettyInt(n));

	COREARRAY_CATCH
}


/// set a working space with selected variant id with a logical vector
JL_DLLEXPORT void SEQ_SetVariant(int file_id,  jl_array_t *var_sel,
	C_BOOL intersect, C_BOOL verbose)
{
	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(file_id);
		TSelection &Sel = File.Selection();
		C_BOOL *pArray = Sel.pVariant();
		size_t Count = File.VariantNum();

		if (!jl_is_nothing(var_sel))
		{
			if (!intersect)
			{
				if (jl_array_len(var_sel) != Count)
					throw ErrSeqArray("Invalid length of 'variant'.");
				memcpy(pArray, (void*)jl_array_data(var_sel), Count);
			} else {
				if (jl_array_len(var_sel) != File.VariantSelNum())
				{
					throw ErrSeqArray(
						"Invalid length of 'variant' "
						"(should be equal to the number of selected variants).");
				}
				// set selection
				C_BOOL *p = (C_BOOL*)jl_array_data(var_sel);
				for (size_t i=0; i < Count; i++)
				{
					if (*pArray)
						*pArray = ((*p++) != 0);
					pArray ++;
				}
			}
		} else {
			// reset the filter
			memset(pArray, 1, Count);
		}

		int n = File.VariantSelNum();
		if (verbose)
			jl_printf(JL_STDOUT, "Number of selected variants: %s\n", PrettyInt(n));

	COREARRAY_CATCH
}


/*
// ================================================================

static bool is_numeric(const string &txt)
{
	char *endptr = (char*)(txt.c_str());
	strtol(txt.c_str(), &endptr, 10);
	return (endptr != txt.c_str()) && (*endptr == 0);
}

/// set a working space flag with selected chromosome(s)
JL_DLLEXPORT PyObject* SEQ_SetChrom(PyObject* gdsfile, PyObject* include,
	PyObject* is_num, PyObject* frombp, PyObject* tobp, PyObject* intersect, PyObject* verbose)
{
	int nProtected = 0;
	int *pFrom=NULL, *pTo=NULL;

	int IsNum = Rf_asLogical(is_num);
	int IsIntersect = Rf_asLogical(intersect);
	if (IsIntersect == NA_INTEGER)
		error("'intersect' should be either FALSE or TRUE.");

	if (Rf_isNull(include))
	{
		if (!Rf_isNull(frombp))
			error("'from.bp' should be NULL.");
		if (!Rf_isNull(tobp))
			error("'to.bp' should be NULL.");
	} else {
		include = PROTECT(AS_CHARACTER(include));
		nProtected ++;
		if (!Rf_isNull(frombp) || !Rf_isNull(tobp))
		{
			if (RLength(include) != RLength(frombp))
				error("'from.bp' should have the same length as 'include'.");
			if (RLength(include) != RLength(tobp))
				error("'to.bp' should have the same length as 'include'.");
			frombp = PROTECT(AS_INTEGER(frombp));
			tobp = PROTECT(AS_INTEGER(tobp));
			pFrom = INTEGER(frombp); pTo = INTEGER(tobp);
			nProtected += 2;
		}
	}

	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();

		vector<C_BOOL> &sel_array = Sel.Variant;
		vector<C_BOOL> tmp_array;
		if (IsIntersect) tmp_array.resize(sel_array.size());

		vector<C_BOOL> &array = IsIntersect ? tmp_array : sel_array;
		memset(&array[0], FALSE, array.size());

		if (Rf_isNull(include))
		{
			// include = NULL
			if (IsNum == NA_INTEGER)
			{
				memset(&array[0], TRUE, array.size());
			} else {
				CChromIndex &Chrom = File.Chromosome();
				map<string, CChromIndex::TRangeList>::iterator it;
				for (it=Chrom.Map.begin(); it != Chrom.Map.end(); it++)
				{
					bool flag = is_numeric(it->first);
					if (((IsNum==TRUE) && flag) || ((IsNum==FALSE) && !flag))
					{
						CChromIndex::TRangeList &rng = it->second;
						vector<CChromIndex::TRange>::iterator it;
						for (it=rng.begin(); it != rng.end(); it++)
						{
							memset(&array[it->Start], TRUE, it->Length);
						}
					}
				}
			}

		} else {
			// include != NULL
			vector<C_Int32> *varPos = NULL;
			if (pFrom && pTo)
				varPos = &File.Position();

			CChromIndex &Chrom = File.Chromosome();
			map<string, CRangeSet> RngSets;

			R_xlen_t n = XLENGTH(include);
			for (R_xlen_t idx=0; idx < n; idx++)
			{
				string s = CHAR(STRING_ELT(include, idx));

				if (IsNum == TRUE)
				{
					if (!is_numeric(s)) continue;
				} else if (IsNum == FALSE)
				{
					if (is_numeric(s)) continue;
				}

				map<string, CChromIndex::TRangeList>::iterator it =
					Chrom.Map.find(s);
				if (it != Chrom.Map.end())
				{
					if (varPos)
					{
						// if from.bp and to.bp
						int from = pFrom[idx], to = pTo[idx];
						if (from == NA_INTEGER) from = 0;
						if (to == NA_INTEGER) to = 2147483647;
						RngSets[s].AddRange(from, to);
					} else {
						// no from.bp and to.bp
						CChromIndex::TRangeList &rng = it->second;
						vector<CChromIndex::TRange>::iterator p;
						for (p=rng.begin(); p != rng.end(); p++)
						{
							memset(&array[p->Start], TRUE, p->Length);
						}
					}
				}
			}

			if (varPos)
			{
				map<string, CRangeSet>::iterator it;
				for (it=RngSets.begin(); it != RngSets.end(); it++)
				{
					CChromIndex::TRangeList &rng = Chrom.Map[it->first];
					CRangeSet &RngSet = it->second;
					vector<CChromIndex::TRange>::const_iterator p;
					for (p=rng.begin(); p != rng.end(); p++)
					{
						size_t i = p->Start;
						size_t n = p->Length;
						C_Int32 *s = &((*varPos)[0]) + i;
						if (!IsIntersect)
						{
							for (; n > 0; n--, i++)
								if (RngSet.IsIncluded(*s++)) array[i] = TRUE;
						} else {
							C_BOOL *b = &sel_array[i];
							for (; n > 0; n--, i++, s++)
							{
								if (*b++)
									if (RngSet.IsIncluded(*s)) array[i] = TRUE;
							}
						}
					}
				}
			}
		}

		if (IsIntersect)
		{
			C_BOOL *p = &sel_array[0];
			C_BOOL *s = &array[0];
			for (size_t n=sel_array.size(); n > 0; n--)
				(*p++) &= (*s++);
		}

		if (Rf_asLogical(verbose) == TRUE)
		{
			int n = GetNumOfTRUE(&sel_array[0], sel_array.size());
			Rprintf("# of selected variants: %s\n", PrettyInt(n));
		}

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}
*/

// ================================================================

/// get the sample/variant filter
JL_DLLEXPORT jl_array_t* SEQ_GetFilter(int file_id, C_BOOL sample)
{
	jl_array_t *rv_ans = NULL;
	COREARRAY_TRY

		CFileInfo &File = GetFileInfo(file_id);
		TSelection &Sel = File.Selection();

		jl_value_t *atype = jl_apply_array_type(jl_bool_type, 1);
		if (sample)
		{
			size_t n = File.SampleNum();
			rv_ans = jl_alloc_array_1d(atype, n);
			memcpy(jl_array_data(rv_ans), Sel.pSample(), n);
		} else {
			size_t n = File.VariantNum();
			rv_ans = jl_alloc_array_1d(atype, n);
			memcpy(jl_array_data(rv_ans), Sel.pVariant(), n);
		}

	COREARRAY_CATCH
	return rv_ans;
}


// ================================================================

/// get the total number of samples
JL_DLLEXPORT C_Int64 SEQ_Attr_NSamp(int file_id)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		return File.SampleNum();
	COREARRAY_CATCH_RET
}

/// get the number of selected samples
JL_DLLEXPORT C_Int64 SEQ_Attr_NSelSamp(int file_id)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		return File.SampleSelNum();
	COREARRAY_CATCH_RET
}

/// get the total number of variants
JL_DLLEXPORT C_Int64 SEQ_Attr_NVar(int file_id)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		return File.VariantNum();
	COREARRAY_CATCH_RET
}

/// get the number of selected variants
JL_DLLEXPORT C_Int64 SEQ_Attr_NSelVar(int file_id)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		return File.VariantSelNum();
	COREARRAY_CATCH_RET
}

/// get the number of selected variants
JL_DLLEXPORT int SEQ_Attr_Ploidy(int file_id)
{
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		return File.Ploidy();
	COREARRAY_CATCH_RET
}



/// get the dimensions of the whole space
JL_DLLEXPORT jl_array_t *SEQ_GetSpace(int file_id)
{
	jl_array_t *rv_ans = NULL;
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		jl_value_t *atype = jl_apply_array_type(jl_int64_type, 1);
		rv_ans = jl_alloc_array_1d(atype, 3);
		C_Int64 *p = (C_Int64*)jl_array_data(rv_ans);
		p[0] = File.Ploidy();
		p[1] = File.SampleNum();
		p[2] = File.VariantNum();
	COREARRAY_CATCH
	return rv_ans;
}


/// get the dimensions of selected space
JL_DLLEXPORT jl_array_t *SEQ_GetSelSpace(int file_id)
{
	jl_array_t *rv_ans = NULL;
	COREARRAY_TRY
		CFileInfo &File = GetFileInfo(file_id);
		jl_value_t *atype = jl_apply_array_type(jl_int64_type, 1);
		rv_ans = jl_alloc_array_1d(atype, 3);
		C_Int64 *p = (C_Int64*)jl_array_data(rv_ans);
		p[0] = File.Ploidy();
		p[1] = File.SampleSelNum();
		p[2] = File.VariantSelNum();
	COREARRAY_CATCH
	return rv_ans;
}


/*
/// set a working space with selected variant id
JL_DLLEXPORT PyObject* SEQ_Summary(PyObject* gdsfile, PyObject* varname)
{
	COREARRAY_TRY

		// the selection
		CFileInfo &File = GetFileInfo(gdsfile);
		TSelection &Sel = File.Selection();
		// the GDS root node
		PdGDSFolder Root = GDS_R_SEXP2FileRoot(gdsfile);
		// the variable name
		string vn = CHAR(STRING_ELT(varname, 0));

		if ((vn=="genotype") || (vn=="phase"))
		{
			PdGDSObj vGeno = GDS_Node_Path(Root, "genotype/data", TRUE);
			if (vGeno == NULL)
			{
				vGeno = GDS_Node_Path(Root, "genotype/~data", FALSE);
				if (vGeno == NULL)
				{
					throw ErrSeqArray(
						"There is no 'genotype/data' or 'genotype/~data'.");
				}
			}

			PROTECT(rv_ans = NEW_LIST(2));

				PyObject* I32 = PROTECT(NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 0, I32);
				C_Int32 Buf[4];
				GDS_Array_GetDim(vGeno, Buf, 3);
				INTEGER(I32)[0] = Buf[2];
				INTEGER(I32)[1] = Sel.Sample.size();
				INTEGER(I32)[2] = Sel.Variant.size();

				PyObject* S32 = PROTECT(NEW_INTEGER(3));
				SET_ELEMENT(rv_ans, 1, S32);
				INTEGER(S32)[0] = Buf[2];
				INTEGER(S32)[1] = GetNumOfTRUE(&Sel.Sample[0], Sel.Sample.size());
				INTEGER(S32)[2] = GetNumOfTRUE(&Sel.Variant[0], Sel.Variant.size());

			PyObject* tmp = PROTECT(NEW_CHARACTER(2));
				SET_STRING_ELT(tmp, 0, mkChar("dim"));
				SET_STRING_ELT(tmp, 1, mkChar("seldim"));
				SET_NAMES(rv_ans, tmp);
			UNPROTECT(4);

		} else {
			PdGDSObj var = GDS_Node_Path(Root, vn.c_str(), TRUE);
			rv_ans = ScalarInteger(GDS_Array_GetTotalCount(var));
		}

	COREARRAY_CATCH
}


/// get a logical vector with selection
JL_DLLEXPORT PyObject* SEQ_SelectFlag(PyObject* select, PyObject* len)
{
	R_len_t n = XLENGTH(select);
	if (XLENGTH(len) != n)
		error("Index variable error.");

	int *p = INTEGER(len);
	R_len_t m = 0;
	for (R_len_t k=n; k > 0; k--, p++)
	{
		if (*p > 0) m += *p;
	}

	PyObject* rv_ans = NEW_LOGICAL(m);
	int *r = INTEGER(rv_ans), *s = INTEGER(select);
	p = INTEGER(len);
	for (; n > 0; n--, s++, p++)
	{
		for (int k=*p; k > 0; k--)
			*r++ = *s;
	}

	return rv_ans;
}


// ===========================================================
// get system configuration
// ===========================================================

JL_DLLEXPORT PyObject* SEQ_IntAssign(PyObject* Dst, PyObject* Src)
{
	INTEGER(Dst)[0] = Rf_asInteger(Src);
	return R_NilValue;
}



// ===========================================================
// get system configuration
// ===========================================================

/// the number of alleles per site
JL_DLLEXPORT PyObject* SEQ_System()
{
	COREARRAY_TRY

		int nProtect = 0;
		rv_ans = PROTECT(NEW_LIST(2));
		PyObject* nm = PROTECT(NEW_CHARACTER(2));
		nProtect += 2;
		SET_NAMES(rv_ans, nm);

		// the number of logical cores
		SET_ELEMENT(rv_ans, 0, ScalarInteger(GDS_Mach_GetNumOfCores()));
		SET_STRING_ELT(nm, 0, mkChar("num.logical.core"));

		// compiler flags
		vector<string> ss;

	#ifdef COREARRAY_SIMD_SSE
		ss.push_back("SSE");
	#endif
	#ifdef COREARRAY_SIMD_SSE2
		ss.push_back("SSE2");
	#endif
	#ifdef COREARRAY_SIMD_SSE3
		ss.push_back("SSE3");
	#endif
	#ifdef COREARRAY_SIMD_SSSE3
		ss.push_back("SSSE3");
	#endif
	#ifdef COREARRAY_SIMD_SSE4_1
		ss.push_back("SSE4.1");
	#endif
	#ifdef COREARRAY_SIMD_SSE4_2
		ss.push_back("SSE4.2");
	#endif
	#ifdef COREARRAY_SIMD_AVX
		ss.push_back("AVX");
	#endif
	#ifdef COREARRAY_SIMD_AVX2
		ss.push_back("AVX2");
	#endif
	#ifdef COREARRAY_SIMD_FMA
		ss.push_back("FMA");
	#endif
	#ifdef COREARRAY_SIMD_FMA4
		ss.push_back("FMA4");
	#endif
		PyObject* SIMD = PROTECT(NEW_CHARACTER(ss.size()));
		nProtect ++;
		SET_ELEMENT(rv_ans, 1, SIMD);
		SET_STRING_ELT(nm, 1, mkChar("compiler.flag"));
		for (int i=0; i < (int)ss.size(); i++)
			SET_STRING_ELT(SIMD, i, mkChar(ss[i].c_str()));

		UNPROTECT(nProtect);

	COREARRAY_CATCH
}
*/

} // extern "C"
