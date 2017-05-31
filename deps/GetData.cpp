// ===========================================================
//
// GetData.cpp: Get data from the GDS file
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
#include "ReadByVariant.h"


using namespace JSeqArray;

extern "C"
{

// ===========================================================
// Get data from a working space
// ===========================================================

static bool is_logical(PdGDSObj Node)
{
	char classname[32];
	classname[0] = 0;
	GDS_Node_GetClassName(Node, classname, sizeof(classname));
	return (strcmp(classname, "dBit1") == 0);
}


// get data
static jl_array_t* VarGetData(CFileInfo &File, const char *name)
{
	static const char *ERR_DIM = "Invalid dimension of '%s'.";

	jl_array_t *rv_ans = NULL;
	TSelection &Sel = File.Selection();

	if (strcmp(name, "sample.id") == 0)
	{
		// ===========================================================
		// sample.id

		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.SampleNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss = Sel.pSample();
		rv_ans = GDS_JArray_Read(N, NULL, NULL, &ss, svCustom);

	} else if (strcmp(name, "position") == 0)
	{
		int n = File.VariantSelNum();
		jl_value_t *atype = jl_apply_array_type(jl_int32_type, 1);
		rv_ans = jl_alloc_array_1d(atype, n);
		if (n > 0)
		{
			const int *base = &File.Position()[0];
			int *p = (int*)jl_array_data(rv_ans);
			C_BOOL *s = Sel.pVariant();
			for (size_t m=File.VariantNum(); m > 0; m--)
			{
				if (*s++) *p++ = *base;
				base ++;
			}
		}

	} else if (strcmp(name, "chromosome") == 0)
	{
		int n = File.VariantSelNum();
		jl_value_t *atype = jl_apply_array_type(jl_string_type, 1);
		rv_ans = jl_alloc_array_1d(atype, n);
		JL_GC_PUSH1(&rv_ans);
		if (n > 0)
		{
			CChromIndex &Chrom = File.Chromosome();
			jl_value_t **p = (jl_value_t**)jl_array_data(rv_ans);
			C_BOOL *s = Sel.pVariant();
			size_t m = File.VariantNum();
			string lastss;
			jl_value_t *last = NULL;
			for (size_t i=0; i < m; i++)
			{
				if (*s++)
				{
					const string &ss = Chrom[i];
					if (ss != lastss)
					{
						lastss = ss;
						last = NULL;
					}
					if (!last)
						last = jl_pchar_to_string(lastss.c_str(), lastss.size());
					*p++ = last;
					jl_gc_wb(rv_ans, last);
				}
			}
		}
		JL_GC_POP();
	
	} else if ( (strcmp(name, "variant.id")==0) ||
		(strcmp(name, "allele")==0) ||
		(strcmp(name, "annotation/id")==0) ||
		(strcmp(name, "annotation/qual")==0) ||
		(strcmp(name, "annotation/filter")==0) )
	{
		// ===========================================================
		// variant.id, allele, annotation/id, annotation/qual, annotation/filter

		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss = Sel.pVariant();
		rv_ans = GDS_JArray_Read(N, NULL, NULL, &ss, svCustom);

	} else if (strcmp(name, "genotype") == 0)
	{
		// ===========================================================
		// genotypic data

		int nSample  = File.SampleSelNum();
		int nVariant = File.VariantSelNum();

		if ((nSample > 0) && (nVariant > 0))
		{
			// initialize GDS genotype Node
			CApply_Variant_Geno NodeVar(File);
			// set
			jl_value_t *atype = jl_apply_array_type(jl_uint8_type, 3);
			rv_ans = jl_alloc_array_3d(atype, File.Ploidy(), nSample, nVariant);
			C_UInt8 *base = (C_UInt8*)jl_array_data(rv_ans);
			ssize_t SIZE = (ssize_t)nSample * File.Ploidy();
			do {
				NodeVar.ReadGenoData(base);
				base += SIZE;
			} while (NodeVar.Next());
		} else {
			jl_value_t *atype = jl_apply_array_type(jl_uint8_type, 3);
			rv_ans = jl_alloc_array_3d(atype, File.Ploidy(), nSample, 0);
		}

	} else if (strcmp(name, "@genotype") == 0)
	{
		static const char *VarName = "genotype/@data";
		PdAbstractArray N = File.GetObj(VarName, TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, VarName);
		// read
		C_BOOL *ss = Sel.pVariant();
		rv_ans = GDS_JArray_Read(N, NULL, NULL, &ss, svInt32);

	} else if (strcmp(name, "#dosage")==0 || strcmp(name, "$dosage")==0)
	{
		// ===========================================================
		// dosage data

		ssize_t nSample  = File.SampleSelNum();
		ssize_t nVariant = File.VariantSelNum();

		if ((nSample > 0) && (nVariant > 0))
		{
			// initialize GDS genotype Node
			CApply_Variant_Dosage NodeVar(File);
			// set
			jl_value_t *atype = jl_apply_array_type(jl_uint8_type, 2);
			rv_ans = jl_alloc_array_2d(atype, nSample, nVariant);
			C_UInt8 *base = (C_UInt8*)jl_array_data(rv_ans);
			do {
				NodeVar.ReadDosage(base);
				base += nSample;
			} while (NodeVar.Next());
		} else {
			jl_value_t *atype = jl_apply_array_type(jl_uint8_type, 2);
			rv_ans = jl_alloc_array_2d(atype, nSample, 0);
		}

	} else if (strcmp(name, "phase") == 0)
	{
		// ===========================================================
		// phase/

		PdAbstractArray N = File.GetObj("phase/data", TRUE);
		// check
		int ndim = GDS_Array_DimCnt(N);
		C_Int32 dim[4];
		GDS_Array_GetDim(N, dim, 3);
		if (ndim<2 || ndim>3 || dim[0]!= File.VariantNum() ||
				dim[1]!=File.SampleNum())
			throw ErrSeqArray(ERR_DIM, name);
		// read
		C_BOOL *ss[3] = { Sel.pVariant(), Sel.pSample(), NULL };
		if (ndim == 3)
			ss[2] = NeedArrayTRUEs(dim[2]);
		rv_ans = GDS_JArray_Read(N, NULL, NULL, ss, svCustom);

	} else if (strncmp(name, "annotation/info/@", 17) == 0)
	{
		if (File.GetObj(name, FALSE) != NULL)
		{
			CIndex &V = File.VarIndex(name);
			rv_ans = V.GetLen_Sel(Sel.pVariant());
		}

	} else if (strncmp(name, "annotation/info/", 16) == 0)
	{
		// ===========================================================
		// annotation/info

		GDS_PATH_PREFIX_CHECK(name);
		PdAbstractArray N = File.GetObj(name, TRUE);
		int ndim = GDS_Array_DimCnt(N);
		if ((ndim!=1) && (ndim!=2))
			throw ErrSeqArray(ERR_DIM, name);

		string name2 = GDS_PATH_PREFIX(name, '@');
		PdAbstractArray N_idx = File.GetObj(name2.c_str(), FALSE);
		if (N_idx == NULL)
		{
			// no index
			C_Int32 dim[4];
			GDS_Array_GetDim(N, dim, 2);
			C_BOOL *ss[2] = { Sel.pVariant(), NULL };
			if (ndim == 2)
				ss[1] = NeedArrayTRUEs(dim[1]);
			C_SVType SV = svCustom;  // is_logical
			rv_ans = GDS_JArray_Read(N, NULL, NULL, ss, SV);

		} else {
			// with index
			CIndex &V = File.VarIndex(name2);
			int var_start, var_count;
			vector<C_BOOL> var_sel;

			jl_array_t *Index = NULL;
			jl_array_t *Dat = NULL;
			JL_GC_PUSH2(&Index, &Dat);
			Index = V.GetLen_Sel(Sel.pVariant(), var_start, var_count, var_sel);

			C_BOOL *ss[2] = { &var_sel[0], NULL };
			C_Int32 dimst[2]  = { var_start, 0 };
			C_Int32 dimcnt[2] = { var_count, 0 };
			if (ndim == 2)
			{
				GDS_Array_GetDim(N, dimcnt, 2);
				dimcnt[0] = var_count;
			}
			Dat = GDS_JArray_Read(N, dimst, dimcnt, ss, svCustom);

			jl_value_t *atype = jl_apply_array_type(jl_any_type, 1);
			rv_ans = jl_alloc_array_1d(atype, 2);
			void **ptr = (void**)jl_array_data(rv_ans);
			ptr[0] = Index; jl_gc_wb(rv_ans, Index);
			ptr[1] = Dat; jl_gc_wb(rv_ans, Dat);
			JL_GC_POP();
		}

	} else if (strncmp(name, "annotation/format/@", 19) == 0)
	{
		string name2(name);
		name2.erase(18, 1).append("/@data");
		if (File.GetObj(name2.c_str(), FALSE) != NULL)
		{
			CIndex &V = File.VarIndex(name2.c_str());
			rv_ans = V.GetLen_Sel(Sel.pVariant());
		}

	} else if (strncmp(name, "annotation/format/", 18) == 0)
	{
		// ===========================================================
		// annotation/format

		GDS_PATH_PREFIX_CHECK(name);
		string name1 = string(name) + "/data";
		string name2 = string(name) + "/@data";
		PdAbstractArray N = File.GetObj(name1.c_str(), TRUE);

		// with index
		CIndex &V = File.VarIndex(name2);
		int var_start, var_count;
		vector<C_BOOL> var_sel;

		jl_array_t *Index = NULL;
		jl_array_t *Dat = NULL;
		JL_GC_PUSH2(&Index, &Dat);
		Index = V.GetLen_Sel(Sel.pVariant(), var_start, var_count, var_sel);

		C_BOOL *ss[2] = { &var_sel[0], Sel.pSample() };
		C_Int32 dimst[2]  = { var_start, 0 };
		C_Int32 dimcnt[2];
		GDS_Array_GetDim(N, dimcnt, 2);
		dimcnt[0] = var_count;
		Dat = GDS_JArray_Read(N, dimst, dimcnt, ss, svCustom);

		jl_value_t *atype = jl_apply_array_type(jl_any_type, 1);
		rv_ans = jl_alloc_array_1d(atype, 2);
		void **ptr = (void**)jl_array_data(rv_ans);
		ptr[0] = Index; jl_gc_wb(rv_ans, Index);
		ptr[1] = Dat; jl_gc_wb(rv_ans, Dat);
		JL_GC_POP();

	} else if (strncmp(name, "sample.annotation/", 18) == 0)
	{
		// ===========================================================
		// sample.annotation

		GDS_PATH_PREFIX_CHECK(name);
		PdAbstractArray N = File.GetObj(name, TRUE);
		// check
		int ndim = GDS_Array_DimCnt(N);
		if ((ndim!=1) && (ndim!=2))
			throw ErrSeqArray(ERR_DIM, name);
		C_Int32 dim[2];
		GDS_Array_GetDim(N, dim, 2);
		if (dim[0] != File.SampleNum())
			throw ErrSeqArray(ERR_DIM, name);

		C_BOOL *ss[2] = { Sel.pSample(), NULL };
		if (ndim == 2)
			ss[1] = NeedArrayTRUEs(dim[1]);
		rv_ans = GDS_JArray_Read(N, NULL, NULL, ss, svCustom);

	} else if (strcmp(name, "#chrom_pos")==0 || strcmp(name, "$chrom_pos")==0)
	{
		// ===========================================================
		// chromosome-position

		PdAbstractArray N1 = File.GetObj("chromosome", TRUE);
		PdAbstractArray N2 = File.GetObj("position", TRUE);
		C_Int64 n1 = GDS_Array_GetTotalCount(N1);
		C_Int64 n2 = GDS_Array_GetTotalCount(N2);
		if ((n1 != n2) || (n1 != File.VariantNum()))
			throw ErrSeqArray("Invalid dimension of 'chromosome' and 'position'.");

		vector<string> chr;
		vector<C_Int32> pos;

		int n = File.VariantSelNum();
		chr.resize(n);
		pos.resize(n);
		C_BOOL *ss = Sel.pVariant();

		GDS_Array_ReadDataEx(N1, NULL, NULL, &ss, &chr[0], svStrUTF8);
		GDS_Array_ReadDataEx(N2, NULL, NULL, &ss, &pos[0], svInt32);

		char buf1[1024] = { 0 };
		char buf2[1024] = { 0 };
		char *p1 = buf1, *p2 = buf2;
		int dup = 0;

		jl_value_t *atype = jl_apply_array_type(jl_string_type, 1);
		rv_ans = jl_alloc_array_1d(atype, n1);
		JL_GC_PUSH1(&rv_ans);
		jl_value_t **p = (jl_value_t**)jl_array_data(rv_ans);

		for (size_t i=0; i < (size_t)n1; i++,p++)
		{
			snprintf(p1, sizeof(buf1), "%s_%d", chr[i].c_str(), pos[i]);
			if (strcmp(p1, p2) == 0)
			{
				dup ++;
				snprintf(p1, sizeof(buf1), "%s_%d_%d", chr[i].c_str(),
					pos[i], dup);
				*p = jl_cstr_to_string(p1);
				jl_gc_wb(rv_ans, *p);
			} else {
				char *tmp;
				tmp = p1; p1 = p2; p2 = tmp;
				*p = jl_cstr_to_string(p2);
				jl_gc_wb(rv_ans, *p);
				dup = 0;
			}
		}

		JL_GC_POP();

	} else if (strcmp(name, "#num_allele")==0 || strcmp(name, "$num_allele")==0)
	{
		// ===========================================================
		// the number of distinct alleles

		ssize_t nVariant = File.VariantSelNum();
		jl_value_t *atype = jl_apply_array_type(jl_int32_type, 1);
		rv_ans = jl_alloc_array_1d(atype, nVariant);
		int *p = (int*)jl_array_data(rv_ans);

		CApply_Variant_NumAllele NodeVar(File);
		for (ssize_t i=0; i < nVariant; i++)
		{
			p[i] = NodeVar.GetNumAllele();
			NodeVar.Next();
		}

/*
	} else if (strcmp(name, "#ref")==0 || strcmp(name, "$ref")==0)
	{
		// ===========================================================
		// the reference allele

		PdAbstractArray N = File.GetObj("allele", TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		size_t n = File.VariantSelNum();
		vector<string> buffer(n);
		C_BOOL *ss = Sel.pVariant();
		GDS_Array_ReadDataEx(N, NULL, NULL, &ss, &buffer[0], svStrUTF8);
		// output
		rv_ans = numpy_new_string(n);
		PyObject **pi = (PyObject**)numpy_getptr(rv_ans);
		for (size_t i=0; i < n; i++)
		{
			const char *p = buffer[i].c_str();
			size_t m = 0;
			for (const char *s=p; *s!=',' && *s!=0; s++) m++;
			numpy_setval(rv_ans, pi, PYSTR_SET2(p, m));
			pi ++;
		}

	} else if (strcmp(name, "#alt")==0 || strcmp(name, "$alt")==0)
	{
		// ===========================================================
		// the reference allele

		PdAbstractArray N = File.GetObj("allele", TRUE);
		// check
		if ((GDS_Array_DimCnt(N) != 1) ||
				(GDS_Array_GetTotalCount(N) != File.VariantNum()))
			throw ErrSeqArray(ERR_DIM, name);
		// read
		size_t n = File.VariantSelNum();
		vector<string> buffer(n);
		C_BOOL *ss = Sel.pVariant();
		GDS_Array_ReadDataEx(N, NULL, NULL, &ss, &buffer[0], svStrUTF8);
		// output
		rv_ans = numpy_new_string(n);
		PyObject **pi = (PyObject**)numpy_getptr(rv_ans);
		for (size_t i=0; i < n; i++)
		{
			const char *p = buffer[i].c_str();
			for (; *p!=',' && *p!=0; p++);
			if (*p == ',') p++;
			numpy_setval(rv_ans, pi, PYSTR_SET(p));
			pi ++;
		}
*/
	} else {
		throw ErrSeqArray(
			"'%s' is not a standard variable name, and the standard format:\n"
			"    sample.id, variant.id, position, chromosome, allele, genotype\n"
			"    annotation/id, annotation/qual, annotation/filter\n"
			"    annotation/info/VARIABLE_NAME, annotation/format/VARIABLE_NAME\n"
			"    sample.annotation/VARIABLE_NAME", name);
	}

	return rv_ans;
}


/// Get data from a working space
COREARRAY_DLL_EXPORT jl_array_t* SEQ_GetData(int file_id, const char *name)
{
	jl_array_t *rv_ans = NULL;
	COREARRAY_TRY
		// File information
		CFileInfo &File = GetFileInfo(file_id);
		// Get data
		rv_ans = VarGetData(File, name);
	COREARRAY_CATCH
	return rv_ans;
}


/// Apply functions over variants in block
COREARRAY_DLL_EXPORT jl_array_t* SEQ_BApply_Variant(int file_id,
	jl_value_t *name, jl_function_t *fun, const char *asis,
	int bsize, C_BOOL verbose, jl_array_t *args)
{
	if (bsize < 1)
		jl_error("'bsize' must be >= 1.");

	jl_array_t *rv_ans = NULL;
	COREARRAY_TRY

		// get a list of variable name
		vector<string> name_list;
		if (jl_is_string(name))
		{
			// String
			name_list.push_back(jl_string_ptr(name));
		} else {
			// Vector{String}
			size_t n = jl_array_len(name);
			jl_value_t **p= (jl_value_t**)jl_array_data(name);
			name_list.resize(n);
			for (size_t i=0; i < n; i++)
				name_list[i] = jl_string_ptr(*p++);
		}
		if (name_list.empty())
			throw ErrSeqArray("'name' should be specified.");

		// File information
		CFileInfo &File = GetFileInfo(file_id);
		// Selection
		TSelection &Selection = File.Selection();

		// the number of selected variants
		int nVariant = File.VariantSelNum();
		if (nVariant <= 0)
			throw ErrSeqArray("There is no selected variant.");

		// the number of data blocks
		int NumBlock = nVariant / bsize;
		if (nVariant % bsize) NumBlock ++;

		// asis
		if (strcmp(asis, "list")==0 || strcmp(asis, "unlist")==0)
		{
			jl_value_t *atype = jl_apply_array_type(jl_any_type, 1);
			rv_ans = jl_alloc_array_1d(atype, NumBlock);
		} else if (strcmp(asis, "none") != 0)
		{
			throw ErrSeqArray("'asis' should be 'none', 'list' or 'unlist'.");
		}

		// the number of variables
		size_t nVar = name_list.size();
		// the number of additional parameters
		size_t nArgs = jl_array_len(args);
		jl_value_t **ArgPtr = (jl_value_t**)jl_array_data(args);

		// protect
		JL_GC_PUSH1(&rv_ans);

		// local selection
		File.SelList.push_back(TSelection());
		TSelection &Sel = File.SelList.back();
		Sel.Sample = Selection.Sample;
		Sel.Variant.resize(File.VariantNum());

		C_BOOL *pBase, *pSel, *pEnd;
		pBase = pSel = Selection.pVariant();
		pEnd = pBase + Selection.Variant.size();

		// progress object
		CProgressStdOut progress(NumBlock, verbose!=0);

		// for-loop
		for (int idx=0; idx < NumBlock; idx++)
		{
			// assign sub-selection
			{
				C_BOOL *pNewSel = Sel.pVariant();
				memset(pNewSel, 0, Sel.Variant.size());
				// for-loop
				for (int bs=bsize; bs > 0; bs--)
				{
					while ((pSel < pEnd) && (*pSel == FALSE))
						pSel ++;
					if (pSel < pEnd)
					{
						pNewSel[pSel - pBase] = TRUE;
						pSel ++;
					} else
						break;
				}
			}

			jl_value_t **list_args;
			JL_GC_PUSHARGS(list_args, nVar+nArgs);

			// load data
			for (size_t i=0; i < nVar; i++)
				list_args[i] = (jl_value_t*)VarGetData(File, name_list[i].c_str());
			for (size_t i=0; i < nArgs; i++)
				list_args[nVar+i] = ArgPtr[i];

			// call Julia function
			jl_value_t *rv = jl_call(fun, list_args, nVar+nArgs);
			// save
			if (rv_ans)
			{
				void **ptr = (void**)jl_array_data(rv_ans);
				ptr[idx] = rv;
				jl_gc_wb(rv_ans, rv);
			}

			JL_GC_POP();

			progress.Forward();
		}

		File.SelList.pop_back();
		JL_GC_POP();

	COREARRAY_CATCH
	if (!rv_ans) rv_ans = (jl_array_t*)jl_nothing;
	return rv_ans;
}

} // extern "C"
