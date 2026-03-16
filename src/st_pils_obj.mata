*! version 1.0.0  06feb2026
*! DoNotDistribute


findfile "pils_macdefs.matah", path(".")
include `"`r(fn)'"'

version 16.0

mata:
mata set matastrict on
					//----------------------------//
					// create object
					//----------------------------//
void st_pils_create_object()
{
	pointer scalar	p

	if ((p = findexternal("`PILS_OBJ'")) == NULL) {
		p = crexternal("`PILS_OBJ'")
	}
	else {
		errprintf("`PILS_OBJ' already exists\n")
		exit(3001)
	}

	*p = `PILS'()
	
}
					//----------------------------//
					// clean object
					//----------------------------//
void st_pils_clean_object()
{
	rmexternal("`PILS_OBJ'")
}

					//----------------------------//
					// get object
					//----------------------------//
void st_pils_get_object()
{
	pointer p

	p = findexternal("`PILS_OBJ'")
	if (p == NULL) {
		errprintf("%s object not found\n", "`PILS_OBJ'")
		exit(3001)
	}

	st_global("s(pils_objname)", `"`PILS_OBJ'"')
}

end
