#   File:       Absorb.make
#   Target:     Absorb
#   Sources:    Mortimer:MPWF:pcabs:absorb.f
#               Mortimer:MPWF:pcabs:abs.f
#               Mortimer:MPWF:pcabs:atomdata.f
#               Mortimer:MPWF:pcabs:cromer.f
#               Mortimer:MPWF:pcabs:crys_abs.f
#               Mortimer:MPWF:pcabs:element.f
#               Mortimer:MPWF:pcabs:gstring.f
#               Mortimer:MPWF:pcabs:matl_abs.f
#               Mortimer:MPWF:pcabs:raycomp.f
#               Mortimer:MPWF:pcabs:rdlat.f
#               Mortimer:MPWF:pcabs:sfcoef.f
#               Mortimer:MPWF:pcabs:str_fact.f
#               Mortimer:MPWF:pcabs:xnrg_b_d.f
#               Mortimer:MPWF:pcabs:xsec_b_d.f
#   Created:    Tuesday, November 22, 1994 11:30:39


OBJECTS = ¶
		Mortimer:MPWF:pcabs:absorb.f.o ¶
		Mortimer:MPWF:pcabs:abs.f.o ¶
		Mortimer:MPWF:pcabs:atomdata.f.o ¶
		Mortimer:MPWF:pcabs:cromer.f.o ¶
		Mortimer:MPWF:pcabs:crys_abs.f.o ¶
		Mortimer:MPWF:pcabs:element.f.o ¶
		Mortimer:MPWF:pcabs:gstring.f.o ¶
		Mortimer:MPWF:pcabs:matl_abs.f.o ¶
		Mortimer:MPWF:pcabs:raycomp.f.o ¶
		Mortimer:MPWF:pcabs:rdlat.f.o ¶
		Mortimer:MPWF:pcabs:sfcoef.f.o ¶
		Mortimer:MPWF:pcabs:str_fact.f.o ¶
		Mortimer:MPWF:pcabs:xnrg_b_d.f.o ¶
		Mortimer:MPWF:pcabs:xsec_b_d.f.o



Absorb ÄÄ Absorb.make {OBJECTS}
	Link -w -f -srt -ad 4 -t APPL -c SSRL ¶
		{OBJECTS} ¶
		"{Libraries}"Runtime.o ¶
		"{Libraries}"Interface.o ¶
		"{FLibraries}"FORTRANlib.o ¶
		"{FLibraries}"IntrinsicLib.o ¶
		"{FLibraries}"FSANELib.o ¶
		"{FLibraries}"FABSORPTIONLib.o ¶
		-o Absorb
	Echo "Include ¶"{FLibraries}Fresources.r¶";" > "{FLibraries}Resource.inc"
	Rez "{FLibraries}Resource.inc" -a -m -o "Absorb"
	 FSIZE "Absorb"
Mortimer:MPWF:pcabs:absorb.f.o Ä Absorb.make Mortimer:MPWF:pcabs:absorb.f
	 FORTRAN Mortimer:MPWF:pcabs:absorb.f -opt=1 
Mortimer:MPWF:pcabs:abs.f.o Ä Absorb.make Mortimer:MPWF:pcabs:abs.f
	 FORTRAN Mortimer:MPWF:pcabs:abs.f -opt=1 
Mortimer:MPWF:pcabs:atomdata.f.o Ä Absorb.make Mortimer:MPWF:pcabs:atomdata.f
	 FORTRAN Mortimer:MPWF:pcabs:atomdata.f -opt=1 
Mortimer:MPWF:pcabs:cromer.f.o Ä Absorb.make Mortimer:MPWF:pcabs:cromer.f
	 FORTRAN Mortimer:MPWF:pcabs:cromer.f -opt=1 
Mortimer:MPWF:pcabs:crys_abs.f.o Ä Absorb.make Mortimer:MPWF:pcabs:crys_abs.f
	 FORTRAN Mortimer:MPWF:pcabs:crys_abs.f -opt=1 
Mortimer:MPWF:pcabs:element.f.o Ä Absorb.make Mortimer:MPWF:pcabs:element.f
	 FORTRAN Mortimer:MPWF:pcabs:element.f -opt=1 
Mortimer:MPWF:pcabs:gstring.f.o Ä Absorb.make Mortimer:MPWF:pcabs:gstring.f
	 FORTRAN Mortimer:MPWF:pcabs:gstring.f -opt=1 
Mortimer:MPWF:pcabs:matl_abs.f.o Ä Absorb.make Mortimer:MPWF:pcabs:matl_abs.f
	 FORTRAN Mortimer:MPWF:pcabs:matl_abs.f -opt=1 
Mortimer:MPWF:pcabs:raycomp.f.o Ä Absorb.make Mortimer:MPWF:pcabs:raycomp.f
	 FORTRAN Mortimer:MPWF:pcabs:raycomp.f -opt=1 
Mortimer:MPWF:pcabs:rdlat.f.o Ä Absorb.make Mortimer:MPWF:pcabs:rdlat.f
	 FORTRAN Mortimer:MPWF:pcabs:rdlat.f -opt=1 
Mortimer:MPWF:pcabs:sfcoef.f.o Ä Absorb.make Mortimer:MPWF:pcabs:sfcoef.f
	 FORTRAN Mortimer:MPWF:pcabs:sfcoef.f -opt=1 
Mortimer:MPWF:pcabs:str_fact.f.o Ä Absorb.make Mortimer:MPWF:pcabs:str_fact.f
	 FORTRAN Mortimer:MPWF:pcabs:str_fact.f -opt=1 
Mortimer:MPWF:pcabs:xnrg_b_d.f.o Ä Absorb.make Mortimer:MPWF:pcabs:xnrg_b_d.f
	 FORTRAN Mortimer:MPWF:pcabs:xnrg_b_d.f -opt=1 
Mortimer:MPWF:pcabs:xsec_b_d.f.o Ä Absorb.make Mortimer:MPWF:pcabs:xsec_b_d.f
	 FORTRAN Mortimer:MPWF:pcabs:xsec_b_d.f -opt=1 
