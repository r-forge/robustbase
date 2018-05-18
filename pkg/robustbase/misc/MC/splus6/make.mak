# Rules for Watcom's WMAKE

!include chapter.mif

.EXTENSIONS
.EXTENSIONS: .dll .lib .res .def .obj .rc .c .cxx .cpp .f .for 

.BEFORE
	set include=$(%WATCOM)\H;$(%WATCOM)\MFC\INCLUDE;$(%WATCOM)\H\NT;$(SHOME)\include
#	set lib=

SPEXPORT = $(SHOME)\cmd\spexport.exe

LINK=wlink
CC=wcc386
CXX=wpp386
FC=wfc386
RC=wrc
LFLAGS = SYSTEM nt_dll &
			LIBPATH $(SHOME)\lib\watcom &
			LIBRARY sqpe.lib,kernel32.lib,user32.lib &
			OPTION DESCRIPTION 'S-PLUS Chapter DLL' &
			OPTION CASEEXACT &
			EXPORT @S.def

CFLAGS = -bt=nt -sg -bd -d0 -zw -5s -dWIN32=1 -zq -fp5
CXXFLAGS = $(CFLAGS)
FFLAGS = -nowarn -on -bd -sc -sg -quiet -stack -5 -fp5 
RFLAGS = /bt=nt /r /fo S.res /DLIBNAME=\"S.dll\"

.rc.res :
	$(RC) $(RFLAGS) $<

.c.obj :
	$(CC) $(CFLAGS) $<

.cxx.obj :
	$(CXX) $(CXXFLAGS) $<

.cpp.obj :
	$(CXX) $(CXXFLAGS) $<

.f.obj :
	$(FC) $(FFLAGS) $<

.for.obj :
	$(FC) $(FFLAGS) $<

S.dll: $(RES) $(OBJ) S.def
	$(LINK) $(LFLAGS) NAME $@ FILE {$(OBJ)}

S.def: $(OBJ)
	$(SPEXPORT) -w -o S.def $(OBJ)

#boot:
#	@if test -s all.Sdata; 	  then (BOOTING_S="TRUE" export BOOTING_S;  echo "terminate('should have been booting S')"| $(SHOME)/S); 	fi


#clean:
#	-del $(OBJ)


