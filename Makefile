CC = clang
GHC = ghc
HSFLAGS = -O2 -threaded -outputdir obj
CFLAGS = -O2
SRCDIR = src
OBJDIR = obj
EXECDIR = bin
TESTDIR = test
C_SOURCES = $(wildcard $(SRCDIR)/*.c)
HS_SOURCES = $(wildcard $(SRCDIR)/*.hs)
C_OBJECTS = $(addprefix $(OBJDIR)/,$(notdir $(C_SOURCES:.c=.o)))
HS_OBJECTS = $(addprefix $(OBJDIR)/,$(notdir $(HS_SOURCES:.hs=.o)))
C_MAIN = $(TESTDIR)/main.c
HS_MAIN = $(TESTDIR)/Main.hs

C_EXECUTABLE = $(EXECDIR)/cmain
HS_EXECUTABLE = $(EXECDIR)/hsmain
STATIC_LIB = lib/libhs071.a
LIBS = -lipopt -llapack -ldl -lcoinhsl -lblas -lgfortran -lm -lquadmath
INCLUDES = -I$(SRCDIR)
LIB_DEPENDENCIES = $(STATIC_LIB) $(LIBS)

all : $(C_OBJECTS) $(STATIC_LIB) $(HS_OBJECTS) $(HS_EXECUTABLE) $(C_EXECUTABLE)

$(HS_OBJECTS) : $(HS_SOURCES)
	$(GHC) -c $(HSFLAGS)  $< -o $@

$(C_OBJECTS) : $(C_SOURCES)
	$(CC) -c $(CFLAGS)  $< -o $@

$(STATIC_LIB) : $(C_OBJECTS)
	ar rcs $(STATIC_LIB) $(C_OBJECTS)

$(HS_EXECUTABLE) : $(HS_MAIN) $(STATIC_LIB) $(HS_OBJECTS)
	$(GHC) $(HSFLAGS) $(INCLUDES) -o $(HS_EXECUTABLE) $(HS_MAIN) $(HS_SOURCES) $(LIB_DEPENDENCIES)

$(C_EXECUTABLE): $(C_MAIN) $(STATIC_LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $(C_MAIN) -o $(C_EXECUTABLE) $(LIB_DEPENDENCIES)


clean:
	rm -f $(C_EXECUTABLE) $(C_OBJECTS) $(HS_EXECUTABLE) $(HS_OBJECTS)
