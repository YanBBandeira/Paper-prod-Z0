# Makefile para compilar múltiplos programas com gfortran

# Compilador e flags
FC = gfortran
FFLAGS = -g -funderscoring -fno-second-underscore -Wall -O6 -ffast-math

# Diretórios e bibliotecas
LIB_DIR = /home/yan/Documents/PhD/Repo/build_tmdlib2211/lib
LIBS = -lTMDlib -lapfelxx -lgsl -lgslcblas -lm `lhapdf-config --cflags --libs` -lstdc++
#LHAPDFS_LIBs = `lhapdf-config --cflags --libs` -lstdc++

# Executáveis e seus arquivos fonte
EXEC1 = kslinear
SRC1 = sigtot_pp_Z0.f DADMUL.f dgauss.f vegas11.f

#EXEC2 = ksnonlinear
#SRC2 = Z0-ksnonlinear.f DADMUL.f dgauss.f CT14Pdf.f

#EXEC3 = gbw
#SRC3 = Z0-gbw.f DADMUL.f dgauss.f CT14Pdf.f

#EXEC4 = ccfmset
#SRC4 = Z0-ccfmset.f DADMUL.f dgauss.f CT14Pdf.f

# Regra padrão
all: $(EXEC1) #$(EXEC2) $(EXEC3) $(EXEC4)

# Regras para compilar cada executável
$(EXEC1): $(SRC1)
	$(FC) $(FFLAGS) $(SRC1) -L$(LIB_DIR) $(LIBS) -o $(EXEC1)

#$(EXEC2): $(SRC2)
#	$(FC) $(FFLAGS) $(SRC2) -L$(LIB_DIR) $(LIBS) -o $(EXEC2)

#$(EXEC3): $(SRC3)
#	$(FC) $(FFLAGS) $(SRC3) -L$(LIB_DIR) $(LIBS) -o $(EXEC3)

#$(EXEC4): $(SRC4)
#	$(FC) $(FFLAGS) $(SRC4) -L$(LIB_DIR) $(LIBS) -o $(EXEC4)

# Limpeza dos arquivos gerados
clean:
	rm -f $(EXEC1) 
# $(EXEC2) $(EXEC3) $(EXEC4)