# Makefile para compilar múltiplos programas com gfortran

# Compilador e flags
FC = gfortran
FFLAGS = -g -funderscoring -fno-second-underscore -Wall -O3 
# Diretórios e bibliotecas
LIB_DIR = /home/yan/Documents/PhD/Repo/build_tmdlib2211/lib
LIBS = -lTMDlib -lapfelxx -lgsl -lgslcblas -lm `lhapdf-config --cflags --libs` -lstdc++
#LHAPDFS_LIBs = `lhapdf-config --cflags --libs` -lstdc++

# Executáveis e seus arquivos fonte
EXEC1 = sigtot_kslinear
SRC1 = sigtot_pp_Z0_kslinear.f DADMUL.f dgauss.f vegas11.f dfint.f

EXEC2 = sigtot_ksnonlinear
SRC2 = sigtot_pp_Z0_ksnonlinear.f DADMUL.f dgauss.f vegas11.f dfint.f

EXEC3 = sigtot_gbw
SRC3 = sigtot_pp_Z0_gbw.f DADMUL.f dgauss.f vegas11.f dfint.f

EXEC4 = sigtot_ccfmset
SRC4 = sigtot_pp_Z0_ccfmset.f DADMUL.f dgauss.f vegas11.f dfint.f

# Regra padrão
all: $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4)

# Regras para compilar cada executável
$(EXEC1): $(SRC1)
	$(FC) $(FFLAGS) $(SRC1) -L$(LIB_DIR) $(LIBS) -o $(EXEC1)

$(EXEC2): $(SRC2)
	$(FC) $(FFLAGS) $(SRC2) -L$(LIB_DIR) $(LIBS) -o $(EXEC2)

$(EXEC3): $(SRC3)
	$(FC) $(FFLAGS) $(SRC3) -L$(LIB_DIR) $(LIBS) -o $(EXEC3)

$(EXEC4): $(SRC4)
	$(FC) $(FFLAGS) $(SRC4) -L$(LIB_DIR) $(LIBS) -o $(EXEC4)

# Limpeza dos arquivos gerados
clean:
	rm -f $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4)