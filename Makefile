# Nom de l'exécutable
EXEC = sousmaillage

# Compilateur Fortran
FC = gfortran

# Nom de l'exécutable
EXEC = sousmaillage

# Compilateur et options
FC = gfortran
FLAGS = -g -O2 -Wall

# Fichiers sources et objets
SRC = main.f90 numeriques.f90 source.f90 lesfonction.f90
OBJ = $(SRC:.f90=.o)

# Règle par défaut
all: $(EXEC)

# Création de l'exécutable
$(EXEC): $(OBJ)
	$(FC) $(FLAGS) -o $@ $(OBJ)

# Compilation des fichiers sources
numeriques.o: numeriques.f90
	$(FC) $(FLAGS) -c numeriques.f90

source.o: source.f90 numeriques.o
	$(FC) $(FLAGS) -c source.f90
	
lesfonction.o: lesfonction.f90 numeriques.o source.o
	$(FC) $(FLAGS) -c lesfonction.f90

main.o: main.f90 numeriques.o lesfonction.o source.o
	$(FC) $(FLAGS) -c main.f90

# Nettoyage des fichiers générés
clean:
	rm -f $(OBJ) *.mod $(EXEC)

# Forcer la recompilation
rebuild: clean all
