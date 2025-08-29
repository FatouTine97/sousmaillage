# Nom de l'exécutable
EXEC = sousmaillage

# Compilateur Fortran
FC = gfortran

# Nom de l'exécutable
EXEC = sousmaillage

# Compilateur et options
FC = gfortran
FLAGS = -O3 -mavx -fbacktrace 

#-Wall

# Fichiers sources et objets
SRC = main.f90 numeriques.f90 source.f90 lesfonctioni.f90  reference.f90# pml.f90 #dispersion.f90
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
 #pml.o : pml.f90 numeriques.o	$(FC) $(FLAGS) -c pml.f90

lesfonctioni.o: lesfonction.f90 numeriques.o source.o #pml.o
	$(FC) $(FLAGS) -c lesfonctioni.f90

reference.o: reference.f90 lesfonctioni.o numeriques.o source.o #pml.o
	$(FC) $(FLAGS) -c reference.f90

#dispersion.o: dispersion.f90 numeriques.o lesfonctioni.o source.o reference.o$(FC) $(FLAGS) -c dispersion.f90

main.o: main.f90 numeriques.o lesfonctioni.o source.o reference.o pml.o #dispersion.o
	$(FC) $(FLAGS) -c main.f90

# Nettoyage des fichiers générés
clean:
	rm -f $(OBJ) *.mod $(EXEC)

# Forcer la recompilation
rebuild: clean all
