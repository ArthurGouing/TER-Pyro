# Compilateur utilisé
CC=g++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -I ../Eigen/Eigen -std=c++11
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -w -I ../Eigen/Eigen -std=c++11

# On choisit comment on compile
CXX_FLAGS = $(OPTIM_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = main.cc MeshAdapt.cpp TimeScheme.cpp FiniteVolume.cpp DataFile.cpp Function.cpp Solution.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
