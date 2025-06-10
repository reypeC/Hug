CC=g++
CFLAGS=-g

OBJ=random.o point.o sources.o dvector.o data.o hull.o hugmodel.o hugmh.o genericinteraction.o straussinteraction.o inconvexhullinteraction.o intruesourcesinteraction.o abcshadow_hug.o\

DIR= ./FGJK/
%.o: $(DIR)%.cpp
	$(CC) $(CFLAGS) -c $(DIR)$*.cpp

all : $(OBJ)  hugHardcore hug_abcHardcore
	$(CC) $(OBJ) hugHardcore.o -o ..\EXEC\hugHardcore.exe
	$(CC) $(OBJ) hug_abcHardcore.o -o ..\EXEC\hug_abcHardcore.exe

hugHardcore :
	$(CC) $(CFLAGS) -c hugHardcore.cpp

hug_abcHardcore :
	$(CC) $(CFLAGS) -c hug_abcHardcore.cpp

clean:
	del /Q *.o
	del /Q ..\EXEC\*.exe
