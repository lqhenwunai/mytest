CXX=g++
Objects=main.o\
	FileRead.o\
	Matrix-tools.o\
	Matrix-container.o\
	Matrix-math.o\
	Util.o\
	Message.o

mycode:$(Objects)
		$(CXX) -DUSE_DEBUG -o mycode $(Objects)
main.o:main.cpp
		$(CXX) -DUSE_DEBUG -c main.cpp
FileRead:FileRead.cpp FileRead.h
		$(CXX) -DUSE_DEBUG -c FileRead.cpp
Matrix-tools:Matrix-tools.cpp Matrix-container.cpp Matrix-tools.h
		$(CXX) -DUSE_DEBUG -c Matrix-tools.cpp
		$(CXX) -DUSE_DEBUG -c Matrix-container.cpp
Matrix-math:Matrix-math.cpp Matrix-math.h
		$(CXX) -DUSE_DEBUG -c Matrix-math.cpp
Util.o:Util.cpp Util.h
		$(CXX) -DUSE_DEBUG -c Util.cpp
Message.o:Message.cpp Message.h
		$(CXX) -DUSE_DEBUG -c Message.cpp

