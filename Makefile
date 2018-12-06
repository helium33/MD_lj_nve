objects_initStruct = init_struct.o tools.o 
objects_md = md.o tools.o
objects_analysis = traj_analy.o tools.o

all: md initStruct trajAnaly

md: $(objects_md)
	gcc -o md $(objects_md) -lm
initStruct: $(objects_initStruct)
	gcc -o initStruct $(objects_initStruct) -lm
trajAnaly: $(objects_analysis)
	gcc -o trajAnaly $(objects_analysis) -lm

tools.o: tools.c tools.h
	gcc -c tools.c  -lm
init_struct.o: init_struct.c tools.h
	gcc -c init_struct.c -lm
md.o: md.c tools.h
	gcc -c md.c -lm
traj_analy.o: traj_analy.c tools.h
	gcc -c traj_analy.c -lm
clean:
	rm *.o
