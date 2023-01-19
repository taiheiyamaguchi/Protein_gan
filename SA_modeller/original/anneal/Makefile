OBJS=dynamics_param.o main.o molecule.o potential.o
TARGET=anneal

$(TARGET): $(OBJS)
	g++ -o $@ $(OBJS)

.C.o:
	g++ -c $<

clean:
	rm -f $(TARGET) $(OBJS)
