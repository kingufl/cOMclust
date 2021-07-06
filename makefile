all: bin/cOMclust

bin/cOMclust:
	mkdir bin
	mkdir bin/results
	g++ src/*.cpp -o bin/cOMclust

clean: 
	rm -rf bin
