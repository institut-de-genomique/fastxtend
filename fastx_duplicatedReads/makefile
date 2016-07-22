OBJLIBS = ../libfastx/libfastx.a
 
all : fastx_duplicatedReads

fastx_duplicatedReads : fastx_duplicatedReads.cpp .deps/fastx_duplicatedReads.Po $(OBJLIBS)
	if g++ -DHAVE_CONFIG_H -I. -I. -I../..   -I../libfastx   -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1 -MT fastx_duplicatedReads.o -MD -MP -MF ".deps/fastx_duplicatedReads.Tpo" -c -o fastx_duplicatedReads.o fastx_duplicatedReads.cpp; \
	then mv -f ".deps/fastx_duplicatedReads.Tpo" ".deps/fastx_duplicatedReads.Po"; else rm -f ".deps/fastx_duplicatedReads.Tpo"; exit 1; fi; \
	g++  -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1   -o fastx_duplicatedReads  fastx_duplicatedReads.o ../libfastx/libfastx.a

install :
	cp -f fastx_duplicatedReads ../bin

clean :
	rm -f fastx_duplicatedReads fastx_duplicatedReads.o
