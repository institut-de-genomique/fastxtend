OBJLIBS = ../libfastx/libfastx.a
 
all : fastx_clean

fastx_clean : fastx_clean.cpp .deps/fastx_clean.Po $(OBJLIBS)
	if g++ -DHAVE_CONFIG_H -I. -I. -I../..   -I../libfastx   -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1 -MT fastx_clean.o -MD -MP -MF ".deps/fastx_clean.Tpo" -c -o fastx_clean.o fastx_clean.cpp; \
	then mv -f ".deps/fastx_clean.Tpo" ".deps/fastx_clean.Po"; else rm -f ".deps/fastx_clean.Tpo"; exit 1; fi; \
	g++  -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1   -o fastx_clean  fastx_clean.o ../libfastx/libfastx.a

install :
	cp -f fastx_clean ../bin

clean :
	rm -f fastx_clean fastx_clean.o
