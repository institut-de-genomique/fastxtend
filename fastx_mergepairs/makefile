OBJLIBS = ../libfastx/libfastx.a
 
all : fastx_mergepairs

fastx_mergepairs : fastx_mergepairs.cpp .deps/fastx_mergepairs.Po $(OBJLIBS)
	if g++ -DHAVE_CONFIG_H -I. -I. -I../..   -I../libfastx   -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1 -MT fastx_mergepairs.o -MD -MP -MF ".deps/fastx_mergepairs.Tpo" -c -o fastx_mergepairs.o fastx_mergepairs.cpp; \
	then mv -f ".deps/fastx_mergepairs.Tpo" ".deps/fastx_mergepairs.Po"; else rm -f ".deps/fastx_mergepairs.Tpo"; exit 1; fi; \
	g++  -g -O2 -Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -DDEBUG -g -O1 -DDEBUG -g -O1   -o fastx_mergepairs  fastx_mergepairs.o ../libfastx/libfastx.a

install :
	cp -f fastx_mergepairs ../bin

clean :
	rm -f fastx_mergepairs fastx_mergepairs.o
