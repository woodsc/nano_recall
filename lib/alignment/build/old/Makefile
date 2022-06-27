#
#OR: mingw32-g++ alignment.cpp -shared -o alignment.so -Ic:/ruby/lib/ruby/1.8/i386-mswin32 -lmsvcrt-ruby18 -lmsvcrt -Lc:/ruby/bin -O3
#
#
clean:
	del alignment.so

alignment:
	mingw32-g++ alignment.cpp -shared -mwindows -static-libgcc -static-libstdc++ -o alignment.so -Ic:/ruby/lib/ruby/1.8/i386-mswin32 -lmsvcrt-ruby18 -lmsvcrt -Lc:/ruby/bin -O3 -Wl,--export-all-symbols 
