if(`test -n "-L/usr/local/Cellar/libxml2/2.9.9_2/lib -lxml2 -lz -lpthread -liconv -lm"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-L/usr/local/Cellar/libxml2/2.9.9_2/lib -lxml2 -lz -lpthread -liconv -lm
else
   setenv LD_LIBRARY_PATH -L/usr/local/Cellar/libxml2/2.9.9_2/lib -lxml2 -lz -lpthread -liconv -lm
endif

endif
