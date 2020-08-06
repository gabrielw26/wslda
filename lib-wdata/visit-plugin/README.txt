Create static library libntg_io.a 

    Go to bsk folder and compile ntg_io.c file using command:
        g++ -O3 -c wdata.c -fPIC
    and next create static library (copy it to folder with visit_plugin)
        ar crf visit-plugin/libwdata.a wdata.o
    
To generate `Makefiles`:
   
    xml2cmake -clobber wdata.xml
    rm CMakeCache.txt
    cmake -DCMAKE_BUILD_TYPE:STRING=Debug
    
To compile and install plugin:
    
    make
    
Making plugin to be public (availabe for everyone) use xml2cmake with additional -public option:

    xml2cmake -public -clobber wdata.xml 
    
    
--------------------------------------------------------------------    
INSTALLING ON DWARF:
--------------------------------------------------------------------
* use node2066
* module load visit/2_13_1
* execute makelinks.sh
* from node2072 execute: install2072.sh
* from node2072 execute: makelinks.sh
    
--------------------------------------------------------------------    
PROBLEMS YOU MAY ENCOUNTER:
--------------------------------------------------------------------

1. If you find problem of type:

Scanning dependencies of target EwdataDatabase_par
...
make[2]: *** No rule to make target '/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/net/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/wci-cl1.llnl.gov/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/vol/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/home/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/brugger1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/buildvisit/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/thirdparty_shared/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/2.13.1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/mpich/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/3.0.4/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/linux-x86_64_gcc-7.3/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/lib/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/libmpich.so', needed by '/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so'.  Stop.
make[1]: *** [CMakeFiles/Makefile2:68: CMakeFiles/EwdataDatabase_par.dir/all] Error 2
make: *** [Makefile:84: all] Error 2
...

you need to modify manually file:
    ./CMakeFiles/EwdataDatabase_par.dir/build.make
    
and remove from it following lines (they may be slightly different on your computer)

/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/net/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/wci-cl1.llnl.gov/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/vol/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/home/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/brugger1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/buildvisit/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/thirdparty_shared/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/2.13.1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/mpich/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/3.0.4/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/linux-x86_64_gcc-7.3/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/lib/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/libmpich.so
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/net/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/wci-cl1.llnl.gov/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/vol/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/home/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/brugger1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/buildvisit/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/thirdparty_shared/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/2.13.1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/mpich/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/3.0.4/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/linux-x86_64_gcc-7.3/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/lib/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/libopa.so
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/net/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/wci-cl1.llnl.gov/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/vol/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/home/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/brugger1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/buildvisit/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/thirdparty_shared/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/2.13.1/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/mpich/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/3.0.4/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/linux-x86_64_gcc-7.3/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/lib/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/libmpl.so
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/usr/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/lib64/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/librt.so
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/usr/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/lib64/opt/visit2_13_2.linux-x86_64/2.13.2/linux-x86_64/lib/libpthread.so

After the modification the file should look like:
...
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /usr/lib64/librt.so
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: /usr/lib64/libpthread.so
/home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so: CMakeFiles/EwdataDatabase_par.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gabrielw/MyProjects/dft/bsk/visit_plugin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX shared library /home/gabrielw/.visit/2.13.2/linux-x86_64/plugins/databases/libEwdataDatabase_par.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/EwdataDatabase_par.dir/link.txt --verbose=$(VERBOSE)
...

Similar operation you may need to apply for other `build.make` files. 
NOTE: if you execute commend
    xml2cmake -clobber wdata.xml
files `build.make` will be restored to original version. 
