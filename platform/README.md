## Info
This folder contains `Makefile` and the `machine.h` files dedicated to various machines.  
You are welcome to contribute to these templates. Once you have a template for a new machine, please send it to us or use the merge request functionality.

## Content
* `legacy`: old template files
* `generic`: use this template as a starting point for your changes
* `simple`: templates that work with standard distributions of Linux installed on desktops (like Fedora)
* `dwarf`: http://datadwarf.if.pw.edu.pl/
* `tsubame`: https://helpdesk.t3.gsic.titech.ac.jp/manuals/handbook.en/
* `eagle`: https://wiki.man.poznan.pl/kdm/index.php?title=Eagle
* `daint`: https://www.cscs.ch/computers/piz-daint/
* `okeanos`: https://kdm.icm.edu.pl/Zasoby/komputery_w_icm.pl/#superkomputer-okeanos
* `lumi`: https://www.lumi-supercomputer.eu/
* `athena`: https://docs.cyfronet.pl/display/~plgpawlik/Athena
* `kamiak`: https://hpc.wsu.edu/kamiak-hpc/what-is-kamiak/
* `frontier`: https://www.olcf.ornl.gov/frontier/
* `tsubame4`: https://www.t4.gsic.titech.ac.jp/en/hardware

## Installing templates
You can use the tool
```
./install-templates.py system
where system is:
        generic
        dwarf
        simple
        tsubame
        ...
        clean <-- use this to clear templates
```
to copy selected machine-dependent files to project folders.

## Note
In older versions, `VERSION<=2026.01.27`, this folder was named `templates`.  
