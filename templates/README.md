## Content
* `legacy`: old template files
* `generic`: use this templates as staring point for you changes
* `simple`: templates that work with standard distributions of linux installed on desktop (like fedora)
* `dwarf`: http://datadwarf.if.pw.edu.pl/
* `tsubame`: https://helpdesk.t3.gsic.titech.ac.jp/manuals/handbook.en/
* `eagle`: https://wiki.man.poznan.pl/kdm/index.php?title=Eagle
* `daint`: https://www.cscs.ch/computers/piz-daint/
* `okeanos`: https://kdm.icm.edu.pl/Zasoby/komputery_w_icm.pl/#superkomputer-okeanos
* `lumi`: https://www.lumi-supercomputer.eu/
* `athena`: https://docs.cyfronet.pl/display/~plgpawlik/Athena
* `kamiak`: https://hpc.wsu.edu/kamiak-hpc/what-is-kamiak/
* `frontier`: https://www.olcf.ornl.gov/frontier/

## Installing templates
You can use tool
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
to copy selected templates to project folders
