## Content
* `legacy`: old template files
* `generic`: use this templates as staring point for you changes
* `simple`: templates that work with standard distributions of linux installed on desktop (like fedora)
* `dwarf`: http://datadwarf.if.pw.edu.pl/
* `tsubame`: https://helpdesk.t3.gsic.titech.ac.jp/manuals/handbook.en/

## Installing templates
You can use tool
```
./install-templates.py system
where system is:
        generic
        dwarf
        simple
        tsubame
        clean <-- use this to clear templates
```
to copy selected templates to project folders
