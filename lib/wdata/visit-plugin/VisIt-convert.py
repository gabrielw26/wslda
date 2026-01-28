import os

def process_files(txt_files, new_dir="plugin-release"):

    base_find="wdata"
    base_repl="WData"

    list_find   =[
        "AVT_wdata",
        "avtwdata",
        "wdataVariable",
        "wdataRealVariable",
        "wdataComplexVariable",
        "wdataVectorVariable",
        " wdata ",
        '"wdata"',
        '"wdata_0.',
        "wdataPlugin",
        "wdataCommon",
        "wdataEngine",
        "_ENTRY(wdata,",
        "wdataMDServer",
        "wdataGeneral",
        "_VERSION(wdata,",
        "wdataDatabase"
        ]
    list_replace=[el.replace(base_find, base_repl) for el in list_find]
    print(list_replace)
    # Create the new directory if it doesn't exist
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    for file_path in txt_files:
        # Read the content of the file
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()

        for find, repl in zip(list_find,list_replace):
            # Find & Replace
            new_content = content.replace(find, repl)
            content=new_content

        # Get the base filename, Find & Replace
        base_name = os.path.basename(file_path)
        new_base_name = base_name.replace(base_find, base_repl)

        # Build the new file path
        new_file_path = os.path.join(new_dir, new_base_name)

        # Write the modified content to the new file
        with open(new_file_path, "w", encoding="utf-8") as f:
            f.write(new_content)
        print("%s --> %s" % (file_path,new_file_path))

# Example usage:
if __name__ == "__main__":
    # List your txt files here
    dest_dir="/home/gabrielw/MyProjects/torm/visit/src/databases/WData"
    txt_files = [
        "avtwdataFileFormat.C",
        "avtwdataFileFormat.h",
        # "wdata.code",
        "wdataCommonPluginInfo.C",
        "wdataEnginePluginInfo.C",
        "wdataMDServerPluginInfo.C",
        "wdataPluginInfo.C",
        "wdataPluginInfo.h"
        # "wdata.xml"
    ]
    process_files(txt_files, new_dir=dest_dir)

print("Update manually: wdata.xml")
print("Update manually [copy body while maintaining license header]: wdata.h")
print("Update manually [copy body while maintaining license header]: wdata.c")
