#!/usr/bin/env python3
import sys

def scan_for_define(file_name):
    """
    Scan the given file for #define statements in .h files.
    If #define statment is found:
        - Extract the name of the defined constant and its value. Keep values as strings.
        - Store the extracted information in a dictionary.
    Examples:
        #define NX (8) : Extract "NX" as the name and "(8)" as the value.
        #define MAX_USER_PARAMS 32 : Extract "MAX_USER_PARAMS" as the name and "32" as the value.
        #define SLDA_FORCE_A1 : Extract "SLDA_FORCE_A1" as the name and "" (empty string) as the value.
    Trim leading spaces and tabs from the extracted name and value. Ignore any comments that may be present on the same line as the #define statement.
    """
    defines = {}
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#define"):
                # Remove comments
                line = line.split('//')[0].split('/*')[0].strip()
                parts = line.split(maxsplit=2)
                if len(parts) >= 2:
                    name = parts[1]
                    value = parts[2] if len(parts) == 3 else ""
                    defines[name] = value
    return defines

def strip_from_value(value, chars):
    """
    Strip the specified characters from the value string.
    Example: If value is "(8)" and chars is "()", the result should be "8".
    """
    return value.strip(chars)

def code_dim(binary_name):
    """
    Determine the dimensionality of the code based on the binary name.
    If the binary name contains "1d", return "1D".
    If the binary name contains "2d", return "2D".
    If the binary name contains "3d", return "3".
    """
    binary_name = binary_name.lower()
    if "1d" in binary_name:
        return "1D"
    elif "2d" in binary_name:
        return "2D"
    elif "3d" in binary_name:
        return "3D"
    else:
        return None

def main(wslda_dir, binary_name):

    # replace in binary name "//" with "/" and remove leading and trailing spaces
    binary_name = binary_name.replace("//", "/").strip()
    code_dim_value = code_dim(binary_name)
    if code_dim_value is None:
        sys.exit(1) # wrong code, do not print anything, just exit with error code

    version_h = wslda_dir + "/VERSION.h"
    d_version_h = scan_for_define(version_h)

    predefines_h =  "./predefines.h"
    d_predefines_h = scan_for_define(predefines_h)

    # printing review report 
    key_width = 8
    key="BINARY"
    print("%8s : %s" % (key, binary_name))
    # enetries from VERSION.h
    for key, value in d_version_h.items():
        value = strip_from_value(value, '"()')
        print("%8s : %s" % (key, value))

    # extract lattice size
    NX=strip_from_value(d_predefines_h.get("NX", ""), '"()')
    NY=strip_from_value(d_predefines_h.get("NY", ""), '"()')
    NZ=strip_from_value(d_predefines_h.get("NZ", ""), '"()')
    key="LATTICE"
    print("%8s : %s x %s x %s" % (key, NX, NY, NZ))
    # extract spacings
    DX=strip_from_value(d_predefines_h.get("DX", ""), '"()')
    DY=strip_from_value(d_predefines_h.get("DY", ""), '"()')
    DZ=strip_from_value(d_predefines_h.get("DZ", ""), '"()')
    key="SPACING"
    print("%8s : %s x %s x %s" % (key, DX, DY, DZ))

    NX_int = int(NX) if NX.isdigit() else None
    NY_int = int(NY) if NY.isdigit() else None
    NZ_int = int(NZ) if NZ.isdigit() else None

    key="MODE"
    warn="Your computation needs caution!" 
    value="UNDEFINED"
    if code_dim_value=="3D" and NZ_int>1 and NY_int>1 and NX_int>1:
        value="FULL-3D"
        warn=""
    elif code_dim_value=="2D" and NZ_int>1 and NY_int>1 and NX_int>1:
        value="QUASI-2D"
        warn=""
    elif code_dim_value=="2D" and NZ_int==1 and NY_int>1 and NX_int>1:
        value="STRICT-2D"
        warn="WARNING: Check Wiki->Strict 2D or 1D mode for more details."
    elif code_dim_value=="1D" and NZ_int>1 and NY_int>1 and NX_int>1:
        value="QUASI-1D"
        warn=""
    elif code_dim_value=="1D" and NZ_int==1 and NY_int==1 and NX_int>1:
        value="STRICT-1D"
        warn="WARNING: Check Wiki->Strict 2D or 1D mode for more details."

    if warn!="":
        value+=" ["+warn+"]"
    print("%8s : %s" % (key, value))

    
# =========================================================================
# Main entry point
# =========================================================================
if __name__ == "__main__":
    # read from command line the first argument as the wslda directory
    if len(sys.argv) != 3:
        print("Usage: review-settings.py <wslda_dir> <binary_name>")
        sys.exit(1)
    wslda_dir = sys.argv[1]
    binary_name = sys.argv[2]

    main(wslda_dir, binary_name)
