import datetime
import getpass
import logging
import numpy as np
import os
import sys
import time

start_time = time.time()
date_time = datetime.datetime.now().strftime("%d.%m.%Y_%H.%M")
#configuration logger in console
#logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

#configuration logger in testsuite.log file

logFilePath = getpass.getuser() + '_logger_' + date_time
logging.basicConfig(filename=logFilePath, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

#disable all logger comunicates
#logging.disable(logging.CRITICAL)


listOfFolders, listOfTags, listOfMakes, listOfRuns, listOfDiffs, listOfChecks, listOfErrors = [], [], [], [], [], [], []

description_file = 'test.desc'
if len(sys.argv) == 2: tlist_file = str(sys.argv[1]) #Giventh value for tlist_file variable
else: tlist_file = 'tlist.txt' #Default variable tlist_file name is "tlist.txt"


# Setting default working directory
def get_global_dir():
    ggd = os.getcwd()
    return ggd

# Reading the text file with a list of directories with given tests to perform
def count_dirs():
    folder_list = []
    try:
        file = open(tlist_file, 'r').read().split('\n')
        logging.debug("File \"" + tlist_file + "\" opened.")
    except OSError:
        logging.critical("\"" + tlist_file + "\" file not found. End of a program.")
        exit()
    for line in file:
        if line.startswith('#') or len(line.strip()) == 0: continue
        else:
            folder_list.append(line)
            logging.debug("Directory \"" + str(line) + "\" found in " + "\"" + tlist_file + "\".")
    folder_list = [x for x in folder_list if x != '']
    logging.info("Directories found: " + str(len(folder_list)))
    return folder_list, len(folder_list)

# Checking if directory with given tests exists in working directory
def check(test_dir):
    if os.path.exists(test_dir):
        return True
    else:
        logging.warning("Directory \"" + test_dir + "\" not found!")
        return False

# Checking if description file exists
def check_desc(dd=description_file):
    if os.path.exists(dd):
        return True
    else:
        return False

# Counting the amount of tests to perform from each description file file
def count_tests():
    try:
        file = open(description_file, 'r').read().split('\n')
    except OSError:
        logging.error("\"" + description_file + "\" file not found!")
        return -1
    amnt = sum([1 for line in file if line.startswith('tag:')])
    return amnt

# Comparing output file to the reference file (value by value)
# If all differences of each values are below given tolerance, test is 'OK'
def run_diff_test(output, reference):
    out_vals, out_tags, ref_vals, ref_tags, ref_tols, comparedValues = [], [], [], [], [], []
    with open(output) as f1:
        for line1 in f1:
            if line1.startswith('#'): continue
            out_tag, out_value = line1.split()
            out_tags.append(str(out_tag))
            out_vals.append(float(out_value.rstrip('\n')))
    with open(reference) as f2:
        for line2 in f2:
            if line2.startswith('#'): continue
            ref_tag, ref_value, ref_tolerance = line2.split()
            ref_tags.append(str(ref_tag))
            ref_vals.append(float(ref_value.rstrip('\n')))
            ref_tols.append(float(ref_tolerance.rstrip('\n')))
    logging.info("Comparing output from file \"" + str(output) + "\" with reference values in file \"" + str(reference) + "\"")
    for i in range(len(out_vals)):
        if abs(out_vals[i] - ref_vals[i]) < ref_tols[i]:
            comparedValues.append(True)
            logging.info("PASS! - Output: " + str(out_tags[i]) + " = " + str(out_vals[i]) + " Reference: " + str(ref_tags[i]) + " = " + str(ref_vals[i]) + " Epsilon = " + str(ref_tols[i]) + " | (out_val - ref_val) = " + str(out_vals[i] - ref_vals[i]))
        else:
            logging.info("FAIL! - Output: " + str(out_tags[i]) + " = " + str(out_vals[i]) + " Reference: " + str(ref_tags[i]) + " = " + str(ref_vals[i]) + " Epsilon = " + str(ref_tols[i]) + " | (out_val - ref_val) = " + str(out_vals[i] - ref_vals[i]))
            comparedValues.append(False)
    if all(values is True for values in comparedValues):
        now_check = 'OK'
    else:
        now_check = 'FAIL'
    return now_check

# Executing Unix commands given for each test in description file
def commands(dir_name):
    tag_test_counter = 0
    tab_checks, tab_makes, tab_runs, tab_diffs = [], [], [], []
    file = open(description_file, "r").read().split('\n')
    for line in file:
        if line.startswith('#'): continue
        if line.startswith('tag:'):
            tagg, add_tag = line.split(': ')
            logging.info(">>>>>>>>>>>> Proceeding with test: " + dir_name + "->" + str(add_tag) + "<<<<<<<<<<<<<<<<<<<<<")
            print(">>>>>>>>>>>>>>>>>> TEST: " + dir_name + "->" + str(add_tag) + "<<<<<<<<<<<<<<<<<<<<<")
            sys.stdout.flush()
            listOfTags.append(add_tag)
            continue
        if line.startswith('exec:'):
            tagg, ex_comm = line.split(': ')
            logging.info("Executing command: " + str(ex_comm))
            os.system(ex_comm)
            continue
        if line.startswith('test:'):
            tag_test_counter += 1
            if tag_test_counter == 1:
                tagg, made_file = line.split(': ')
                logging.info("Checking if file \"" + str(made_file) + "\" exists.")
                if os.path.exists(made_file):
                    logging.info("File \"" + str(made_file) + "\" found.")
                    now_make = 'OK'
                else:
                    logging.error("FAIL - File \"" + str(made_file) + "\" not found.")
                    now_make = 'FAIL'
                continue
            if tag_test_counter == 2:
                tagg, out_file = line.split(': ')
                logging.info("Checking if file \"" + str(out_file) + "\" exists.")
                if os.path.exists(out_file):
                    logging.info("File \"" + str(out_file) + "\" found.")
                    now_run = 'OK'
                else:
                    logging.error(" FAIL - File \"" + str(out_file) + "\" not found.")
                    now_run = 'FAIL'
                tag_test_counter = 0
                continue
        if line.startswith('diff:'):
            tagg, diff_comm = line.split(': ')
            files_to_compare = diff_comm.split(' ')
            file1, file2 = files_to_compare[0], files_to_compare[1]
            logging.info("Checking if files \"" + str(file1) + "\" and \"" + str(file2) + "\" exist.")
            if os.path.exists(file1) and os.path.exists(file2):
                logging.info("Comparing files \"" + str(file1) + "\" and \"" + str(file2) + "\"")
                check_diff = run_diff_test(file1, file2)
                now_check = str(check_diff)
                tab_checks.append(now_check)
            else:
                logging.error("FAIL - File \"" + str(file1) + "\" or/and file \"" + str(file2) + "\" not found.")
                now_check = 'FAIL'
                tab_checks.append(now_check)
            tab_makes.append(str(now_make))
            tab_runs.append(str(now_run))
        if len(line) == 0: continue
    return tab_makes, tab_runs, tab_checks

# Generating final report
def runReport(__Folder, __Tag, __Make, __Run, __Check):
    #Settings of the report output file
    filepath = S + '/' + getpass.getuser() + '_report_' + date_time
    if os.path.exists(filepath): os.remove(filepath)
    max_arr = []
    tmp_arr = [__Folder, __Tag, __Make, __Run, __Check]
    for i in tmp_arr:
        max_arr.append(len(max(i, key=len)))
    spc = max(max_arr) + 2

    #Summarise total amount of tests, OK-tests and FAIL-tests.
    __tests = len(__Tag)
    __oks = 0
    __fails = 0
    for i in range(len(__Folder)):
        if __Make[i] == 'OK' and __Run[i] == 'OK' and __Check[i] == 'OK':
            __oks += 1
        else:
            __fails += 1
    __arr_title = ['--- Summary ---', '']
    __arr_test = ['TESTS:', __tests]
    __arr_ok = ['OK:', __oks]
    __arr_fail = ['FAIL:', __fails]
    summs = ([__arr_title, __arr_test, __arr_ok, __arr_fail])

    #Main part of report
    __Folder.insert(0, 'Folder')
    __Tag.insert(0, 'Tag')
    __Make.insert(0, 'Make')
    __Run.insert(0, 'Run')
    __Check.insert(0, 'Check')
    arrs = np.transpose([__Folder, __Tag, __Make, __Run, __Check])

    #Generating output file
    ff = open(filepath, 'a')
    np.savetxt(ff, summs, '%-7s', '\t')
    ff.write('\n')
    if __fails != 0: ff.write('Some test FAILED. For further information check your .log file.\n\n')
    else: ff.write('SUCCESS. No FAIL occured!\n\n')
    np.savetxt(ff, arrs, '%-{}s'.format(spc), '\t')
    ff.close()
    
    print("\n\nmore %s\n" % filepath)
    sys.stdout.flush()
    
    return __tests, __oks, __fails

# Cute print to the Terminal window
def cute_print():
    print('\n')
    printable_array = [listOfFolders, listOfTags, listOfMakes, listOfRuns, listOfChecks]
    col_width = max(len(word) for row in printable_array for word in row) + 2
    for row in printable_array:
        print("".join(word.ljust(col_width) for word in row))
    print('\n')
    sys.stdout.flush()

# Main function
S = get_global_dir()  # Set head directory
F, N = count_dirs()  # Amount and names of tests to perform
if N == 0:
    logging.critical("No directories found in \"" + tlist_file + "\". End of a program.")
    exit()
for folder in F:
    os.chdir(S)  # Return to the head directory after each iteration
    s = os.getcwd()  # Entering to the dir with tests inside to perform
    f = str(folder)
    check_if_dir_exists = check(f)  # Checking if dir with tests exists
    if check_if_dir_exists:
        os.chdir(s + '/{}'.format(f))
        check_desc_file = check_desc()  # Checking if description file exists
        if check_desc_file:
            logging.info("\"" + description_file + "\" opened in \"" + f + "\" directory.")
            amount = count_tests()  # Counting tags starting with "tag" from description file.
            # Returns the amount of tests in each dir
            logging.info("Number of tests in \"" + str(f) + "\\" + description_file + "\" : " + str(amount))
            if amount == 0:
                logging.error("There is no test in \"" + description_file + "\"! Skipping " + f + " directory.")
                continue
            if amount == -1:
                logging.error("There is no \"" + description_file + "\" file! Skipping " + f + " directory.")
                continue
            for test in range(amount):
                listOfFolders.append(f)
            makes, runs, checks = commands(f)  # EXECUTING TESTS
            for m in makes: listOfMakes.append(m)
            for r in runs: listOfRuns.append(r)
            for c in checks: listOfChecks.append(c)
        else:
            logging.warning("\"" + description_file + "\" not found in \"" + f + "\" directory. Skipping \"" + f + "\" directory.")
    else:
        logging.warning("\"" + f + "\" directory not found. Skipping \"" + f + "\" directory.")
        continue
logging.info("No more tests to do.")
logging.info("Creating the report.")
xx, yy, zz = runReport(listOfFolders, listOfTags, listOfMakes, listOfRuns, listOfChecks)

#cute_print()

sys.stdout.flush()
print('\n--- Summary ---')
print('TESTS: ' + str(xx))
print('OK: ' + str(yy))
print('FAIL: ' + str(zz))
logging.info("End of program.")
print('Done.')
print("--- %.8s seconds ---" % (time.time() - start_time))
sys.stdout.flush()
