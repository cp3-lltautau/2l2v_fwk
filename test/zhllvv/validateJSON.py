import json
import os

invalid_json_files = []
read_json_files = []

def parse():
    for files in os.listdir(os.getcwd()):
        if ((os.path.isdir(files) == False) and files.endswith('.json'))  :
            #print files
            with open(files) as json_file:
                try:
                    json.load(json_file)
                    read_json_files.append(files)
                except ValueError, e:
                    print ("file %s, JSON object issue: %s") % (files,e)
                    invalid_json_files.append(files)

    print "\n\n\n"
    print "The following files are not valid JSON: \n"
    listToDelete=""
    for ex in invalid_json_files:
        listToDelete+=" "
        listToDelete+=ex

    print listToDelete,"\n\n"
    print "there are ",len(read_json_files),"valid JSON files"

def main():
    # my code here
    parse()

if __name__ == "__main__":
    main()
