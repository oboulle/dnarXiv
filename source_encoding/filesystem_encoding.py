#!/usr/bin/python3
import os
import argparse

"""
conversion of a file system (directories and files) into a single text file
"""


def directory_to_text(dir_path : str) -> str:
    """
    convert a directory and all its subdirectories/files into text 
    """
    text_dir =  "***DIR***"+dir_path+"***DIR***\n"
    
    for file in os.listdir(dir_path):
        file_path = os.path.join(dir_path, file)
        if os.path.isfile(file_path):
            text_dir += file_to_text(file_path)
        elif os.path.isdir(file_path):
            text_dir += directory_to_text(file_path)
    
    text_dir += "***ENDDIR***"+dir_path+"***ENDDIR***\n"
    return text_dir

    
def file_to_text(file_path: str) -> str:
    """
    convert a file into text
    """
    text_file = "***FILE***"+file_path+"***FILE***\n"
    file_input = open(file_path, 'r')
    try :
        line = file_input.readline()
    except UnicodeDecodeError:
        #print(file_path,"is not convertible")
        return ""
    
    while line != "":
        text_file += line
        line = file_input.readline()  
    file_input.close()
    
    # avoid error case when the file doesn't end with a newline
    if not text_file.endswith("\n"): text_file += "\n"
    
    text_file += "***ENDFILE***"+file_path+"***ENDFILE***\n"

    return text_file


def convert_filesystem_to_file(filesystem_path: str, output_path: str) -> str:
    """
    convert some filesystem located in a dir into a textfile
    """
    filesystem_abspath = os.path.abspath(filesystem_path)
    output_abspath = os.path.abspath(output_path)
    
    if os.path.isfile(output_abspath):
        print("error filesystem_encoder : ",output_abspath,"already exists")
        return
    
    prevdir = os.getcwd() # save current working dir
    os.chdir(os.path.expanduser(os.path.dirname(filesystem_abspath))) # set working dir to parent of the filesystem
    try:
        # turn the filesystem into text and save it
        filesystem_text = directory_to_text(os.path.basename(filesystem_abspath))
        file_output = open(output_abspath, "w")
        file_output.write(filesystem_text)
        file_output.close()
    finally:
        os.chdir(prevdir) # return to previous working dir
        
        
def text_to_filesystem(text_lines: list) -> None:
    """
    recursive function to extract all the documents from the text
    """
    
    # stop condition
    if text_lines == []:
        return
    
    # get the first document
    doc_type = text_lines[0].split("***")[1]
    doc_path = text_lines[0].split("***")[2]
    
    end_doc_index = text_lines.index("***END"+doc_type+"***"+doc_path+"***END"+doc_type+"***\n")
    
    # extract the lines associated to the content
    doc_content = text_lines[1:end_doc_index]
    
    if doc_type == "FILE":
        # write the file
        file_output = open(doc_path, "w")
        file_output.write("".join(doc_content))
        file_output.close()
        
    elif doc_type == "DIR":
        # create the dir and call recursively the function to its content
        if not os.path.isdir(doc_path): 
            os.mkdir(doc_path)
        text_to_filesystem(doc_content)
    
    # continue recursively with rest of the lines
    text_to_filesystem(text_lines[end_doc_index+1:])
    

def convert_file_to_filesystem(encoded_file_path: str, output_path: str) -> str:
    """
    convert an encoded textfile into a filesystem
    """
    file_abspath = os.path.abspath(encoded_file_path)
    output_abspath = os.path.abspath(output_path)
    
    if not os.path.isfile(file_abspath):
        print("error filesystem_encoder : ",file_abspath,"not found")
        return
    
    if not os.path.isdir(output_abspath):
        print("error filesystem_encoder : ",output_abspath,"not found")
        return
    
    prevdir = os.getcwd() # save current working dir
    os.chdir(os.path.expanduser(output_abspath)) # set working dir to parent of the filesystem
    try:
        # turn the text into a filesystem
        text_file = open(encoded_file_path)
        text_lines = text_file.readlines()
        dir_path = text_lines[0].split("***")[2]
        if not os.path.isdir(dir_path): os.mkdir(dir_path)
        text_to_filesystem(text_lines[1:-1])
    finally:
        os.chdir(prevdir) # return to previous working dir
    
    
# =================== main ======================= #
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='convert a file system into a single file')
    parser.add_argument('-i', action='store', dest='input_path', required=True,
                        help='path to the directory/file to convert')
    parser.add_argument('-o', action='store', dest='output_path', required=True,
                        help='path to the output single file/path to save the filesystem')

    # ---------- input list ---------------#
    arg = parser.parse_args()
    
    #convert_filesystem_to_file(arg.input_path, arg.output_path)

    convert_file_to_filesystem(arg.input_path, arg.output_path)
    
    # just import the functions
    