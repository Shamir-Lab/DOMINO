import os 
import random

def format_script(file_path, uid=True, **kwargs):
    formatted_script = open(file_path+".format").read().format(**kwargs)
    if uid:
        exec_file_name="{}_{}{}".format(os.path.splitext(file_path)[0] ,random.random(), os.path.splitext(file_path)[1]) #
    else:
        exec_file_name = file_path
    # int "exec_file: "+exec_file_name 
    open(exec_file_name, "w+").write(formatted_script)
    return exec_file_name


