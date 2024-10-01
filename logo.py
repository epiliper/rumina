import shutil

Y = "\033[33m"
r = "\033[0m"
R = "\033[31m"
G = "\033[32m"
B = "\033[34m"
C = "\033[36m"
M = "\033[35m"

LOGO = f"""{B}
██████  ██    ██ ███    ███ ██ ███    ██  █████  
██   ██ ██    ██ ████  ████ ██ ████   ██ ██   ██ 
██████  ██    ██ ██ ████ ██ ██ ██ ██  ██ ███████ 
██   ██ ██    ██ ██  ██  ██ ██ ██  ██ ██ ██   ██ 
██   ██  ██████  ██      ██ ██ ██   ████ ██   ██ 
{r}
"""


def print_logo():
    term_width = shutil.get_terminal_size().columns
    logo_lines = LOGO.splitlines()

    for line in logo_lines:
        padding = (term_width - len(line)) // 2
        print(" " * padding + line)


def print_file_info(num_files, file_num, file_name):
    print(f"\n{C}FILE {file_num}/{num_files}:{r} {file_name.split("/")[-1]}")
    print(f"{C}============================={r}")


def print_file_end():
    print(f"{C}============================={r}")
