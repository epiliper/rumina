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
