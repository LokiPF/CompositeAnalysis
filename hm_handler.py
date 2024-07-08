import subprocess


def run_tcl_script(hypermesh_path, tcl_script_path):
    command = f'"{hypermesh_path}" -tcl "{tcl_script_path}"'
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error running TCL script: {stderr.decode()}")
    else:
        print(f"TCL script output: {stdout.decode()}")


def run_solver(optistruct_path, solver_deck_path):
    command = f'"{optistruct_path}" "{solver_deck_path}"'
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error running solver: {stderr.decode()}")
    else:
        print(f"Solver output: {stdout.decode()}")


# Paths to executables, TCL scripts, and solver deck
hypermesh_path = 'C:\\Program Files\\Altair\\2023.1\\hwdesktop\\hm\\bin\\win64\\hmbatch.exe'
optistruct_path = 'C:\\Program Files\\Altair\\2023.1\\hwsolvers\\scripts\\optistruct.bat'
modify_tcl_script = './modify_young_modulus.tcl'
export_tcl_script = './export_solver_deck.tcl'
solver_deck_path = 'output.fem'

# Run the TCL script to modify the .hm file
run_tcl_script(hypermesh_path, modify_tcl_script)

# Run the TCL script to export the solver deck
run_tcl_script(hypermesh_path, export_tcl_script)

# Run the solver
run_solver(optistruct_path, solver_deck_path)
