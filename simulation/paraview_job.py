from paraview.simple import *
import glob
import os

vtk_dir = "_vtk"
files = sorted(glob.glob(os.path.join(vtk_dir, "*.vtr")))

if not files:
    raise RuntimeError("Nenhum ficheiro .vtr encontrado.")

# Abre os ficheiros como grupo
data = OpenDataFile(files)

# Mostra os dados (abre na vista ativa)
Show(data)

# Renderiza a vista para mostrar os dados
Render()

print("Dados abertos com sucesso!")

