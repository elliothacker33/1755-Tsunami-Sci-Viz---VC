from paraview.simple import *
import os

# ==========================================
# ⚙️ CONFIGURAÇÕES (Ajusta aqui se precisares)
# ==========================================
PVD_PATH = "_vtk/tsunami_1755.pvd"

# Exagero Vertical: Aumenta para veres a onda a subir paredes
WARP_SCALE = 1500.0

# Cores da Onda (Em metros)
MIN_WAVE = -2.0   # Cava da onda (Azul escuro)
MAX_WAVE = 10.0   # Crista da onda (Vermelho vivo)

def create_viz():
    # ----------------------------------------
    # 1. PREPARAR O AMBIENTE
    # ----------------------------------------
    # Apagar visualizações antigas para não sobrepor
    files = GetSources()
    for f in files.values():
        Delete(f)

    view = GetActiveViewOrCreate('RenderView')
    view.ViewSize = [1280, 720]
    view.OrientationAxesVisibility = 1

    # Fundo Gradiente "Cinemático"
    view.BackgroundColorMode = 'Gradient'
    view.Background = [0.15, 0.15, 0.25]  # Azul acinzentado escuro
    view.Background2 = [0.05, 0.05, 0.05] # Quase preto no topo

    if not os.path.exists(PVD_PATH):
        print(f"❌ ERRO: Ficheiro não encontrado: {PVD_PATH}")
        return

    print("🌊 A carregar dados...")
    data_source = PVDReader(FileName=PVD_PATH)
    data_source.UpdatePipeline()
"""
    # USAR MERGE BLOCKS (Mais robusto que CleanToGrid para GeoClaw)
    # Isto funde as várias caixas do AMR numa só malha
    merged = MergeBlocks(Input=data_source)

    # ----------------------------------------
    # 2. CAMADA DA TERRA (Cinzenta e Sólida)
    # ----------------------------------------
    print("🌍 A construir Portugal...")

    # Filtro: Apenas Terra (Batimetria > 0)
    thresh_terra = Threshold(Input=merged)
    thresh_terra.Scalars = ['POINTS', 'bathymetry']
    thresh_terra.LowerThreshold = 0.0
    thresh_terra.UpperThreshold = 99999.0

    # Dar volume (Exagero vertical)
    warp_terra = WarpByScalar(Input=thresh_terra)
    warp_terra.Scalars = ['POINTS', 'bathymetry']
    warp_terra.ScaleFactor = WARP_SCALE

    # Visualização
    disp_terra = Show(warp_terra, view)
    disp_terra.Representation = 'Surface'
    disp_terra.MapScalars = 0  # Desligar cores automáticas
    disp_terra.DiffuseColor = [0.6, 0.55, 0.5] # Cor de "Terra Seca" (Bege/Cinzento)
    disp_terra.Specular = 0.0  # Sem brilho (mate)

    # ----------------------------------------
    # 3. CAMADA DA ÁGUA (Colorida e Brilhante)
    # ----------------------------------------
    print("🌊 A encher o oceano...")

    # Filtro: Apenas Água (Profundidade > 0.1m)
    # Ignoramos zonas muito rasas para limpar o visual
    thresh_agua = Threshold(Input=merged)
    thresh_agua.Scalars = ['POINTS', 'water_depth']
    thresh_agua.LowerThreshold = 0.1
    thresh_agua.UpperThreshold = 99999.0

    # Dar volume (Atenção: Usamos surface_elevation para a água ficar por cima da terra)
    warp_agua = WarpByScalar(Input=thresh_agua)
    warp_agua.Scalars = ['POINTS', 'surface_elevation']
    warp_agua.ScaleFactor = WARP_SCALE

    # Visualização
    disp_agua = Show(warp_agua, view)

    # Pintar pela elevação (η)
    ColorBy(disp_agua, ('POINTS', 'surface_elevation'))

    # Configurar Tabela de Cores (Azul -> Branco -> Vermelho)
    lut = GetColorTransferFunction('surface_elevation')
    lut.ApplyPreset('Cool to Warm', True)

    # FIXAR O INTERVALO DE CORES (CRÍTICO!)
    # Isto garante que a água parada é branca/azul clara e o Tsunami é Vermelho
    lut.RescaleTransferFunction(MIN_WAVE, MAX_WAVE)

    # Estilo "Líquido"
    disp_agua.Opacity = 0.85      # Ligeira transparência
    disp_agua.Specular = 0.8      # Brilho forte
    disp_agua.SpecularPower = 80.0

    # ----------------------------------------
    # 4. FINALIZAR
    # ----------------------------------------
    # Barra de cores
    bar = GetScalarBar(lut, view)
    bar.Title = "Altura da Onda (m)"
    bar.ComponentTitle = ""
    bar.WindowLocation = 'LowerRightCorner'

    view.ResetCamera()
    GetActiveCamera().Elevation(20) # Levantar a câmara um pouco

    # Atualizar animação
    GetAnimationScene().UpdateAnimationUsingDataTimeSteps()

    Render()
    print("✅ Visualização Pronta! Carrega no PLAY.")
    """

if __name__ == "__main__":
    create_viz()
