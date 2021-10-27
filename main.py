import argparse


def leeGrafoArchivo(file_path):
    archivo = open(file_path)
    lista = archivo.readlines()
    cantidadDeVertices = int(lista.pop(0))
    listaVertices = []
    for x in range(cantidadDeVertices):
        listaVertices.append(lista[x].rstrip("\n"))
    ARITAS = lista[cantidadDeVertices:]
    aristas = []
    for par in ARITAS:
        aristas.append((par.split()[0], par.split()[1]))
    return listaVertices, aristas

def main():
    # Definimos los argumentos de linea de comando que aceptamos
    parser = argparse.ArgumentParser()

    # Verbosidad, opcional, False por defecto
    parser.add_argument(
        '-v', '--verbose',
        action = 'store_true',
        help = 'Muestra mas informacion al correr el programa'
    )
    # Archivo del cual leer el grafo
    parser.add_argument(
        'file_name',
        help = "Archivo del cual leer el grafo a dibujar"
    )
    # Cantidad de iteraciones
    parser.add_argument(
        'iters',
        type = int,
        help = "Cantidad de iteraciones del algoritmo"
    )
    # Temperatura inicial
    parser.add_argument(
        '--temp',
        type = float,
        help = 'Temperatura inicial',
        default = 1000.0
    )
    # Cantidad de iteraciones entre actualizaciones de pantalla
    parser.add_argument(
        '--refresh',
        type = int,
        help = 'Cantidad de iteraciones entre actualizaciones de pantalla',
        default = 1
    )
    # Constante temperatura
    parser.add_argument(
        '--ctemp',
        type = float,
        help = 'Constante con la cual baja la temperatura cada step',
        default = 0.95
    )

    parser.add_argument(
        '--crepul',
        type = float,
        help = 'constante usada para calcular la repulsion entre vértices',
        default = 20
    )
    
    parser.add_argument(
        '--catrac',
        type = float,
        help = 'constante usada para calcular la atraccion de aristas',
        default = 3
    )
    
    parser.add_argument(
        '--cgrav',
        type = float,
        help = 'constante usada para calcular la fuerza de gravedad',
        default = 0.1
    )
    
    args = parser.parse_args()

    # Creamos nuestro objeto LayoutGraph
    layout_gr = LayoutGraph(
        leeGrafoArchivo(args.file_name),
        iters = args.iters,
        refresh = args.refresh,
        temperaturaInicial = args.temp,
        constanteTemperatura = args.ctemp,
        gravedad = args.cgrav,
        repultionConstant = args.crepul,
        attractionConstant = args.catrac,
        verbose = args.verbose
    )

    # Ejecutamos el layout
    layout_gr.layout()
    return


if __name__ == '__main__':
    main()
